#include "utils/profiler.hpp"
#include <sirius.hpp>
#include "filesystem.hpp"
#include <utils/json.hpp>
#include <cfenv>
#include <fenv.h>

using namespace sirius;
using json = nlohmann::json;

const std::string aiida_output_file = "output_aiida.json";

enum class task_t : int
{
    ground_state_new     = 0,
    ground_state_restart = 1,
    k_point_path         = 2
};

void json_output_common(json& dict__)
{
    dict__["git_hash"] = sirius::git_hash();
    dict__["comm_world_size"] = Communicator::world().size();
    dict__["threads_per_rank"] = omp_get_max_threads();
}

void rewrite_relative_paths(json& dict__, fs::path const &working_directory = fs::current_path())
{
    // the json.unit_cell.atom_files[] dict might contain relative paths,
    // which should be relative to the json file. So better make them
    // absolute such that the simulation context does not have to be
    // aware of paths.
    if (!dict__.count("unit_cell"))
        return;

    auto &section = dict__["unit_cell"];

    if (!section.count("atom_files"))
        return;

    auto &atom_files = section["atom_files"];

    for (auto& label : atom_files.items()) {
        label.value() = working_directory / std::string(label.value());
    }
}

nlohmann::json preprocess_json_input(std::string fname__)
{
    if (fname__.find("{") == std::string::npos) {
        // If it's a file, set the working directory to that file.
        auto json = utils::read_json_from_file(fname__);
        rewrite_relative_paths(json, fs::path{fname__}.parent_path());
        return json;
    } else {
        // Raw JSON input
        auto json = utils::read_json_from_string(fname__);
        rewrite_relative_paths(json);
        return json;
    }
}

std::unique_ptr<Simulation_context> create_sim_ctx(std::string fname__,
                                                   cmd_args const& args__)
{

    auto json = preprocess_json_input(fname__);

    auto ctx_ptr = std::make_unique<Simulation_context>(json.dump(), Communicator::world());
    Simulation_context& ctx = *ctx_ptr;

    auto& inp = ctx.cfg().parameters();
    if (inp.gamma_point() && !(inp.ngridk()[0] * inp.ngridk()[1] * inp.ngridk()[2] == 1)) {
        TERMINATE("this is not a Gamma-point calculation")
    }

    ctx.import(args__);

    return ctx_ptr;
}


double ground_state(Simulation_context& ctx,
                    task_t              task,
                    cmd_args const&     args,
                    int                 write_output)
{
    ctx.print_memory_usage(__FILE__, __LINE__);

    auto& inp = ctx.cfg().parameters();

    std::string ref_file = args.value<std::string>("test_against", "");
    /* don't write output if we compare against the reference calculation */
    bool write_state = (ref_file.size() == 0);

    bool const reduce_kp = ctx.use_symmetry() && ctx.cfg().parameters().use_ibz();
    K_point_set kset(ctx, ctx.cfg().parameters().ngridk(), ctx.cfg().parameters().shiftk(), reduce_kp);
    DFT_ground_state dft(kset);

    ctx.print_memory_usage(__FILE__, __LINE__);

    auto& potential = dft.potential();
    auto& density = dft.density();

    if (task == task_t::ground_state_restart) {
        if (!utils::file_exists(storage_file_name)) {
            TERMINATE("storage file is not found");
        }
        density.load();
        potential.load();
    } else {
        dft.initial_state();
    }

    /* launch the calculation */
    auto result = dft.find(inp.density_tol(), inp.energy_tol(), ctx.cfg().iterative_solver().energy_tolerance(),
            inp.num_dft_iter(), write_state);
    /* compute forces and stress */
    if (ctx.cfg().control().print_stress() && !ctx.full_potential()) {
        Stress& s       = dft.stress();
        auto stress_tot = s.calc_stress_total();
        s.print_info();
        result["stress"] = std::vector<std::vector<double>>(3, std::vector<double>(3));
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result["stress"][i][j] = stress_tot(j, i);
            }
        }
    }
    if (ctx.cfg().control().print_forces()) {
        Force& f         = dft.forces();
        auto& forces_tot = f.calc_forces_total();
        f.print_info();
        result["forces"] = std::vector<std::vector<double>>(ctx.unit_cell().num_atoms(), std::vector<double>(3));
        for (int i = 0; i < ctx.unit_cell().num_atoms(); i++) {
            for (int j = 0; j < 3; j++) {
                result["forces"][i][j] = forces_tot(j, i);
            }
        }
    }

    if (write_state && write_output) {
        json dict;
        json_output_common(dict);

        dict["task"] = static_cast<int>(task);
        dict["context"] = ctx.serialize();
        dict["ground_state"] = result;
        //dict["timers"] = utils::timer::serialize();
        dict["counters"] = json::object();
        dict["counters"]["local_operator_num_applied"] = ctx.num_loc_op_applied();
        dict["counters"]["band_evp_work_count"] = ctx.evp_work_count();

        if (ctx.comm().rank() == 0) {
            std::string output_file = args.value<std::string>("output", std::string("output_") +
                                                              ctx.start_time_tag() + std::string(".json"));
            std::ofstream ofs(output_file, std::ofstream::out | std::ofstream::trunc);
            ofs << dict.dump(4);
        }

        //if (args.exist("aiida_output")) {
        //    json dict;
        //    json_output_common(dict);
        //    dict["task"] = static_cast<int>(task);
        //    if (result >= 0) {
        //        dict["task_status"] = "converged";
        //        dict["num_scf_iterations"] =  result;
        //    } else {
        //        dict["task_status"] = "unconverged";
        //    }
        //    dict["volume"] = ctx.unit_cell().omega() * std::pow(bohr_radius, 3);
        //    dict["volume_units"] = "angstrom^3";
        //    dict["energy"] = dft.total_energy() * ha2ev;
        //    dict["energy_units"] = "eV";
        //    if (ctx.comm().rank() == 0) {
        //        std::ofstream ofs(aiida_output_file, std::ofstream::out | std::ofstream::trunc);
        //        ofs << dict.dump(4);
        //    }
        //}
    }


    if (ctx.cfg().control().verification() >= 1) {
        dft.check_scf_density();
    }

    auto repeat_update = args.value<int>("repeat_update", 0);
    if (repeat_update) {
        auto lv = ctx.unit_cell().lattice_vectors();
        auto a = std::pow(ctx.unit_cell().omega(), 1.0 / 3);
        for (int i = 0; i < repeat_update; i++) {
            double t = static_cast<double>(i) / repeat_update;
            auto lv1 = lv;
            for (int x: {0, 1, 2}) {
                lv1(x, 0) = lv(x, 0) + 0.15 * a * std::sin(t * twopi);
                lv1(x, 1) = lv(x, 1) + 0.15 * a * std::cos(t * twopi);
            }
            ctx.unit_cell().set_lattice_vectors(lv1);
            dft.update();
            auto r1 = dft.find(inp.density_tol(), inp.energy_tol(), ctx.cfg().iterative_solver().energy_tolerance(),
                    inp.num_dft_iter(), write_state);
            if (ctx.cfg().control().print_stress() && !ctx.full_potential()) {
                Stress& s       = dft.stress();
                auto stress_tot = s.calc_stress_total();
                auto elem = std::vector<std::vector<double>>(3, std::vector<double>(3));
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        elem[i][j] = stress_tot(j, i);
                    }
                }
                r1["stress"] = elem;
            }
            if (ctx.cfg().control().print_forces()) {
                Force& f         = dft.forces();
                auto& forces_tot = f.calc_forces_total();
                auto elem = std::vector<std::vector<double>>(ctx.unit_cell().num_atoms(), std::vector<double>(3));
                for (int i = 0; i < ctx.unit_cell().num_atoms(); i++) {
                    for (int j = 0; j < 3; j++) {
                        elem[i][j] = forces_tot(j, i);
                    }
                }
                r1["forces"] = elem;
            }
        }
    }

    //dft.print_magnetic_moment();

    if (ref_file.size() != 0) {
        json dict_ref;
        std::ifstream(ref_file) >> dict_ref;

        double e1 = result["energy"]["total"].get<double>();
        double e2 = dict_ref["ground_state"]["energy"]["total"].get<double>();

        if (std::abs(e1 - e2) > 1e-5) {
            std::printf("total energy is different: %18.7f computed vs. %18.7f reference\n", e1, e2);
            ctx.comm().abort(1);
        }
        if (result.count("stress") && dict_ref["ground_state"].count("stress")) {
            double diff{0};
            auto s1 = result["stress"].get<std::vector<std::vector<double>>>();
            auto s2 = dict_ref["ground_state"]["stress"].get<std::vector<std::vector<double>>>();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    diff += std::abs(s1[i][j] - s2[i][j]);
                }
            }
            if (diff > 1e-5) {
                std::printf("total stress is different!");
                std::cout << "  reference: " << dict_ref["ground_state"]["stress"] << "\n";
                std::cout << "  computed: " << result["stress"] << "\n";
                ctx.comm().abort(2);
            }
        }
        if (result.count("forces") && dict_ref["ground_state"].count("forces")) {
            double diff{0};
            auto s1 = result["forces"].get<std::vector<std::vector<double>>>();
            auto s2 = dict_ref["ground_state"]["forces"].get<std::vector<std::vector<double>>>();
            for (int i = 0; i < ctx.unit_cell().num_atoms(); i++) {
                for (int j = 0; j < 3; j++) {
                    diff += std::abs(s1[i][j] - s2[i][j]);
                }
            }
            if (diff > 1e-6) {
                std::printf("total force is different!");
                std::cout << "  reference: " << dict_ref["ground_state"]["forces"] << "\n";
                std::cout << "  computed: " << result["forces"] << "\n";
                ctx.comm().abort(3);
            }
        }
    }

    /* wait for all */
    ctx.comm().barrier();

    return result["energy"]["total"].get<double>();
}

/// Run a task based on a command line input.
void run_tasks(cmd_args const& args)
{
    /* get the task id */
    task_t task = static_cast<task_t>(args.value<int>("task", 0));

    /* get the input file name */
    auto fpath = args.value<fs::path>("input", "sirius.json");

    if (fs::is_directory(fpath)) {
        fpath /= "sirius.json";
    }

    if (!fs::exists(fpath)) {
        if (Communicator::world().rank() == 0) {
            std::printf("input file does not exist\n");
        }
        return;
    }

    auto fname = fpath.string();

    if (task == task_t::ground_state_new || task == task_t::ground_state_restart) {
        auto ctx = create_sim_ctx(fname, args);
        ctx->initialize();
        if (ctx->comm().rank() == 0) {
            auto dict = ctx->serialize();
            std::ofstream ofs("setup.json", std::ofstream::out | std::ofstream::trunc);
            ofs << dict.dump(4);
        }
        //if (ctx->full_potential()) {
        //    ctx->gk_cutoff(ctx->aw_cutoff() / ctx->unit_cell().min_mt_radius());
        //}
        ground_state(*ctx, task, args, 1);
    }

    if (task == task_t::k_point_path) {
        auto ctx = create_sim_ctx(fname, args);
        ctx->cfg().iterative_solver().energy_tolerance(1e-12);
        ctx->gamma_point(false);
        ctx->initialize();
        //if (ctx->full_potential()) {
        //    ctx->gk_cutoff(ctx->aw_cutoff() / ctx->unit_cell().min_mt_radius());
        //}

        Potential potential(*ctx);

        Density density(*ctx);

        K_point_set ks(*ctx);

        json inp;
        std::ifstream(fname) >> inp;

        /* list of pairs (label, k-point vector) */
        std::vector<std::pair<std::string, std::vector<double>>> vertex;

        auto labels = inp["kpoints_path"].get<std::vector<std::string>>();
        for (auto e: labels) {
            auto v = inp["kpoints_rel"][e].get<std::vector<double>>();
            vertex.push_back({e, v});
        }

        std::vector<double> x_axis;
        std::vector<std::pair<double, std::string>> x_ticks;

        /* first point */
        x_axis.push_back(0);
        x_ticks.push_back({0, vertex[0].first});
        ks.add_kpoint(&vertex[0].second[0], 1.0);

        double t{0};
        for (size_t i = 0; i < vertex.size() - 1; i++) {
            vector3d<double> v0 = vector3d<double>(vertex[i].second);
            vector3d<double> v1 = vector3d<double>(vertex[i + 1].second);
            vector3d<double> dv = v1 - v0;
            vector3d<double> dv_cart = dot(ctx->unit_cell().reciprocal_lattice_vectors(), dv);
            int np = std::max(10, static_cast<int>(30 * dv_cart.length()));
            for (int j = 1; j <= np; j++) {
                vector3d<double> v = v0 + dv * static_cast<double>(j) / np;
                ks.add_kpoint(&v[0], 1.0);
                t += dv_cart.length() / np;
                x_axis.push_back(t);
            }
            x_ticks.push_back({t, vertex[i + 1].first});
        }

        ks.initialize();

        //density.initial_density();
        density.load();
        potential.generate(density, ctx->use_symmetry(), true);
        Band band(*ctx);
        Hamiltonian0<double> H0(potential);
        if (!ctx->full_potential()) {
            band.initialize_subspace(ks, H0);
            if (ctx->hubbard_correction()) {
                TERMINATE("fix me");
                //potential.U().compute_occupation_matrix(ks); // TODO: this is wrong; U matrix should come form the saved file
                //potential.U().calculate_hubbard_potential_and_energy(potential.U().occupation_matrix());
            }
        }
        band.solve<double, double>(ks, H0, true, ctx->cfg().iterative_solver().energy_tolerance());

        ks.sync_band<double, sync_band_t::energy>();
        if (Communicator::world().rank() == 0) {
            json dict;
            dict["header"] = {};
            dict["header"]["x_axis"] = x_axis;
            dict["header"]["x_ticks"] = std::vector<json>();
            dict["header"]["num_bands"] = ctx->num_bands();
            dict["header"]["num_mag_dims"] = ctx->num_mag_dims();
            for (auto& e: x_ticks) {
                json j;
                j["x"] = e.first;
                j["label"] = e.second;
                dict["header"]["x_ticks"].push_back(j);
            }
            dict["bands"] = std::vector<json>();

            for (int ik = 0; ik < ks.num_kpoints(); ik++) {
                json bnd_k;
                bnd_k["kpoint"] = std::vector<double>(3, 0);
                for (int x = 0; x < 3; x++) {
                    bnd_k["kpoint"][x] = ks.get<double>(ik)->vk()[x];
                }
                std::vector<double> bnd_e;

                for (int ispn = 0; ispn < ctx->num_spinors(); ispn++) {
                    for (int j = 0; j < ctx->num_bands(); j++) {
                        bnd_e.push_back(ks.get<double>(ik)->band_energy(j, ispn));
                    }
                }
                //ks.get_band_energies(ik, bnd_e.data());
                bnd_k["values"] = bnd_e;
                dict["bands"].push_back(bnd_k);
            }
            std::ofstream ofs("bands.json", std::ofstream::out | std::ofstream::trunc);
            ofs << dict.dump(4);
        }
    }
}

int main(int argn, char** argv)
{
    std::feclearexcept(FE_ALL_EXCEPT);
    cmd_args args;
    args.register_key("--input=", "{string} input file name");
    args.register_key("--output=", "{string} output file name");
    args.register_key("--task=", "{int} task id");
    args.register_key("--aiida_output", "write output for AiiDA");
    args.register_key("--test_against=", "{string} json file with reference values");
    args.register_key("--repeat_update=", "{int} number of times to repeat update()");
    args.register_key("--fpe", "enable check of floating-point exceptions using GNUC library");
    args.register_key("--control.processing_unit=", "");
    args.register_key("--control.verbosity=", "");
    args.register_key("--control.verification=", "");
    args.register_key("--control.mpi_grid_dims=","");
    args.register_key("--control.std_evp_solver_name=", "");
    args.register_key("--control.gen_evp_solver_name=", "");
    args.register_key("--control.fft_mode=", "");
    args.register_key("--control.memory_usage=", "");
    args.register_key("--parameters.ngridk=", "");
    args.register_key("--parameters.gamma_point=", "");
    args.register_key("--parameters.pw_cutoff=", "");
    args.register_key("--iterative_solver.orthogonalize=", "");
    args.register_key("--iterative_solver.early_restart=", "{double} value between 0 and 1 to control the early restart ratio in Davidson");
    args.register_key("--mixer.type=", "{string} mixer name (anderson, anderson_stable, broyden2, linear)");
    args.register_key("--mixer.beta=", "{double} mixing parameter");

    args.parse_args(argn, argv);
    if (args.exist("help")) {
        std::printf("Usage: %s [options]\n", argv[0]);
        args.print_help();
        return 0;
    }

#if defined(_GNU_SOURCE)
    if (args.exist("fpe")) {
        feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
    }
#endif

    sirius::initialize(1);

    run_tasks(args);

    int my_rank = Communicator::world().rank();

    sirius::finalize(1);

    if (my_rank == 0)  {
        //auto timing_result = ::utils::global_rtgraph_timer.process().flatten(1).sort_nodes();
        auto timing_result = ::utils::global_rtgraph_timer.process();
        std::cout << timing_result.print({rt_graph::Stat::Count, rt_graph::Stat::Total, rt_graph::Stat::Percentage,
                                          rt_graph::Stat::SelfPercentage, rt_graph::Stat::Median, rt_graph::Stat::Min,
                                          rt_graph::Stat::Max});
        std::ofstream ofs("timers.json", std::ofstream::out | std::ofstream::trunc);
        ofs << timing_result.json();
    }
    if (std::fetestexcept(FE_DIVBYZERO)) {
        std::cout << "FE_DIVBYZERO exception\n";
    }
    if (std::fetestexcept(FE_INVALID)) {
        std::cout << "FE_INVALID exception\n";
    }
    if (std::fetestexcept(FE_UNDERFLOW)) {
        std::cout << "FE_UNDERFLOW exception\n";
    }
    if (std::fetestexcept(FE_OVERFLOW)) {
        std::cout << "FE_OVERFLOW exception\n";
    }

    return 0;
}
