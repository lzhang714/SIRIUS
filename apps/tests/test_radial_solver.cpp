#include <sirius.hpp>
#include "radial/radial_solver.hpp"

using namespace sirius;

void test_radial_solver()
{
    Radial_grid_lin_exp<double> rgrid(1500, 1e-7, 2.0);
    int zn{38};
    std::vector<double> v(rgrid.num_points());
    for (int ir = 0; ir < rgrid.num_points(); ir++) {
        v[ir] = -zn * rgrid.x_inv(ir);
    }

    Radial_solver rsolver(zn, v, rgrid);
    std::vector<double> p, rdudr;
    std::array<double, 2> uderiv;
    rsolver.solve(relativity_t::iora, 1, 2, -0.524233, p, rdudr, uderiv);

    std::stringstream s;
    s << "radial_functions.dat";
    FILE* fout = fopen(s.str().c_str(), "w");

    for (int ir = 0; ir < rgrid.num_points(); ir++) {
        fprintf(fout, "%f ", rgrid[ir]);
        fprintf(fout, "%f ", p[ir]);
        fprintf(fout, "%f ", rdudr[ir]);
        fprintf(fout, "\n");
    }
    fclose(fout);
}

void test_radial_solver_v2()
{
    auto rgrid = Radial_grid_factory<double>(radial_grid_t::power, 10000, 1e-7, 222.0, 2);

    int zn{11};
    std::vector<double> v(rgrid.num_points());
    for (int ir = 0; ir < rgrid.num_points(); ir++) {
        v[ir] = -zn * rgrid.x_inv(ir);
    }

    Radial_solver rsolver(zn, v, rgrid);
    std::vector<double> p(rgrid.num_points());
    std::vector<double> dpdr(rgrid.num_points());
    Spline<double> chip(rgrid);
    std::vector<double> q(rgrid.num_points());
    std::vector<double> dqdr(rgrid.num_points());
    Spline<double> chiq(rgrid);

    double enu{-2.42};
    double de{0.1};

    int l{4};
    int n{5};
    //rsolver.integrate_forward_rk4<relativity_t::none, true>(-0.125, l, 0, chip, chiq, p, dpdr, q, dqdr);
    int nn = rsolver.integrate_forward_gsl(enu, l, 0, chip, chiq, p, dpdr, q, dqdr, true);
    std::cout << "enu="<<enu<<", nn="<<nn<<std::endl;

    //double emin, emax;
    //if (nn > n - l - 1) {
    //    emax = enu;
    //    while (nn > n - l - 1)  {
    //        enu -= de;
    //        de *= 1.3;
    //        nn = rsolver.integrate_forward_gsl(enu, l, 0, chip, chiq, p, dpdr, q, dqdr, true);
    //        std::cout << "enu="<<enu<<", nn="<<nn<<std::endl;
    //    }
    //    emin = enu;
    //    std::cout << "#1 : emin, emax, enu = " << emin <<", " << emax << ", " << enu << std::endl;
    //}
    //if (nn <= n - l - 1) {
    //    emin = enu;
    //    while (nn <= n - l - 1)  {
    //        enu += de;
    //        de *= 1.3;
    //        nn = rsolver.integrate_forward_gsl(enu, l, 0, chip, chiq, p, dpdr, q, dqdr, true);
    //    }
    //    emax = enu;
    //    std::cout << "#2 : emin, emax, enu = " << emin <<", " << emax << ", " << enu << std::endl;
    //}

    //std::cout << "emin, emax = " << emin <<", "<<emax << std::endl;

    //while (true) {
    //    enu = (emin + emax) / 2.0;
    //    nn = rsolver.integrate_forward_gsl(enu, l, 0, chip, chiq, p, dpdr, q, dqdr, true);
    //    if (nn <= n - l - 1) {
    //        emin = enu;
    //    } else {
    //        emax = enu;
    //    }
    //    std::cout << "new emin=" << emin << std::endl;
    //    std::cout << "new emax=" << emax << std::endl;
    //    std::cout << "nn="<<nn<<", std::abs(emin - emax)=" <<std::abs(emin - emax)<< std::endl;
    //    if (std::abs(emin - emax) < 1e-10 && nn == n - l - 1) {
    //        break;
    //    }
    //}
    //std::cout << "emin, emax = " << emin <<", "<<emax << std::endl;




    //Spline<double> sp(rgrid);
    //for (int i = 0; i < rgrid.num_points(); i++) {
    //    sp(i) = p[i] * p[i];
    //}
    //double norm = sp.interpolate().integrate(0);
    //norm = std::pow(norm, -0.5);
    //for (int i = 0; i < rgrid.num_points(); i++) {
    //    p[i] *= norm;
    //    dpdr[i] *= norm;
    //}
    //

    Bound_state bs(relativity_t::none, zn, n, l, 0, rgrid, v, enu);


    std::stringstream s;
    s << "radial_functions.dat";
    FILE* fout = fopen(s.str().c_str(), "w");

    for (int ir = 0; ir < rgrid.num_points(); ir++) {
        fprintf(fout, "%18.12f %18.12f %18.12f\n", rgrid[ir], p[ir], dpdr[ir]);
    }
    fclose(fout);
}

int main(int argn, char** argv)
{
    cmd_args args;

    args.parse_args(argn, argv);
    if (args.exist("help"))
    {
        printf("Usage: %s [options]\n", argv[0]);
        args.print_help();
        return 0;
    }

    sirius::initialize(1);
    test_radial_solver_v2();
    sirius::finalize();
}
