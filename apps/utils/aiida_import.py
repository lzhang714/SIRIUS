import sys
import math
import json

au2angs = 1.889725989

def main():

    if len(sys.argv) != 3:
        print("Usage: python %s input.json length_units"%(sys.argv[0]))
        print("       where length_units = 'ang' or 'au'")
        sys.exit(0)

    jin = json.load(open(sys.argv[1], 'r'))

    lu = sys.argv[2]
    if lu not in ('ang', 'au'):
        print("wrong units of length")
        sys.exit(0)

    StructureData = DataFactory('structure')

    L = jin["unit_cell"]["lattice_vectors"]

    if lu == 'au':
        for i in range(3):
            for j in range(3):
                L[i][j] /= au2angs

    print("lattice vectors: ", L)

    s = StructureData(cell=L)

    for atom in jin["unit_cell"]["atoms"]:
        for pos in jin["unit_cell"]["atoms"][atom]:
            pos_cart = [0] * 3
            for x in range(3):
                pos_cart[x] = L[0][x] * pos[0] + L[1][x] * pos[1] + L[2][x] * pos[2]
            s.append_atom(position=pos_cart, symbols=atom)

    s.store()

    print "created structure with uuid='{}' and PK={}".format(s.uuid,s.pk)

    return

if __name__ == "__main__":
    main()

