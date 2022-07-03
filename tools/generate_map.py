# python 3 script to generate the map file for fix-phonon
#
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.lammps.data import LammpsData
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if __name__ == "__main__":
    """
    Script to generate the map file needed by fix-phonon from a data file or POSCAR.

    Pymatgen (https://pymatgen.org) is required, and the data file or POSCAR should correspond
    to the ideal crystal structure for the code to work correctly. The generated map file might
    differ from the one by latgen (https://github.com/lingtikong/latgen) but should work as well.

    For help, please use:
        python generate_map -h

    (C) LT Kong (konglt@sjtu.edu.cn)
    July 2022

    """
    import os
    import numpy as np
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="The structure file to analyze, must be either a vasp POSCAR or a LAMMPS data file for point particle.")
    parser.add_argument("-s", "--style", default="poscar", choices=['poscar', 'atomic', 'charge', 'bond', 'angle', 'full', 'molecular'],
                        help="To define style of the input structure file; by default: poscar.")
    parser.add_argument("-o", "--output", default="map.in", help="To define the output file name; by default: map.in.")
    parser.add_argument("-p", "--symprec", type=float, default=0.01, help="To define the precision for spacegroup analysis; by default: 0.01.")
    parser.add_argument("-a", "--angtol", type=float, default=5, help="To define the angle tolerance for spacegroup analysis; by default: 5.")
    args = parser.parse_args()
    
    if not os.path.exists(args.file):
        print(f"\nFile {args.file} required but not found!\n")
   
    elif args.style.lower() == 'poscar':
        poscar = Poscar.from_file(args.file)
        structure = poscar.structure

    else:
        try:
            lmp = LammpsData.from_file(args.file, atom_style = args.style, sort_id = True)
            structure = lmp.structure

        except:
            raise(f"Failed to read {args.file} with {args.style} format.")

    sa = SpacegroupAnalyzer(structure, symprec= args.symprec, angle_tolerance= args.angtol)
    primitive = sa.find_primitive()

    print(f"\nThe lattice parameters of the supercell are:\n  {structure.lattice.parameters}")
    print(f"\nThe lattice parameters of the primitive cell are:\n  {primitive.lattice.parameters}")

    ang_diff = np.array(structure.lattice.angles) - np.array(primitive.lattice.angles)
    if np.sqrt(np.sum(ang_diff * ang_diff)) > args.angtol:
       print("\nWARNING: It seems that the supercell and the identified primitive cell differ\ntoo much, the result might not be reliable!\n")

    nx = int((structure.lattice.a + args.symprec)/primitive.lattice.a)
    ny = int((structure.lattice.b + args.symprec)/primitive.lattice.b)
    nz = int((structure.lattice.c + args.symprec)/primitive.lattice.c)
    if nx < 1 or ny < 1 or nz < 1:
       print(f"Your supercell does not seem to be a simple repetition of the primitive one, the mission cannot be proceeded!\n")
       exit()

    nucell = primitive.num_sites
    if nx*ny*nz*nucell != structure.num_sites:
       print(f"Your supercell does not seem to be a simple repetition of the primitive one, the mission cannot be proceeded!\n")
       exit()

    print(f"\nYour supercell contains {nx} x {ny} x {nz} primitive cells, each of {nucell} atoms.")

    map2prim = sa.get_symmetry_dataset()['mapping_to_primitive']
    id_prim = dict()
    for i,imap in enumerate(map2prim):
        if imap == 0:
           origin = i
           break

    id_prim[0] = origin
    dist2orig = [-1]*nucell
    for i,site in enumerate(structure.sites):
        dist = structure.get_distance(origin, i)
        iu = map2prim[i]
        if dist2orig[iu] < 0 or dist < dist2orig[iu]:
           dist2orig[iu] = dist
           id_prim[iu] = i

    map_string = list()
    map_string.append(f"{nx} {ny} {nz} {nucell}")
    map_string.append(f"Map file for the supercell in {args.file}")

    for i,site in enumerate(structure.sites):
        iu = map2prim[i]
        ip = id_prim[iu]
        rel_coords = np.mod(site.frac_coords - structure.sites[ip].frac_coords, 1.)
        ix, iy, iz = map(int, np.mod(rel_coords*[nx, ny, nz] + args.symprec, [nx, ny, nz]))
        map_string.append(f"{ix} {iy} {iz} {iu} {i+1}")

    fp = open(args.output, "w")
    fp.write("\n".join(map_string))
    fp.close()
    print(f"\nMission completed, the map file is written to {args.output}.\n\n")
