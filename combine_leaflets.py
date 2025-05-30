import sys
import os
import argparse
from Bio import PDB 
import numpy as np
import string

class Lipid: 
    def __init__(self, name, id, max_z, min_z, lower_leaflet): 
        self.name = name
        self.id = id
        self.max_z = max_z
        self.min_z = min_z
        self.lower_leaflet = lower_leaflet

def is_bottom_leaflet(residue): 
    pass

def is_lipid(line, lipids):
    if line[0:4] == "ATOM": 
        residue_name = line[17:21].strip()
        if residue_name in lipids + ["DPPC", "DIPC", "CHOL", "POPC", "POPE", "POPG", "POPS", "DOPC", "DOPG", "DOPS", "DOPE", "DOPS", "DOPG", "DOPC", "DIPC", "DPPC", "CHOL"]:
            return True
    return False

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", help="Input pdb file")
    parser.add_argument("--output", help="Output pdb file")
    parser.add_argument("--lipids", help="List of lipid names to include (non standard)", nargs='+', default=[])
    parser.add_argument("--invert", help="If your system is inverted", action="store_true")

    args = parser.parse_args()

    lipid_lines = {}
    total_z = 0
    total_lipids = 0

    with open(args.pdb, 'r') as infile: 
        for line in infile: 
            if is_lipid(line, args.lipids):
                lipid_id = line[22:26].strip()
                if lipid_id not in lipid_lines.keys():
                    lipid_lines[lipid_id] = []
                    lipid_lines[lipid_id].append(line.strip())
                    total_lipids += 1
                    total_z += float(line[46:54])
                else:
                    lipid_lines[lipid_id].append(line.strip())
                    total_lipids += 1
                    total_z += float(line[46:54])

    avg_z = total_z / total_lipids

    lipids = {}
    max_z_lower = -np.inf
    min_z_upper = np.inf

    for lipid_id, lipid in lipid_lines.items(): 
        lipid_type = lipid[0][17:21].strip()
        z_max = -np.inf
        z_min = np.inf
        z_avg = 0
        for line in lipid: 
            z = float(line[46:54])
            if z > z_max:
                z_max = z
            if z < z_min:
                z_min = z
            z_avg += z

        z_avg /= len(lipid)

        if z_avg > avg_z: 
            lower_leaflet = False
            if z_min < min_z_upper:
                min_z_upper = z_min
        else:
            lower_leaflet = True
            if z_max > max_z_lower:
                max_z_lower = z_max

        if args.invert:
            lower_leaflet = not lower_leaflet

        lipids[lipid_id] = Lipid(lipid_type, lipid_id, z_max, z_min, lower_leaflet)

    z_gap = max_z_lower - min_z_upper + 1.5

    with open(args.pdb, 'r') as infile, open(args.output, 'w') as outfile: 
        for line in infile: 
            if is_lipid(line, args.lipids):
                if lipids[line[22:26].strip()].lower_leaflet:
                    line = line[:46] + str(float(line[46:54]) + z_gap) + line[54:]

            outfile.write(line)



if __name__ == "__main__":
    main(sys.argv[1:])
