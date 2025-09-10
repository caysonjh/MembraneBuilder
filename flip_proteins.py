import sys

pdb_file = sys.argv[1]

with open(pdb_file) as f, open(pdb_file[:-4] + "_flipped.pdb", "w") as g:
    for line in f:
        if line.startswith("ATOM"):
            new_z = -float(line[46:54])
            new_line = line[:46] + "{:8.3f}".format(new_z) + line[54:]
            g.write(new_line)
        else:
            g.write(line)