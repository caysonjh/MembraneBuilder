import sys 

def main(args):
    if len(args) < 2:
        print("Usage: {} <filename>".format(args[0]))
        return 1

    filename = args[1]

    lipid_sizes = {'POPE': 125, 'POPG': 127, 'POPC': 134}

    with open("/Users/caysonhamilton/Biophysics/T4SS/INSANE_TESTING/charmm36.ff/lipid.rtp") as infile: 
        pope_line_count = lipid_sizes['POPE'] + 2
        popg_line_count = lipid_sizes['POPG'] + 2
        popc_line_count = lipid_sizes['POPC'] + 2
        correct_atoms = {'POPE': [], 'POPG': [], 'POPC': []}
        for line in infile: 
            if line.startswith("[ POPE ]") or pope_line_count < lipid_sizes['POPE']+2 and not pope_line_count < 0: 
                if pope_line_count >= lipid_sizes['POPE']: 
                    pass
                else: 
                    correct_atoms['POPE'].append(line.split()[0])
                pope_line_count -= 1
            elif line.startswith("[ POPG ]") or popg_line_count < lipid_sizes['POPG']+2 and not popg_line_count < 0:
                if popg_line_count >= lipid_sizes['POPG']: 
                    pass
                else: 
                    correct_atoms['POPG'].append(line.split()[0])
                popg_line_count -= 1
            elif line.startswith("[ POPC ]") or popc_line_count < lipid_sizes['POPC']+2 and not popc_line_count < 0:
                if popc_line_count >= lipid_sizes['POPC']: 
                    pass
                else: 
                    correct_atoms['POPC'].append(line.split()[0])
                popc_line_count -= 1
            else: 
                continue

    #[print(f'{lipid}: {atoms}\n {len(atoms)}') for lipid, atoms in correct_atoms.items()]

    with open(filename) as infile, open("{}_cleaned.pdb".format(filename.partition('.')[0]), "w") as outfile: 
        lipid_id_lines = {}
        for line in infile:
            if not line.startswith("HETATM"):
                outfile.write(line)
                continue
            else: 
                res_id = line[22:26].strip()
                lipid_type = line[17:21].strip()
                lipid_id = lipid_type + res_id
                if lipid_id not in lipid_id_lines:
                    lipid_id_lines[lipid_id] = []
                lipid_id_lines[lipid_id].append(line)

        for lipid_id, lines in lipid_id_lines.items():
            if len(lines) < lipid_sizes[lipid_id[:4]]:
                print(f'Too Short --> {lipid_id}: {len(lines)}')
                continue
            else: 
                all_correct = True
                for line in lines: 
                    lipid_type = line[17:21].strip()
                    if line[12:16].strip() not in correct_atoms[lipid_type]: 
                        print(f'Incorrect atom --> {lipid_type}: {line[12:16].strip()}')
                        all_correct = False
                        break
                    
                if all_correct:
                    outfile.writelines(lines)
                


def remove_incomplete_lipids(filename):
    main([None,filename])

if __name__ == "__main__":
    sys.exit(main(sys.argv))
