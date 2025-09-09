import sys
import os
import argparse
from Bio import PDB 
import numpy as np
import string
import time
#from remove_incomplete_lipids import remove_incomplete_lipids

def get_lipid_files(lipids):
    """ Get the lipid files from the current directory

    Args:
        lipids (_list_): Iterates through the lipid types specified in the input and tries to find the corresponding pdb file in the current directory

    Returns:
        _list_: List of lipid files found in the current directory
    """
    lipid_files = []
    for lipid in lipids:
        lipid_file = str(lipid) + ".pdb"
        if lipid_file not in os.listdir():
            print("Lipid file not found: " + lipid_file)
            sys.exit(1)

        os.system("gmx editconf -f " + lipid_file + " -o " + lipid_file + " -center 0 0 0" + " >> gmx_output.log 2>&1")
        lipid_files.append(lipid_file)

    return lipid_files


 
def constrict_lipids(lipid_files, xy_scale=0.6, z_scale=1.2):
    """ Constrict the lipids by scaling the x and y coordinates by xy_scale and the z coordinates by z_scale

    Args:
        lipid_files (_list_): List of lipid files to be constricted
        xy_scale (float, optional): Value to scale x and y coordinates by. Defaults to 0.8.
        z_scale (float, optional): Value to scale z coordinates by. Defaults to 1.5.

    Returns:
        _type_: List of new lipid files that have been constricted
    """
    new_lipid_files = []
    for lipid in lipid_files:
        constricted_file = "constricted_" + lipid
        new_lipid_files.append(constricted_file) 
        with open(lipid, 'r') as f, open(constricted_file, 'w') as g:
            lines = f.readlines()
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    x = x * xy_scale
                    y = y * xy_scale
                    z = z * z_scale
                    line = line[:30] + "{:8.3f}".format(x) + "{:8.3f}".format(y) + "{:8.3f}".format(z) + line[54:]
                    g.write(line)
                else: 
                    g.write(line)

    return new_lipid_files

def get_dims_at_z(pdb, z, tolerance=10):
    """ Get the dimensions of the system at a given z coordinate

    Args:
        pdb (_str_): System to get dimensions of
        z (_int_): Z coordinate to get dimensions at
        tolerance (int, optional): Distance on either size of z to include atoms for when finding the max and min coordinates. Defaults to 10.

    Returns:
        _float, float_: Tuple of the x and y dimensions of the system at the given z coordinate
    """
    parser = PDB.PDBParser(QUIET=True) 
    structure = parser.get_structure("protein", pdb)
    atom_coords = np.array([atom.coord for atom in structure.get_atoms()])

    atoms_near_z = np.array([atom for atom in atom_coords if abs(atom[2]-z) < tolerance]).reshape(-1, 3)
    if len(atoms_near_z) == 0:
        return [0, 0, 0, 0]

    min_x = np.min(atoms_near_z[:,0])
    max_x = np.max(atoms_near_z[:,0])
    max_y = np.max(atoms_near_z[:,1])
    min_y = np.min(atoms_near_z[:,1])

    return [min_x, max_x, min_y, max_y]

# Return size of system in pdb file
def get_dimensions(pdb):
    """ Get the dimensions of the system in the pdb file, measured by the max and min coordinates in the x, y, and z directions

    Args:
        pdb (_str_): PDB file to get dimensions of

    Returns:
        _float, float, float_: Tuple of the x, y, and z dimensions of the system in the pdb file
    """
    parser = PDB.PDBParser(QUIET=True) 
    structure = parser.get_structure("protein", pdb)
    atom_coords = np.array([atom.coord for atom in structure.get_atoms()])

    dim_x = np.max(atom_coords[:,0]) - np.min(atom_coords[:,0])
    dim_y = np.max(atom_coords[:,1]) - np.min(atom_coords[:,1])
    dim_z = np.max(atom_coords[:,2]) - np.min(atom_coords[:,2])

    return dim_x, dim_y, dim_z

# Turn into to string if int is more than 5 digits 
def numberToBase(n, b):
    """ Convert a number to a string in a given base

    Args:
        n (_int_): Number to convert
        b (_int_): Base to convert the number to

    Returns:
        _str_: String representation of the number in the given base
    """
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    s = "0123456789" + string.ascii_uppercase[:b-10] 
    return "".join([s[x] for x in digits[::-1]])

def check_clash(x, y, protein_min_x, protein_max_x, protein_min_y, protein_max_y):
    """ Check if a lipid is in the protein

    Args:
        x (_float_): X coordinate of the lipid
        y (_float_): Y coordinate of the lipid

    Returns:
        _bool_: True if the lipid is in any protein, False otherwise
    """
    
    if protein_min_x < x < protein_max_x and protein_min_y < y < protein_max_y:
        return True
    else:
        return False

def insert_lipids(lipid_files, box_size, lipid_ratios, protein_dims, outfile, z=0, buffer=0, atom_number=1, res_number=1, z_buffer=0.1):
    """ Insert lipids into the system

    Args:
        lipid_files (_list_): List of all lipid files to insert
        box_size (_int_): Size of the box to insert the lipids into
        lipid_ratios (_list_): List of ratios of each lipid to insert
        protein_dims (_dict_): Dictionary containing the dimensions of the protein at each z coordinate
        outfile (_str_): Output file to write the new system to
        z (int, optional): Z coordinate to insert lower leaflet at. Defaults to 0.
        buffer (int, optional): Buffer to put between lipids. Defaults to 0.
        atom_number (int, optional): Number to start the atom id counting at. Defaults to 1.
        res_number (int, optional): Number to start the residue id counting at. Defaults to 1.
        z_buffer (float, optional): Buffer to put between the two leaflets. Defaults to 0.5.

    Returns:
        _int, int_: Tuple of the atom number and residue number of the last atom inserted
    """

    total_inserted = 0

    # Create a list of lipids based on the ratios
    lipids = [lipid_files[i] for i in range(len(lipid_files)) for _ in range(int(lipid_ratios[lipid_files[i]]))]
    #print(lipids)
    # Start at the bottom left corner of the box
    x, y = -box_size//2, -box_size//2
    first_leaflet = True
    with open(outfile, 'a') as out: 
        # Iterate until the box is full
        while(True):
            # Choose a random lipid to insert based on ratios
            curr_lipid = np.random.choice(lipids)
            lip_x, lip_y, lip_z = get_dimensions(curr_lipid)

            with open(curr_lipid, 'r') as f: 
                lines = []
                bad_pos = False
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        new_x = float(line[30:38]) + x
                        new_y = float(line[38:46]) + y
                        new_z = float(line[46:54]) + z
                        # If this is the first leaflet, flip the lipid upside down
                        if first_leaflet:
                            new_z = -new_z + 2*z

                        # Check if the lipid is in the protein at the given z coordinate
                        for protein_dim in protein_dims:
                            lipid_height = lip_z 
                            dim = protein_dim[int(2 * round(z/2))]
                            dim_top = protein_dim[int(2 * round((z + lipid_height/2)/2))]
                            dim_bottom = protein_dim[int(2 * round((z - lipid_height/2)/2))]
                            if check_clash(new_x, new_y, dim[0], dim[1], dim[2], dim[3]) or check_clash(new_x, new_y, dim_top[0], dim_top[1], dim_top[2], dim_top[3]) or check_clash(new_x, new_y, dim_bottom[0], dim_bottom[1], dim_bottom[2], dim_bottom[3]): 
                                bad_pos = True
                                break
                            # If lipid is in protein, skip this lipid
                        if bad_pos:
                                break
                        else: 
                            # If the atom number is more than 5 digits, convert to base 20
                            if len(str(atom_number)) > 5: 
                                a_num = numberToBase(atom_number, 20)
                                lines.append(line[:6] + a_num.rjust(5) + line[11:22] + "{:4d}".format(res_number) + line[26:30] + "{:8.3f}".format(new_x) + "{:8.3f}".format(new_y) + "{:8.3f}".format(new_z) + line[54:])
                            else:
                                # Otherwise, just write the atom number
                                lines.append(line[:6] + "{:5d}".format(atom_number) + line[11:22] + "{:4d}".format(res_number) + line[26:30] + "{:8.3f}".format(new_x) + "{:8.3f}".format(new_y) + "{:8.3f}".format(new_z) + line[54:])

                        atom_number += 1

                # Add the lipids to the system that aren't in the protein
                if not bad_pos:
                    for line in lines:
                        out.write(line)
                    res_number += 1
                    total_inserted += 1

            # Move to the next position
            x += lip_x + buffer
            # If we reach the end of the box, move to the next row
            if x > box_size//2:
                x = -box_size//2
                y += lip_y + buffer

                if y > box_size//2:
                    # If we finish the first leaflet, move to the second leaflet
                    if first_leaflet:
                        x = -box_size//2
                        y = -box_size//2
                        z = z + np.mean([get_dimensions(lip)[2] for lip in lipid_files]) + z_buffer 
                        first_leaflet = False
                    else:
                        # Once we've finished both leaflets, we're done 

                        print(f"Inserted {total_inserted} lipids")
                        break

    return atom_number, res_number

def insert_protein(protein_file, outfile, atom_number, res_number):
    """ Insert the protein into the system

    Args:
        protein_file (_str_): PDB file containing the protein to insert
        outfile (_str_): Output file to write the new system to
        atom_number (_int_): Atom number to start counting at
        res_number (_int_): Residue number to start counting at

    Returns:
        _int, int_: Tuple of the atom number and residue number of the last atom inserted
    """
    with open(protein_file, 'r') as f, open(outfile, 'a') as out:
        for line in f:
            if line.startswith("ATOM"):
                new_atom = atom_number
                new_res = res_number
                out.write(line[:6] + "{:5d}".format(new_atom) + line[11:22] + "{:4d}".format(new_res) + line[26:])
                new_atom += 1
                
            if line.startswith("TER"):
                new_res += 1
                out.write(line[:6] + "{:5d}".format(new_atom) + line[11:22] + "{:4d}".format(new_res) + line[26:])
                new_atom += 1
    
    return new_atom, new_res


def run_builder(proteins, lipids, lipid_ratios, output, box_size, z=0, buffer=2, z_buffer=0.1, xy_constrict=0.7, z_constrict=1):
    """ Run the builder

    Args:
        proteins (_list_): List of protein files to insert
        lipids (_list_): List of lipid files to insert
        lipid_ratios (_list_): List of ratios of each lipid to insert
        output (_str_): Output file to write the new system to
        box_size (_float_): Size of the box to insert the lipids into
        z (_float_): Z coordinate to insert lower leaflet at.
        buffer (_float_): Buffer to put between lipids.
        z_buffer (_float_): Buffer to put between the two leaflets.
        xy_constrict (_float_): Value to scale x and y coordinates by.
        z_constrict (_float_): Value to scale z coordinates by.
    """

    start_time = time.time()

    if len(lipids) != len(lipid_ratios):
        print("Number of lipid ratios must match number of lipids")
        sys.exit(1)

    with open(output, 'w') as f:
        f.write("REMARK    THIS IS A SIMULATION BOX\n")

    os.system("rm -f gmx_output.log")
    lipid_files = get_lipid_files(lipids)
    constricted_lipid_files = constrict_lipids(lipid_files, xy_scale=xy_constrict, z_scale=z_constrict)
    lipid_ratios = {constricted_lipid_files[i]: lipid_ratios[i] for i in range(len(constricted_lipid_files))}
    print(f'Getting protein dimensions for {len(proteins)} proteins')
    protein_dims = [{i: get_dims_at_z(protein, i) for i in range(int(2*round(int(-box_size)/2)), int(2*round(int(box_size)+2/2)), 2)} for protein in proteins]

    print('Inserting proteins')
    atom_num, res_num = 1, 1
    for protein in proteins:
        atom_num, res_num = insert_protein(protein, output, atom_number=atom_num, res_number=res_num)

    print('Inserting lipids')
    atom_num, res_num = insert_lipids(constricted_lipid_files, box_size, lipid_ratios, protein_dims, output, z, buffer, atom_num, res_num, z_buffer)
    os.system(r'rm \#*')

    print("Prepping for chimera-x")
    os.system("gmx editconf -f " + output + " -o " + output + " >> gmx_output.log 2>&1")
    os.system(f"sed -i '' '/OXT/d' {output}")

    print(f"Finished in {time.time() - start_time:.2f} seconds")
    
