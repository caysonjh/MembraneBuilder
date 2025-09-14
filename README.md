# MembraneBuilder
Library for embedding proteins in a membrane, especially for large systems.

## Installation 


```bash 
pip install git+https://github.com/caysonjh/MembraneBuilder.git#egg=MembraneBuilder
```

## Dependencies

GROMACS must be installed for the code to work correctly. Installation instructions can be found [here](https://manual.gromacs.org/current/install-guide/index.html)

## Usage 

The main function is `run_builder()` with the following arguments: 

 `proteins`: The name of the protein pdb files -- **Required unless generating a membrane only system**. 
  - These files should be in the same directory as the script. 
  - This option is configured to accept multiple proteins, and can thus embed multiple proteins into a membrane together. 
  - The protein pdb files *must* be in your desired position (i.e. centered vs located at other coordinates) prior to running the script. If you are unsure of the configuration, run `gmx editconf -f your_protein.pdb -center 0`. The output will include a line that indicates the center of your protein. It will automatically center the protein and output a file called `out.pdb`, which you can delete if centering the protein is not needed for what you are doing. 
  - If you are including multiple proteins, or a single protein with a hollow space meant to be filled with membrane, you will need to split the proteins into multiple files and provide them separately to the `--protein` option. The script is not sophisticated enough to recognize empty space within a protein's area and will not place lipids there. 

`lipids`: The names of the lipids the membrane will be composed of -- **Required**.
  - The lipid files corresponding to the names of provided here should be in the same directory as the script (i.e. if `--lipid POPC` is specified then `POPC.pdb` should be in the directory)
  - Any number of lipids can be specified, and can be inserted in any ratio, specified by `--lipid_ratios`.
  - The script will automatically center and prepare the lipid file, so no additional preparation is needed here.

`lipid_ratios`: The ratios of the lipids specified in `--lipids` -- **Required**.
  - The ratios should be specified in the same order as the lipids in `--lipids`.
  - The ratios should be integers, and there should be the same number specified as there are lipids in `--lipids`.
  - This will not affect the total number of lipids inserted into the membrane, but rather allows for membranes composed of multiple lipids.

`output`: The name of the output file -- **Required**.

`box_size`: The size of the simulation box to be created (in &Aring;ngstr&ouml;ms) -- **Required**.
  - This primarily determines the size of the membrane, and should be large enough to encompass the entirety of the protein. It will also directly determine how many lipids are inserted.

`z`: The z-coordinate where to begin inserting the membrane -- **Default: 0**. 
  - This is the z-coordinate where lipid insertion will begin, specifically the lower leaflet.
  - This option can be adjusted to change where in relation to your protein the membrane is placed. The default will insert the membrane at the center, but often membranes should be inserted at the polar region of the protein. 
  - This is often best figured out by trial and error, as the script will not automatically determine the best z-coordinate for the membrane.

`buffer`: The amount of buffer space to leave between each lipid during insertion (in &Aring;ngstr&ouml;ms) -- **Default: 2**. 
  - This will directly affect the density of the membrane. 
  - Adding a smaller buffer will increase membrane density, but will also increase the risk of atoms overlapping and causing infinite forces and other issues during simulation. 

`z-buffer`: The amount of buffer space to leave between the upper and lower leaflet in the membrane (in &Aring;ngstr&ouml;ms) -- **Default: 0.1**. 
  - This will directly impact how close the two membrane leaflets are to each other.
  - The script will automatically separate the two leaflets by the average height of each of the lipids included, but this is often not enough to separate the two leaflets. 
  - Because we want there to be as little gap between the leaflets as possible, this is often best figured out by trial and error. 

`xy_constrict`: The factor by which to constrict the bonds between atoms in the xy-plane -- **Default: 0.7**. 
  - This is used to reduce the size taken up by each of the lipids in the xy-plane, and thus increase the number of lipids that can be included without introducing atom overlap and infinite forces. 
  - This is often best determined by trial and error, until the desired lipid density is achieved.

`z_constrict`: The factor by which to constrict the bonds between atoms in the z-direction -- **Default: 1 (no constriction)**. 
  - This is used to reduce the size taken up by each of the lipids in the z-direction, and thus increase the number of lipids that can be included without introducing atom overlap and infinite forces.
  - This is often best determined by trial and error, until the desired lipid density is achieved.
