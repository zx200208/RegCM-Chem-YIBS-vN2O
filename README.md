# RegCM-Chem-YIBS-vN2O
This is a modified version of the coupled regional climate-chemical-ecological model RegCM-Chem-YIBs with the addition of nitrous oxide species.
The RegCM-Chem-YIBs-vN2O model enables coupling between climate, chemistry and ecology, which add the specie of nitrous oxide.

Following the steps below, the user can install the model and start using it.

code replace and add

Replace the source code under RegCM-Chem-YIBs/Main with the file source code mod_lm_interface.F90, mod_ncio.F90, mod_output.F90 in the code folder Replace the source code under RegCM-Chem-YIBs/Main/chemlib with the file source code in the code/chemlib folder Replace the source code under RegCM-Chem-YIBs/Main/mpplib Add the yibs directory (the directory contains all the YIBs program source code)in code folder to the RegCM-Chem-YIBs/Main directory install

./configure --enable-clm45 CC=icc FC=ifort --with-netcdf=/your netcdf directory ***CC and GCC depend on your compiler make make install Running

Create a run directory ../bin/terrainCLM45 namelist_example.in ../bin/sstCLM45 namelist_example.in ../bin/icbcCLM45 namelist_example.in ../bin/emcre_gridCLM45 namelist_example.in ../bin/interp_emissions namelist_example.in ../bin/chem_icbcCLM45 namelist_example.in ../bin/mksurfdataCLM45 namelist_example.in Add CO2 and N2O species to CHBC file ***CO2 species can be added to the CHB(file generated in step 7) using the cdo command (which you need to install first) Create a CO2 and N2O source input file ***You can use the cdo command (which you need to install first) in adding CO2 and N2O emission source files,which you can make with reference to the sample data Create the yibs input file ***You can refer to the example data yibs_driver.nc to create yibs input files using ncl ../bin/regcmMPICLM45 namelist_example.in
