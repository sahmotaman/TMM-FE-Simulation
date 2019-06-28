#TMM-FE-Simulation
This repository is created for sharing Abaqus material subroutines and finite element (FE) models for thermo-micro-mechanical (TMM) FE simulation of metal forming processes.


For TMM-FE Simulation:

- Create a (working) directory for your simulation.

- Copy the suitable "*.cae" and/or "*.inp" file(s) from "Abaqus FE Models" folder to the working directory.

- Copy "material_data.inp" from "TMM Material Properties - 20MnCr5" folder to the working directory.

- Here, Fortran scripts including Abaqus user material subroutines have been written using free-format Fortran 90 syntax. Therefore, for error-free compilation, flag "free" needs to be added to the Fortran compiler options in the Abaqus default environment file. The other way is to copy suitable "abaqus_v6.env" file from the folder "Abaqus Environment Files" to the working directory.

- In case of using Linux operating system, the extension of corresponding Abaqus material subroutine file located in "Abaqus User Material Subroutines" folder from "*.for" needs to renamed to "*.f".

- In case of using Abaqus CAE interface for running the simulation, open the "*.cae" file in the working directory. Edit the existing job in the "Job Manager". In the tab "General", point to the location (path) of the corresponding material user subroutine file which for Abaqus Explicit and Abaqus Standard are named "umat_subroutine.*" and "vumat_subroutine.*", respectively. Now, the job is ready for run.

- In case of using Abaqus Command interface, first change the directory to the working directory (cd command), and then execute the following command for job submission:
abaqus interactive job=[name of *.inp file] user=[path of material subroutine file]


Further notes:

- The consistent units systems used in this work is the "ton.mm.s" system. For obtaining more information regarding the units of various quantities in different consistent systems of units, refer to "units_system.pdf" file in "Abaqus FE Models" folder.

- The solution-dependent variables (SDVs) are described in "Abaqus FE Models/SDVs.txt".

- The FE models and material subroutines are created for Abaqus 2017. Thus, "*.cae" files may only be opened by Abaqus CAE 2017 or later versions of Abaqus CAE. Moreover, in order to use Abaqus subroutines, suitable version of Intel Fortran Compiler must be installed and linked to Abaqus.


S. Amir H. Motaman
