# DFT-Turbomole

## Description

This WaNo allows the users to perform periodic DFT calculations using the [VASP](https://www.vasp.at/) code without requiring a deep understanding of VASP's functionalities and input files. The user can choose between single point calculations, structure optimisations, ab-initio MD and NEB calculations. Electronic DOS and band structure calculations are also possible, for which it is recommended to use this WaNo several times in a corresponding workflow. The mandatory input of this WaNo consists only of a structure file in POSCAR format (except for NEB calculations), and all other input files can be created automatically. 

## Server setup

The code is based on python and the necessary virtual environment on the server is provided by Simstack. In addition, a configuration file called ```vasp.config``` needs to be provided in ```$NANOMATCH/$NANOVER/configs``` on the server where the environment variable $NANOVER can be adjusted in the xml file. The config file is responsible for setting up the VASP environment on the cluster which can be done by loading a corresponding module or by defining all necessary environment variables. In addition some WaNo-specific variables must be set to define the necessary commands and paths (see ```vasp.config``` for an example).

## Required input

The only required input is an initial structure in POSCAR format.

## WaNo Settings

- **VASP version**  
The chosen version must be available on the compute cluster. Please adjust the xml and ```vasp.config``` files accordingly.

- **Structure file**
This is the only mandatory input (except for NEB calculations; see above), and the structure must be given in POSCAR format.

- **Plane waves file**  
The user can choose between uploading a suitable POTCAR file or to rely on of the default. If the chosen PAW type is not available for a given element, the WaNo automatically looks for an alternative.

- **K-sampling**  
Different options are possible, with default being to define the minimum distance between two k-points in reciprocal space.

- **Electronic structure settings**  
Here, the parameters of DFT calculation are set up. Please check the [VASP manual](https://www.vasp.at/wiki/index.php) for further information on the different options. If a vDW functional is chosen, the ```vdw_kernel.bindat``` file must additionally be provided.

- **NEB settings**  
In case of a NEB calculation the user must provide either the initial and final structure of the path, so that the WaNo can linearly interpolate in order to create the starting images. Alternatively, a tar file can be provided which must contain folders ```00``` to ```xx``` containg each a starting image file named ```POSCAR```. 

- **Output files*  
Check the corresponding boxes if the CHGCAR and WAVECAR files are to be reused for analysis or restarting a calculation, otherwise empty files will be written.

## Output

The output of this WaNo consists of the VASP ouput files CONTCAR, CHGCAR, WAVECAR as well as IBZKPT. In addition, an archive file with all relevant output files is written. The file ```output_dict.yml``` contains information about snapshots of an MD calculation and is empty otherwise.

