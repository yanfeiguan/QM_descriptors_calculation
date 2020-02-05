# QM_descriptors_calculation
This repository contains a QM descriptor generator which automatically calculate atomic/bond 
descriptors, e.g. partial charges, fukui indices, and bond index.

## Requirements
While it is possible to run all of the code on a normal desktop, HPC makes the 
calculation significantly faster. To run the code, you will need:
1. python 3.7
2. rdkit 2019
3. pandas
4. GFN2-xtb by Grimme https://xtb-docs.readthedocs.io/en/latest/contents.html
5. NBO for population analysis https://nbo6.chem.wisc.edu/, for NBO configurations see:
https://nbo6.chem.wisc.edu/INSTALL.gaussian
6. G16 for SCF calculation.

The first three can be installed through the following conda environment. 

## Installation
1. install Miniconda from https://conda.io/miniconda.html
2. create conda environment for QM_descriptors by:
    ```bash
   conda env create -f environment.yml
    ```
   ```bash
   source activate QM_descriptors
   ```

## Usage
To calculate the QM descriptors, run: ```python main.py --smiles <input.csv>```

For a list of optional flag, run: ```python main.py -h```

Two paths need to be specified in the main.py:

    XTB_PATH = '$GFN_XTB_PATH'
    G16_PATH = '$G16_PATH'
    
### Input
The code takes a .csv of smiles as input, e.g.:

    id,smiles
    0,CHEMBL231079,C
    1,CHEMBL1521196,C1=CC=CC=C1
    
The .csv must have these two columns with the same head. The id 
can be either str or int.

### Output
The code will generate three folders, as specified by the 
arguments '--MMFF_conf_folder', '--xtb_folder', '--DFT_folder', which hold results for 
MMFF conformer searching, semi-empirical optimization, and DFT electronic structure calculations, 
respectively.

The QM descriptors from DFT calculations will be parsed automatically and saved as a dataframe in the 
.pickle file specified by '--output' argument.

### Use on the HPC
See submit.sh for an example of submitting scripts. The following parameters in the 
submit.sh need to be replaced:

    $CONDA_PYTHON_PATH: python path in your conda environment, e.g.: /home/yanfeig/miniconda3/envs/QM_descriptors/bin
    $NBOPATH: path to the NBO bin folder, e.g.: /home/yanfeig:/home/yanfeig/nbo6/bin
    $CONDA_PACKAGE_PATH: package path in your conda environment, e.g.: /home/yanfeig/miniconda3/envs/QM_descriptors/lib/python3.7/site-packages
    $G16ROOT: path to the G16 root
    $SCRATCH_FOLDER_G16: Scratch folder holding G16 calculation 
    
## Known limitations
Allowed atom types currently are limited to the MMFF94 force field, which is used to search conformers for the 
input molecules. Allowed atoms and ions are: C, H, N, O, F, Si, P, S, Cl, Br, I,
 Fe+2, Fe+3, F-, Cl-, Br-, Li+, Na+, K+, Zn+2, Ca+2, Cu+1, Cu+2, and Mg+2

## Contributors
Yanfei Guan, Duminda Ranasinghe, Oscar Wu