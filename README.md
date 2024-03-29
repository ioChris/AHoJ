# Apo-Holo Juxtaposition

![example workflow](https://github.com/ioChris/AHoJ/actions/workflows/ci-conda.yml/badge.svg)


AHoJ is a structural bioinformatics tool that allows automated search and alignment of APO and HOLO structures for a given HOLO structure (and vice versa) in the PDB.
It can be used with single or multiple searches for creating customized Apo-Holo datasets.

In AHoJ the user can define the ligand of interest in a structure and find apo or holo structures according to this/these ligands.
The ligand can be defined either by selecting a binding residue or the ligand itself. AHoJ will mark the binding site of the ligand(s), find other structures that belong to the same protein and feature the same binding site(s), superimose them, and search for ligands. It will then decide whether each structure is apo or holo.

### Getting started
The main mode of search in AHoJ is "focused search". It is based on a user-specified point of interest (ligand, water molecule, modified residue) and it looks for apo and holo variations (structures) for this chunk of protein, based on the presence or absence of the specified point of interest at the particular binding site. Specifying the point of interest is therefore central to the workflow. Additional parameters can be controlled in by the user.

#### Specifying the ligand/point of interest
This is done by providing a series of arguments in the form of a query: the i)structure, ii)chain(s), iii)name of ligand or residue name (using the PDB 3-character code) and iv)position of ligand or residue in the sequence (PDB residue index is used). Specifying all four arguments is optional but recommended as it avoids ambiguity in the case that more than one ligand molecules of the same type exist within the same chain. When specifying a binding residue, it is obligatory to specify the position.
However, the minimum number of arguments is one - the structure (as PDB code). In such case, AHoJ will automatically detect the chains, the ligands and their positions. This mode offers flexibility but it is not recommended as any ligand will be detected in the query structure (or chain -if specified), which might be undesirable, and also the eventual apo-holo classification of the candidates will not depend on a single ligand, but possibly multiple, which can complicate the interpretation of results. AHoJ can work with ligands that are designated as heteroatoms in the PDB, which means small and medium-sized ligands (including post-translationally modified residues and water molecules) but not protein subunits.

Therefore, the main search mode starts with a holo structure, where the user knows the ligand that will be used as a starting point. This user-specified ligand will define the binding site as well as the search and annotation of the results. Any other (non-specified) ligands in the query structure will be ignored. In the case that the user does not know the ligand or the binding residues, AHoJ will automatically detect available ligands in the query structure.
If the query structure however does not bind any ligands (apo), AHoJ will still look for structures that belong to the same protein with the query, but it will not focus on a particular binding site. Instead it will report ligands that bind anywhere in the candidate structures.


##  Objective
**Given the holo form of a protein** (the complex of a protein structure **and its user-specified ligand(s)**):
1. **find** all unbound forms of that protein (**apo-proteins**) that lack the specified ligand(s) in the superimposed binding site(s) of the specified ligand(s).
2. **find** all other bound forms of that protein (**holo-proteins**) that include the specified ligand(s) or other ligands in the superimposed binding site(s) of the specified ligand(s).
3. **align both apo- and holo-proteins** to the query holo protein

**Given the apo form of a protein** (a protein structure **without any ligands**): 
1. **find** all other forms of that protein that bind at least one ligand (**holo-proteins**).
2. **find** all other unbound forms of that protein (**apo-proteins**) that completely lack ligands.
3. **align both apo- and holo-proteins** to the query apo protein

Note: ligands are confined to chemical components that are not part of the protein according to the PDB 
(however modified amino acids (e.g. phosphorylated), and water molecules can also be specified and considered as ligands). 
Ligands are named and defined according to their 1-3 character code from PDB (i.e. ATP, HEM, ZN).

##  Methodology
The application is conducting a search within the experimentally determined protein structures in the Protein Data Bank. 
It retrieves identical proteins with known structures, and then looks whether the specified ligand(s) are present or absent in the particular binding site(s). 

Depending on the user’s input, the application can look for the specified ligand(s) in a single specified chain, or the whole protein (all chains).
It can find apo-proteins that i) simply lack the specified ligand in a given binding site (but may bind a different ligand at the same site for example), or ii) lack any ligand in any part of the structure and thus constitute universal apo structures.

##  Requirements
The application was built and initially ran in a Windows 10 computer. The python packages were all installed through the Anaconda package manager.
The python version and packages used to run the application appear below:

* Python  3.9
* PyMOL		 2.4.1
* pywget		0.31
* wget		  3.2
* tmalign 20170708 
* pymol-psico 3.4.19
* lxml     4.8.0


Note: Newer versions of these packages should be functional as long as they are inter-compatible, including the open-source version of PyMOL.

Aside of the package dependencies, the application is using two precompiled files that are based on the SIFTS residue-level mapping between the UniProt sequence and the PDB structure.

* `uniprot_segments_observed_REVERSE_SPnum.bin`
* `uniprot_segments_observed_dict.bin`
                                 
These files need to be generated by `prepare.py` before running the main script `apoholo.py`.


## Installation instructions

### Windows

for Windows 10 64bit OS

##### 1.  Install Anaconda or Miniconda (Anaconda comes with several packages preinstalled).
For Windows 10 64-bit version, miniconda3-py3.9.7 (here, Conda 4.11.0 Python 3.9.7) was installed.
https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe


##### 2.  Create a new environment
Click the Start menu and run “Anaconda Prompt (miniconda3)” “as administrator”.
~~~sh
conda create -n ahoj      # Create a virtual environment named ahoj
conda activate ahoj       # Activate the new environment
~~~

##### 3.  Install supporting packages
~~~sh
conda install -c conda-forge pip numpy pmw
conda install -c anaconda pywget
conda install -c speleo3 tmalign
conda install -c speleo3 pymol-psico
conda install -c conda-forge rich pytest-sugar
conda install -c conda-forge lxml 
~~~     

##### 4.  Install Open-source PyMOL
There are several pre-compiled Open-source PyMOL distributions from "Christoph Gohlke of the Laboratory for Fluorescence Dynamics, University of California, Irvine":
https://www.lfd.uci.edu/~gohlke/pythonlibs/#pymol-open-source
Here, we use the following files:
~~~
pymol‑2.5.0‑cp39‑cp39‑win_amd64.whl
pymol_launcher-2.1-cp39-cp39-win_amd64.whl
~~~

In the conda `ahoj` environment, navigate to the download directory (e.g., `C:\Downloads`)
* Enter the downloads directory (`cd Downloads`)
* Install pymol_launcher (PyMOL will be installed automatically)
~~~sh
pip install --no-index --find-links="%CD%" pymol_launcher-2.1-cp39-cp39-win_amd64.whl
~~~

To update PyMOL later, run
~~~sh
pip install --upgrade --no-deps pymol-2.5.0-cp39-cp39-win_amd64.whl
~~~

### Windows (alternative)

Clean installation on Windows using `conda` (Miniconda3) package manager. 
       
1. Install conda manually with installer https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
2. Execute the following commands in cmd or PowerShell:
~~~sh
# clone repo
git clone git@github.com:ioChris/AHoJ.git && cd AHoJ

# init conda env
conda env create -f conda-env.yml 
conda activate ahoj

# prepare
python prepare.py  

# test
pytest -v 
 
# use
python apoholo.py --query "1a73 A,B ZN"
~~~

Tested on Windows 10.

### Linux
   
Clean installation on Linux using `conda` (Miniconda3) package manager. 

~~~sh
# install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh
bash Miniconda3-py39_4.11.0-Linux-x86_64.sh -b && bash

# clone repo
git clone git@github.com:ioChris/AHoJ.git && cd AHoJ

# init conda env
# to remove env: conda deactivate && conda env remove -n ahoj
conda env create -f conda-env.yml 
conda activate ahoj

# prepare
python prepare.py   

# test
pytest -v 

# use
python apoholo.py --query "1a73 A,B ZN"
~~~

Tested on Ubuntu 21.10 and Debian 9.        
      

## Setup and usage instructions

The application consists of the main script `apoholo.py`
and a helper script that downloads a UniProt-PDB chain mapping file and compiles 2 helper files: `prepare.py`.

Before the main script is run, the helper script must be run so that the UniProt file is downloaded, 
and the helper files are compiled and saved to the designated folder of the server. 
After the first run, the helper script should be run once a week to update the existing UniProt files. 
The script replaces the old files automatically.

The helper script needs no arguments to run, other than setting the root global working folder, 
according to which it will create a subdirectory with the helper files (e.g. root/SIFTS).

The main script has a number of user parameters that have preset defaults within the script, 
and currently needs a single line input (= query) from the user to run. 

             
### Query format

The maximum arguments within the single line input are of this form:
~~~
<pdb_code> <chains> <ligand_name> <ligand_position>
~~~  

* `pdb_code`: This is the 4-character code of a PDB protein structure (case-insensitive). This argument is obligatory and only 1 PDB code can be input per line. (i.e. “1a73” or “3fav” or “3FAV”). If it is the only argument (e.g. because the user does not know the ligand that binds to the structure), it will trigger automatic detection of ligands in the structure.
* `chains`: A single chain or multiple chains separated by commas (without whitespace), or “!” in the case of ligand-binding-only chains, or “\*” in the case of all chains (i.e. “A” or “A,C,D” or “!” or “\*”). This argument is case-sensitive and it is obligatory if the user intends to provide any argument after that (i.e. ligands or position). Note that "!" is preffered to "*", when the user does not know the chain of the ligand, as it allows the application to detect and process only the ligand-binding chains of the query structure.
* `ligand_name`: This argument is case-insensitive. A single ligand, multiple ligands separated by commas (without whitespace), or no ligands can be input per line (i.e. “HEM” or “hem” or “ATP” or “ZN” or “HEM,ATP,ZN”) or “\*” for the automatic detection of all ligands in the specified chain(s). Besides specifying the ligand directly by its name (and optionally, its position), the user can also specify a residue that binds the ligand by its position (i.e. “HIS 260”) and AHoJ will detect the ligand (as long as it is within 4.5 Angstroms of the residue). This approach however can lead to the selection of more than one ligands, if they are within this radius from the specified residue. This argument is non-obligatory, if omitted or specified as “\*”, AHoJ will automatically detect the ligands in the structure. If there are no ligands in the query structure, it will be characterised as apo and the search for candidates will continue in a non-binding-site-specific manner. A water molecule can also be specified as a ligand (i.e. “HOH”) but in such cases, its position must be specified as well. Note: when specifying the `position` argument, the user can only specify one ligand per query.
* `ligand_position`: This argument is an integer (i.e. “260” or “1”). It refers to the PDB index of the previously specified ligand, binding residue or water molecule. When this argument is specified, only one ligand or residue can be specified in the previous argument.
                       
## Usage examples

Example of a user query:
~~~sh
# consider ZN ligand in position 201 in chain A of 1a73 
python apoholo.py --query '1a73 A ZN 201'
~~~

The application will fetch the structure 1a73, get chain A, and look for zinc+2 (ZN) ligand in position 201 of the sequence to validate the input. 
If ZN is found in chain A and position 201 of 1a73 (1a73A), it will retrieve all other known chains that belong to the same protein with 1a73A, 
it will align them with 1a73A and look for ZN (and also other ligands) at the superimposed binding site of ZN in 1a73A. If it finds protein chains with ZN, 
it will list them as HOLO, if the superimposed site is empty of ligands, the chain will be listed as APO. 
If another ligand is detected on that site instead of ZN, 
the chain will be listed as APO or HOLO, depending on the value of `--lig_free_sites` parameter (if the user wants APO with no other ligands there, 
it will be listed as HOLO, and if the user does not mind other ligands in this binding site, it will be listed as APO).

Example of an alternative query that leads to the same result with the previous example:
~~~sh
# consider ligands near residue HIS134 in chain A of 1a73 (the detected ligand will be ZN 201 in chain A)
python apoholo.py --query '1a73 A HIS 134'
~~~

### More examples of user queries

~~~sh 
# consider ZN ligands in chains A and B of 1a73
python apoholo.py --query '1a73 A,B ZN'

# consider ZN ligands in all chains of 1a73
python apoholo.py --query '1a73 * ZN'
# or
python apoholo.py --query '1a73 ! ZN'   # "!" is preferred to "*" as it only processes ligand-binding query chain(s)

# find and consider all ligands in all chains of 1a73
python apoholo.py --query '1a73'

# find and consider all ligands in chain A of 1a73
python apoholo.py --query '1a73 A'

# consider ZN and MG ligands in chain A of 1a73
python apoholo.py --query '1a73 A ZN,MG'

# consider ZN ligands in all chains of 3fav 
python apoholo.py --query '3fav ZN'
~~~

##  Parameters

### Basic

**`--res_threshold`** : resolution threshold [default = `3.8`]

Floating point number that represents angstroms and is applied as a cutoff point when assessing candidate structures that are resolved by scattering methods (X-ray crystallography, electron microscopy, neutron diffraction). It applies at the highest resolution value, when this is available in the PDB structure file. It can take any value, suggested min/max = 1.5/8. Condition is <=

**`--include_nmr`** : include NMR structures [default = `1`]

0 or 1. When set to 1 (ON), NMR structures are considered as candidates. In the case of multiple states for a certain structure, the first one is considered.

**`--xray_only`** x-ray structures only [default = `0`]

0 or 1. When set to 1 (ON), only X-ray structures are considered. This overrides the NMR setting.

**`--lig_free_sites`** : ligand-free sites [default = `1`]

0 or 1. This pertains to the apo/holo classification of candidate structures. When set to 1 (ON), it does not tolerate any ligands (in addition to the user-specified one(s)) in the superimposed binding sites of the candidate apo structures. When set to 0 (OFF), it tolerates ligands other than the user-specified one(s) in the same superimposed binding site(s). If the user wants to find apo structures that don't bind any ligands in the superimposed binding site(s) of the query ligand(s), this value should be set to 1 (default).

[//]: # (**`--autodetect_lig`** : auto-detect ligands [default = `0`])

[//]: # ()
[//]: # (0 or 1. When set to 1 &#40;ON&#41;, the algorithm will auto-detect any possible ligands &#40;and their respective binding sites&#41; in the query protein. It will then conduct the apo-protein search according to these ligands and binding sites.)

[//]: # ()
[//]: # (**`--reverse_search`** [default = `0`])

[//]: # ()
[//]: # (0 or 1. When set to 1 &#40;ON&#41;, the algorithm will consider that the input is an apo structure without any ligands, essentially reversing the search. It will try to find structures &#40;chains&#41; of the same protein with and without ligands. It will list those as holo and apo respectively. This reverse search is broader in scope because there is no starting ligand &#40;and thus binding site&#41; as a reference point. Any structures that belong to the same protein and bind at least one ligand, will be characterised as holo.)

### Advanced

**`--bndgrsds_threshold`** : binding residues threshold [default = `1.0`, min/max = 1/100]

Floating point number that represents a percentage (%) and is applied as a minimum cut-off upon the percentage of the number of successfully mapped binding residues in the candidate chain out of the total number of binding residues in the query chain. The binding residues are mapped between query and candidate by converting PDB to UniProt numbering. "1%" translates to at least 1% percent of the query residues being present in the candidate structure, for the structure to be considered as apo or holo.

**`--save_apo`** : save aligned Apo chains [default = `1`]

0 or 1. When set to 1 (ON), saves the structure files of the aligned APO chains (mmCIF). Disabling this is only recommended in multiple queries if visualizations are not needed (reduces download size). This setting does not affect the search for apo or holo chains or the final result reports.

**`--save_holo`** : save aligned Holo chains [default = `1`]

0 or 1. When set to 1 (ON), saves the structure files of the aligned HOLO chains (mmCIF). Disabling this is only recommended in multiple queries if visualizations are not needed (reduces download size). This setting does not affect the search for apo or holo chains or the final result reports.

**`--overlap_threshold`** : sequence overlap threshold [default = `0`, min/max = 0/100]

Floating point number that represents a percentage (%) and is applied as a cutoff point when comparing the sequence overlap between the query and the candidate chain. It applies to the percentage of sequence overlap between query and candidate chains, and it is calculated from the query's perspective according to the UniProt residue numbering. If set to 100 (%), it means that the candidate chain has to completely cover the query chain. It can be longer than the query, but not shorter.

Note: "100" guarantees complete coverage, but it is the strictest setting. If the user wants a more lenient filtration, they can lower the value, or even set it to 0 and rely on the template-modeling score (TM-score) by using the default value (0.5) or setting their own TM-score cutoff with the "min_tmscore" parameter.

**`--lig_scan_radius`** : ligand scanning radius [default = `4.5`]

Floating point number that represents angstroms and is applied as a scanning radius when looking for ligands in the candidate structures. This scanning radius is applied on the positions of the atoms of the superimposed query ligands to the aligned candidate structure, to scan for ligands. The resulting scanning space is a "carved" surface that has the shape of the query ligand, extended outward by the given radius. If the candidate structure binds ligands outside of this superimposed area, they will be ignored, and the candidate will be characterised as an apo-protein.

**`--min_tmscore`** : minimum TM-score [default = `0.5`, min/max = 0/1]

Floating point number that is applied as a minimum accepted template-modeling score between the query and the candidate chain. Value 1 indicates a perfect match, values higher than 0.5 generally assume the same fold in SCOP/CATH.

**`--water_as_ligand`** : [default = `0`]

0 or 1. When set to 1 (ON), allows the detection of water molecules (i.e., 'HOH') as ligands in the superposition of the query ligand(s) in the candidate chains. If this setting is enabled and at least one water molecule is detected within the scanning radius, that would warrant a holo classification for the candidate chain. When a water molecule is defined in the user query, this setting is automatically enabled.

**`--nonstd_rsds_as_lig`** : non-standard residues as ligands [default = `0`]

0 or 1. When set to 1 (ON), allows the detection of non-standard -or modified- residues (e.g., 'TPO', 'SEP') as ligands in the superposition of the query ligand(s) in the candidate chains. If this setting is enabled and at least one modified residue is detected within the scanning radius, that would warrant a holo classification for the candidate chain. When a modified residue is defined in the user query, this setting is automatically enabled.

Note: The current list of non-standard residues includes the following residue names: 'SEP TPO PSU MSE MSO 1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG PYG PYL SEC PHA'.

**`--d_aa_as_lig`** : D-amino acids as ligands [default = `0`]

0 or 1. When set to 1 (ON\), allows the detection of D-residues (e.g., 'DAL', 'DSN'\) as ligands in the superposition of the query ligand(s\) in the candidate chains. If this setting is enabled and at least one D-residue is detected within the scanning radius, that would warrant a holo classification for the candidate chain. When a D-residue is defined in the user query, this setting is automatically enabled.

Note: The current list of D-residues includes the following residue names: 'DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA'.


##  Results

The results are visualised in the browser through Mol* and they can be downloaded as a zip file after a run has completed.

#### Files

In a successful run, AHoJ should generate the following files:

i) PDB structure files (cif.gz format) for the query structure (whole structure) and the successfully processed apo and holo candidate chains, aligned to the respective query chain(s).

Note: a given candidate chain could be a match for more than one query chains, and could thus appear more than once, in each case aligned to the respective query chain.

ii) 2 CSV files with the successfully processed candidate chains for apo and holo chains respectively [results_apo.csv, results_holo.csv]. These CSV files contain the following information for each found chain: **`query_chain, apo_chain, Resolution, R-free, %UniProt_overlap, Mapped_bndg_rsds, %Mapped_bndg_rsds, RMSD, TM_score, iTM_score, ligands`**

Note: UniProt overlap refers to the percentage of sequence overlap that each candidate chain has over the query chain (higher is better, max is 100%). The ligands listed in the files refer to the ligands that were detected in the superimposed positions of the specified query ligands, thus they might not include ligands that bind elsewhere in the candidate chains. If the CSV file for apo chains includes ligands (which seems contradicting), it indicates that the user set the parameter **`--lig_free_sites`** to 0 (OFF), and thus any other ligands besides the query ligand were allowed in the superimposed binding sites of candidate structures.

iii) 1 CSV file with the positions of the relevant ligands that were detected in both query and resulting candidate structures [ligands.csv]. This file is needed to load ligand selections when loading the results into the PyMOL with the included script.

iv) Console log file with information from the standard output [console.log]. This file can be used for reference and for better understanding the mechanism of action of AHoJ.

v) A PyMOL script file for loading the results into PyMOL [load_results_into_PyMOL.pml]. This is useful for viewing the results locally. The script has to be opened through PyMOL. The resulting session can then be saved by the user as a PyMOL session (.pse).

#### Visualization

The web application allows the visualization of the results in the browser with the molstar (Mol*) viewer. Web application: https://github.com/rdk/AHoJ-webapp

The results can also be downloaded and visualized locally by loading the PyMOL script that is included in the results folder through PyMOL [load_results_into_PyMOL.pml]. The script has to be loaded from within the results folder. After downloading and unpacking the results into a folder, start a new PyMOL session and open the .pml file through it. A PyMOL installation is needed for this to work (Incentive or Open-Source).

## Notes / Troubleshooting

### Case sensitivity in Windows

When running the application locally in Windows systems, it is possible to run into a situation where a structure has 2 chains of the same letter but in different case (i.e., "5cdwC" and "5cdwc"); this is typical in large EM assemblies that have numerous chains. In such situations, a saved aligned chain can be overwritten by the last alphanumerically-identical chain. This is because Windows by default treats file and directory names as case-insensitive. This can be changed for a particular directory in Windows 10 or later. To learn how change these settings, visit:
https://learn.microsoft.com/en-us/windows/wsl/case-sensitivity

Quick fix

To change case-sensitivity of a directory in the Windows file system, run PowerShell as Administrator and use the following commands:
~~~sh
# Enable case sensitivity in a Windows directory
fsutil.exe file setCaseSensitiveInfo <path> enable
# Disable case sensitivity in a Windows directory
fsutil.exe file setCaseSensitiveInfo <path> disable
~~~
