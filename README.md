# Apo-Holo Juxtaposition

Python tool for search and alignment of apo structures for a given holo structure (and vice versa).

Tool for creating customized apo-holo datasets.


##  Objective
**Given the holo form of a protein** (the complex of a protein structure **and its user-specified ligand(s)**):
1. **find** all unbound forms of that protein (**apo-proteins**) that lack the specified ligand(s) in the superimposed binding site(s) of the specified ligand(s).
2. **find** all other bound forms of that protein (**holo-proteins**) that include the specified ligand(s) or other ligands in the superimposed binding site(s) of the specified ligand(s).
3. **align both apo- and holo-proteins** to the query holo protein

**Given the apo form of a protein** (a protein structure **without any ligands**): 
1. **find** all other forms of that protein that bind at least one ligand (**holo-proteins**).
2. **find** all other unbound forms of that protein (**apo-proteins**) that lack any ligands.
3. **align both apo- and holo-proteins** to the query apo protein

Note: ligands are confined to chemical components that are not part of the protein according to the PDB 
(phosphorylated amino acids can also be considered and specified as ligands). 
Ligands are named and defined according to their (up to) 3-character code from PDB (i.e. ATP, HEM, ZN).

##  Methodology
The application is conducting a search within the experimentally determined protein structures in the Protein Data Bank. 
It retrieves identical proteins with known structures, and then looks whether the specified ligand(s) are present or absent. 
Depending on the user’s preferences, the application can look for the specified ligand(s) in a single specified chain, or the whole protein (all chains). 
It can find apo-proteins that i) simply lack the specified ligand in a given binding site (but may bind a different ligand at the same site), or ii) 
lack any known ligand in a given binding site and thus constitute universal apo- sites.

##  Requirements
The application was built and initially ran in a Windows 10 computer. The python packages were all installed through the Anaconda package manager.
The python version and packages used to run the application appear below:

* Python  3.8.11
* PyMOL		 2.4.1
* pywget		0.31
* wget		  3.2
* tmalign 20170708 
* pymol-psico 3.4.19


Note: Newer versions of these packages should be functional as long as they are inter-compatible, including the open-source version of PyMOL.

Aside of the package dependencies, the application is using two precompiled files that are based on the SIFTS UniProt – PDB chain mapping.

* `pdb_chain_uniprot_REVERSE_SPnum.txt`
* `pdb_chain_uniprot_dict.txt`
                                 
These files need to be generated by `download_n_modify_SIFTS.py` before running the main script `apoholoR_web.py`.


## Installation instructions

### Windows

for Windows 10 64bit OS
Modified guide from https://omicx.cc/posts/2021-04-20-install-pymol-windows/

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
~~~     

##### 4.  Install Open-source PyMOL
There are several pre-compiled Open-source PyMOL distributions from "Christoph Gohlke of the Laboratory for Fluorescence Dynamics, University of California, Irvine":
https://www.lfd.uci.edu/~gohlke/pythonlibs/#pymol-open-source
Here, we download the following files:
~~~
pymol‑2.5.0‑cp39‑cp39‑win_amd64.whl
pymol_launcher-2.1-cp39-cp39-win_amd64.whl
~~~

In the conda ahoj environment, switch to download directory (e.g., D:\Downloads)
* Switch to Drive D:
* Enter the downloads directory (`cd Downloads`)
* Install pymol_launcher (PyMOL will be installed automatically)
~~~sh
pip install --no-index --find-links="%CD%" pymol_launcher-2.1-cp39-cp39-win_amd64.whl
~~~

To update PyMOL later, run
~~~sh
pip install --upgrade --no-deps pymol-2.6.0a0-cp39-cp39-win_amd64.whl
~~~


### Linux
   
#### Using conda

Clean installation on Linux using *conda* package manager. Tested on Ubuntu 21.10 and Debian 9.

~~~sh
# clone repo
git clone git@github.com:ioChris/AHoJ.git && cd AHoJ

# install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh
bash Miniconda3-py39_4.11.0-Linux-x86_64.sh && bash
conda update conda

# init conda env
conda env create -f conda-env.yml
conda activate ahoj

# prepare
python download_n_modify_SIFTS.py --work_dir '../ahoj_workdir'  
# test
python apoholo_J.py --work_dir '../ahoj_workdir' --single_line_input '1a73 A,B ZN'
~~~

#### Without conda

*TODO*

## Setup and usage instructions

The application consists of the main script
~~~
apoholo_J.py
~~~

and a helper script that downloads a UniProt-PDB chain mapping file and compiles 2 helper files
~~~
download_n_modify_SIFTS.py
~~~

Before the main script is run, the helper script must be run so that the UniProt file is downloaded, 
and the helper files are compiled and saved to the designated folder of the server. 
After the first run, the helper script should be run once a week to update the existing UniProt files. 
The script replaces the old files automatically.

The helper script needs no arguments to run, other than setting the root folder for the server, 
according to which it will create a subdirectory with the helper files (e.g. root/SIFTS).
As soon as it finishes running, the main script can be run.
The main script has a number of user parameters that have preset defaults within the script, 
and currently needs a single line input from the user to run. 
This single line input can carry anywhere from 1 to 3 arguments (right now), 
which are separated by whitespace characters. 

The maximum arguments within the single line input are of this form:
~~~
<pdb_id> <chains> <ligands>
~~~  

* `pdb_id`: This is the 4-character code of a PDB protein structure. This argument is obligatory and only 1 PDB ID can be input per line. (i.e. “1a73” or “3fav” or “3FAV”)
* `chains`: A single chain or multiple chains separated by commas (without whitespace), or “ALL” in the case of all chains (i.e. “A” or “A,C,D” or “ALL”). This argument is non-obligatory, if omitted, all chains will be considered.
* `ligands`: A single ligand, multiple ligands separated by commas (without whitespace), or no ligands can be input per line (i.e. “HEM” or “hem” or “ATP” or “ZN” or “HEM,ATP,ZN”). This argument is non-obligatory, if omitted, the user should activate the automatic detection of the ligands in the structure from the available option, unless the user is starting with an apo structure, in which case they will need to activate the reverse mode (search for holo from apo).

Example of an input argument (query):
~~~
python apoholo_J.py --single_line_input '1a73 A ZN'
~~~
-considers ZN ligands in chain A of 1a73


The application will fetch the structure 1a73, get chain A, and look for zinc+2 (ZN) ligands to verify the input argument. 
If ZN is found in chain A of 1a73 (1a73A), it will retrieve all other known chains that belong to the same protein with 1a73A, 
it will align them with 1a73A and look for ZN ligands at the superimposed binding site of ZN in 1a73A. If it finds protein chains with ZN, 
it will list them as HOLO, if the superimposed site is empty of ligands, the chain will be listed as APO. If another ligand is detected on that site instead of ZN, 
the chain will be listed as APO or HOLO, depending on the user’s preferences (if the user wants APO with no other ligands there, 
it will be listed as HOLO, and if the user does not mind other ligands in this binding site, it will be listed as APO).

More examples of user queries.
~~~
python apoholo_J.py --single_line_input '1a73 A,B ZN'
~~~
-considers ZN ligands in chains A and B of 1a73
~~~
python apoholo_J.py --single_line_input '1a73 ALL ZN'
~~~
or
~~~
python apoholo_J.py --single_line_input '1a73 ZN' (with ligand auto-detection OFF)
~~~
-considers ZN ligands in all chains of 1a73
~~~
python apoholo_J.py --single_line_input '1a73' (with ligand auto-detection ON or OFF)
~~~
-finds and considers all ligands in all chains of 1a73
~~~
python apoholo_J.py --single_line_input '1a73 A' (with ligand auto-detection ON)
~~~
-finds and considers all bound ligands in chain A of 1a73
~~~
python apoholo_J.py --single_line_input '1a73 A ZN,MG'
~~~
-considers ZN and MG ligands in chain A of 1a73
~~~
python apoholo_J.py --single_line_input '3fav ZN'
~~~
-considers ZN ligands in all chains of 3fav

##  Parameters

### Basic

**res_threshold** : resolution threshold [default = 3.5]

Floating point number that represents angstroms and is applied as a cutoff point when assessing the apo-holo candidate structures that are resolved by scattering methods (X-ray crystallography, electron microscopy, neutron diffraction). It applies at the highest resolution value, when this is available in the PDB structure file. It can take any value, suggested min/max = 1.5/8. Condition is <=

**NMR** [default = 1]

0 or 1. When set to 1 (ON), NMR structures are considered as candidates. In the case of multiple states for a certain structure, the first one is considered.

**x-ray only** [default = 0]

0 or 1. When set to 1 (ON), only X-ray structures are considered. It overrides the NMR setting.

**lig_free_sites** : ligand-free sites [default = 0]

0 or 1. When set to 1 (ON), it does not tolerate any ligands (in addition to the user-specified one(s)) in the superimposed binding sites of the candidate apo-proteins. When set to 0 (OFF), it tolerates ligands other than the user-specified one(s) in the same superimposed binding site(s).

**autodetect_lig** : auto-detect ligands [default = 0]

0 or 1. When set to 1 (ON), the algorithm will auto-detect any possible ligands (and their respective binding sites) in the query protein. It will then conduct the apo-protein search according to these ligands and binding sites.

**reverse_search** [default = 0]

0 or 1. When set to 1 (ON), the algorithm will consider that the input is an apo structure without any ligands, essentially reversing the search. It will try to find structures (chains) of the same protein with and without ligands. It will list those as holo and apo respectively. This reverse search is broader in scope because there is no starting ligand (and thus binding site) as a reference point. Any structures that belong to the same protein and bind at least one ligand, will be characterised as holo.

### Advanced

**save_separate** [default = 1]

0 or 1. When set to 1 (ON), the server will save all aligned chains that are in the opposite category from the starting query (apo/holo). In a regular search where the query is a holo-protein (searching for apo from holo), it will save any apo chains that it will find. In a reverse search, it would save all holo chains. If the user wishes to save both apo and holo chains, they can turn on the next parameter, "save_oppst".

**save_oppst** : save opposite [default = 1]

0 or 1. When set to 1 (ON), the server will not only find, but also save chains that are in the same category with the starting query (apo or holo). In a regular search where the query is a holo-protein (searching for apo from holo), it will also save any holo chains that it will find. In a reverse search (starting with an apo-protein and looking for holo), it will also save the apo chains that it will find. This setting is dependant on the previous parameter "save_separate" which has to be ON for this parameter to work. This setting does not affect the text output of the server which always includes both apo and holo chains.

**overlap_threshold** [default = 0, min/max = 0/100]

Floating point number that represents a percentage (%) and is applied as a cutoff point when comparing the sequence overlap between the query and the candidate chain. It applies to the percentage of sequence overlap between query and candidate chains, and it is calculated from the query's perspective according to the UniProt residue numbering. If set to 100 (%), it would mean that the candidate chain has to completely overlap with the query chain. It can be longer than the query, but not shorter.

Note: "100" guarantees a complete overlap, but it is the strictest setting. If the user wants a more lenient filtration, they can lower the value, or even set it to 0 and rely on the template-modeling score (TM-score) by using the default value (0.5) or setting their own TM-score cutoff with the "min_tmscore" parameter.

**lig_scan_radius** : ligand scanning radius [default = 5]

Floating point number that represents angstroms and is applied as a scanning radius for ligands, from the point of the superimposition of the query ligand(s) in the respective superimposed binding sites on the candidate chains. If the candidate has ligands bound outside of this sphere, they will be tolerated, and the candidate will be characterised as apo-protein.

**min_tmscore** : minimum TM-score [default = 0.5, min/max = 0/1]

Floating point number that is applied as a minimum cutoff point for template-modeling score between the query and the candidate chain. Value 1 indicates a perfect match, values higher than 0.5 generally assume the same fold in SCOP/CATH.  default = 0.5
