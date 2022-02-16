# Server installation instructions

for Windows 10 64bit OS
Modified guide from https://omicx.cc/posts/2021-04-20-install-pymol-windows/

1.	Install Anaconda or Miniconda (Anaconda comes with several packages preinstalled).
For Windows 10 64-bit version, miniconda3-py3.9.7 (here, Conda 4.11.0 Python 3.9.7) was installed.
https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

2.	Create PyMOL environment
Click the Start menu and run “Anaconda Prompt (miniconda3)” “as administrator”.
Create a virtual environment named pymol:
conda create -n pymol

Activate the new environment
conda activate pymol
	
3.	Install supporting packages
conda install -c conda-forge pip numpy pmw
conda install -c anaconda pywget
conda install -c speleo3 tmalign
conda install -c speleo3 pymol-psico

4.	Download PyMOL whl files
There are several pre-compiled Open-Source PyMOL distributions from "Christoph Gohlke of the Laboratory for Fluorescence Dynamics, University of California, Irvine":
https://www.lfd.uci.edu/~gohlke/pythonlibs/#pymol-open-source
Here, we download the following files (both “pymol.whl” and “pymol_launcher.whl”):
pymol-2.6.0a0-cp39-cp39-win_amd64.whl
pymol_launcher-2.1-cp39-cp39-win_amd64.whl

5.	Install whl files
•	In the conda pymol environment, switch to download directory (e.g., D:\Downloads)
* Switch to Drive D:
d:
* Enter the downloads directory
cd Downloads

•	Install pymol_launcher first
pip install --no-index --find-links="%CD%" pymol_launcher-2.1-cp39-cp39-win_amd64.whl

PyMOL is installed automatically.
To update PyMOL later, run
pip install --upgrade --no-deps pymol-2.6.0a0-cp39-cp39-win_amd64.whl

6.	Launch PyMOL
In an activated pymol environment, run
pymol

# Setup and usage instructions

The application consists of the main script
•	apoholoR_web.py
and a helper script that downloads a UniProt-PDB chain mapping file and compiles 2 helper files
•	download_n_modify_SIFTS.py
Before the main script is run, the helper script must be run so that the UniProt file is downloaded, and the helper files are compiled and saved to the designated folder of the server. After the first run, the helper script should be run once a week to update the existing UniProt files. The script replaces the old files automatically.

The helper script needs no arguments to run, other than setting the root folder for the server, according to which it will create a subdirectory with the helper files (e.g. root/SIFTS).
As soon as it finishes running, the main script can be run.
The main script has a number of user parameters that have preset defaults within the script, and currently needs a single line input from the user to run. This single line input can carry anywhere from 1 to 3 arguments (right now), which are separated by whitespace characters. The maximum arguments within the single line input are of this form:
<pdb_id>  <chains>  <ligands>
  
pdb_id:
this is a 4-character code of a PDB protein structure. This argument is obligatory and only 1 can be input per line. (i.e. “1a73” or “3fav” or “3FAV”)
chains:
A single chain or multiple chains separated by commas or “ALL” in the case of all chains (i.e. “A” or “A,C,D” or “ALL”). This argument is non-obligatory, if omitted, all chains will be considered.
ligands:
A single ligand, multiple ligands separated by commas but without whitespace, or no ligands can be input per line (i.e. “HEM” or “hem” or “ATP” or “ZN” or “HEM,ATP,ZN”). This argument is non-obligatory, if omitted, the user should activate the automatic detection of the ligands in the structure from the available option, unless the user is starting with an apo structure, in which case they will need to activate the reverse mode (search for holo from apo).

Example of an input argument (query):
“1a73 A ZN”
The application will fetch the structure 1a73, get chain A, and look for zinc+2 (ZN) ligands to verify the input argument. If ZN is found in chain A of 1a73, it will retrieve all other known chains that are the same with chain A of 1a73, it will align them with 1a73_A and look for ZN ligands at the superimposed binding site of ZN on 1a73_A. If it finds protein chains with ZN, it will list them as HOLO, if the superimposed site is empty of ligands, the chain will be listed as APO. If another ligand is detected on that site instead of ZN, the chain will be listed as APO or HOLO, depending on the user’s preferences (if the user wants APO with no other ligands there, it will be listed as HOLO, and if the user does not mind other ligands in this binding site, it will be listed as APO).
