# apo_holo_webserver
 Webserver for the detection of apo structures from holo

##  Objective
Given the holo form of a protein (the complex of a protein structure and its specified ligand(s)), detect an unbound version of the protein (apo-protein) that is lacking the specified ligand(s).
Note: ligands are confined to chemical components that are not part of the protein according to the PDB*(phosphorylated amino acids can be considered and specified as ligands). They are named and defined according to their 3-character code from PDB (i.e. ATP, HEM, ZN).

##  Methodology
The application is conducting a search within the experimentally determined protein structures in the Protein Data Bank. It retrieves identical proteins with known structures, and then it looks whether the specified ligand(s) is present or absent. Depending on the user’s preferences, the application can look for the specified ligand(s) in a single specified chain, or the whole protein (all chains). It can find apo-proteins that i) simply lack the specified ligand in a given binding site (but may bind a different ligand at the same site), or ii) lack any known ligand in a given binding site and thus constitute universal apo- sites.

##  Server setup requirements
The application was built and initially ran in a Windows 10 OS. The python packages were all installed through the Anaconda package manager.
The python version and packages used to run the application appear below:

>Python  3.8.11

>PyMOL		 2.4.1

>pywget		0.31

>wget		  3.2

>tmalign 20170708

>pymol-psico 3.4.19


Note: Newer versions of these packages should be functional as long as they are intercompatible, including the open-source version of PyMOL.

Aside of the package dependencies, the application is using two precompiled files that are based on the SIFTS UniProt – PDB chain mapping.

>pdb_chain_uniprot_REVERSE_SPnum.txt

>pdb_chain_uniprot_dict.txt

