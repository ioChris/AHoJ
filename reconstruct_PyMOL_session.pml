# -*- coding: utf-8 -*-
"""
@author: Christos
"""

''' This script builds the PyMOL session by
    i) loading the structures in the results folder (*.cif.gz)
    ii) annotating the ligand selections (by parsing the ligands.csv file in the results)

    Note: the script should be loaded from the root directory of the files through PyMOL
'''

import glob
import sys

# Start python session
python

# Specify ligands file
ligands_file = 'ligands.csv'

# Find query file & struct
query_file = None
try:
    for filename in glob.glob('query_*'):
        query_file = filename
        query_object = query_file[:10]
        query_struct = query_file[6:10]
except Exception():
    sys.exit(1)

# Put ligands into list (bulk)
ligands_list = list()
if os.path.exists(ligands_file):
    with open (ligands_file, 'r') as ligs_in:
        ligs_header = ligs_in.readline()
        for line in ligs_in:
            ligands_list.append(line[:-1])
else:
    print('Ligands file [ligands.csv] not found')

# Get all structure filenames in the results folder
all_files = glob.glob('*.cif.gz')

# Put ligands into dict (struct:positions)
if ligands_list:
    ligands_dict = dict()
    ligands_dict_bulk = dict()

    for lig_chain in ligands_list:

        structchain = lig_chain.split(',')[0]
        positions = lig_chain.split(',')[1]
        positions = positions.split('-')

        ligands_dict_bulk.setdefault(structchain[:4], []).append(positions)  # This groups all ligands under same PDB code, not per chain
        ligands_dict.setdefault(structchain, []).append(positions)


    # Remove duplicate values from ligands_dict_bulk (preserve the order of values)
    for struct, bulk_positions in ligands_dict_bulk.items():
        unpack_positions = sorted(sum(bulk_positions, [])) # Unpack list of lists into a single list
        ligands_dict_bulk[struct] = list(ligands_dict_bulk.fromkeys(unpack_positions))

    # Remove duplicate values from ligands_dict (preserve the order of values)
    for struct, bulk_positions in ligands_dict.items():
        unpack_positions = sorted(sum(bulk_positions, [])) # Unpack list of lists into a single list
        ligands_dict[struct] = list(ligands_dict.fromkeys(unpack_positions))

elif os.path.exists(ligands_file) and not ligands_list:
	print('Ligands file [ligands.csv] does not contain any ligands')

# Load all structures into PyMOL session
for file in all_files:
    cmd.load(file)
    object_name = file[:-7]
    lig_positions = None

    if ligands_list: # If ligands file was found and was not empty

        if object_name.startswith('query'):  # Query
            struct = object_name[6:]

            # Get ligands for this object
            try:
                lig_positions = ligands_dict_bulk[struct]
            except Exception:
                continue

            if lig_positions is not None:
                for lig_position in lig_positions:
                    index = lig_position.split('_')[0]
                    chain = lig_position.split('_')[1]
                    resname = lig_position.split('_')[2]

                    # Restructure ligand selection name
                    lig_sel = resname + '_' + chain + index

                    # Create selection
                    cmd.select(lig_sel + '_query', object_name + ' and chain ' + chain + ' and resi ' + index + ' and resn ' + resname)

        else:  # Result

            structchain_result = object_name[-22:-17]
            structchain_query = object_name[-5:]

            # Get ligands for this object
            try:
                lig_positions = ligands_dict[structchain_result]
            except Exception:
                continue

            if lig_positions is not None:
                for lig_position in lig_positions:
                    index = lig_position.split('_')[0]
                    chain = lig_position.split('_')[1]
                    resname = lig_position.split('_')[2]

                    # Restructure ligand selection name
                    lig_sel = resname + '_' + chain + index

                    # Create selection
                    cmd.select(lig_sel + '_' + structchain_result + '-' + structchain_query, object_name + ' and chain ' + chain + ' and resi ' + index + ' and resn ' + resname)

cmd.deselect()
cmd.orient()

# End python session
python end
