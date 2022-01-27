# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:24:57 2021

@author: ChrisX
"""
# Apo finder/ Apofind
'''Given a holo structure, with specified chain(s) and ligand(s), find apo structures
    Use PyMOL parser
    '''

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
import pymol
import pymol.cmd as cmd
pymol.finish_launching()
import psico.fitting

import ast
import gzip
import os
import wget
import time
import argparse
#from os import walk
import sys
#import csv


## User input
#multiline_input = '3fav all zn\n1a73 a zn,MG,HEM\n5ok3 all tpo'

#single_line_input = '3fav all zn'
#single_line_input = '1a73 a zn,MG,HEM'
#single_line_input = '5ok3 all tpo'

#'2ZB1 all gk4'
#'7l1f all F86'
single_line_input ='3CQV all hem,f86,mg,tpo,act' #hem
#single_line_input = '5gss all gsh' # slow

# Create the parser, add arguments
parser = argparse.ArgumentParser()
# User arguments
parser.add_argument('--res_threshold', type=float, default=3, help='resolution cut-off for apo chains (angstrom), condition is <=')
parser.add_argument('--NMR', type=int, default=1, help='0/1: discard/include NMR structures')
parser.add_argument('--lig_free_sites', type=int, default=1, help='0/1: resulting apo sites will be free of any other known ligands in addition to specified ligands')
parser.add_argument('--water_as_ligand', type=int, default=0, help='0/1: consider HOH atoms as ligands (can be used in combination with lig_free_sites)(strict)')
parser.add_argument('--save_session', type=int, default=1, help='0/1: save each result as a PyMOL ".pse" session (zipped, includes annotations -recommended)')
parser.add_argument('--multisave', type=int, default=0, help='0/1: save each result in a .pdb file (unzipped, no annotations -not recommended)')
# Internal/advanced arguments
parser.add_argument('--overlap_threshold', type=int, default=100, help='% of overlap between apo and holo chain (w UniProt numbering), condition is ">="')
parser.add_argument('--ligand_scan_radius', type=str, default='5', help='angstrom radius to look around holo ligand(s) superposition')
parser.add_argument('--apo_chain_limit', type=int, default=999, help='limit number of apo chains to consider when aligning (for fast test runs)')
parser.add_argument('--min_tmscore', type=float, default=0.5, help='minimum acceptable TM score for apo-holo alignments (condition is "<" than)')
parser.add_argument('--beyond_hetatm', type=int, default=0, help='0/1: when enabled, does not limit holo ligand detection to HETATM records for specified ligand/residue [might need to apply this to apo search too #TODO]')
args = parser.parse_args()



## User options
res_threshold = 3       # resolution cut-off for apo chains (angstrom), condition is '<='
NMR = 1                 # 0/1: discard/include NMR structures
lig_free_sites = 1      # 0/1: resulting apo sites will be free of any other known ligands in addition to specified ligands
water_as_ligand = 0     # 0/1: consider HOH atoms as ligands (can be used in combination with lig_free_sites)(strict)
save_session = 1        # 0/1: save each result as a PyMOL ".pse" session (zipped, includes annotations -recommended)
multisave = 0           # 0/1: save each result in a .pdb file (unzipped, no annotations -not recommended)
save_separate = 1       # 0/1: save each chain object in a separate file
look_in_archive = 0     # 0/1: search if the same query has been processed in the past (can give very fast results)

## Internal variables
#job_id = '0007'
overlap_threshold = 0  # % of overlap between apo and holo chain (w UniProt numbering), condition is ">="
ligand_scan_radius = '5' # angstrom radius to look around holo ligand(s) superposition
apo_chain_limit = 999    # limit number of apo chains to consider when aligning (for fast test runs)
min_tmscore = 0.5        # minimum acceptable TM score for apo-holo alignments (condition is "<" than)
beyond_hetatm = 0        # 0/1: when enabled, does not limit holo ligand detection to HETATM records for specified ligand/residue [might need to apply this to apo search too #TODO]
nonstd_rsds_as_lig = 0   # 0/1: ignore/consider non-standard residues as ligands
d_aa_as_lig = 0          # 0/1: ignore/consider D-amino acids as ligands

# Pass settings to a string
settings_param = 'res' + str(res_threshold) + '_NMR' + str(NMR) + '_ligfree' + str(lig_free_sites) + '_h2olig' + str(water_as_ligand) + '_overlap' + str(overlap_threshold) + '_ligrad' + str(ligand_scan_radius) + '_tmscore' + str(min_tmscore) + '_beyondhet' + str(beyond_hetatm) + '_nonstdrsds' + str(nonstd_rsds_as_lig) + '_drsds' + str(d_aa_as_lig)

# 3-letter names of amino acids and h2o (inverted selection defines ligands)
nolig_resn = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
if water_as_ligand == 0:    nolig_resn.append('HOH')
# Non-standard residues [SEP TPO PSU MSE MSO][1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG]
nonstd_rsds = "SEP TPO PSU MSE MSO 1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG PYG PYL SEC PHA".split()
if nonstd_rsds_as_lig == 0:    nolig_resn.extend(nonstd_rsds)
# D-amino acids
d_aminoacids = "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA".split()
if d_aa_as_lig == 0:    nolig_resn.extend(d_aminoacids)



## Define functions
def root_path():
    npath = os.path.normpath(os.getcwd())   # Normalize the path string for the OS
    path0 = os.path.join(npath.split(os.sep)[0], '\\', npath.split(os.sep)[1], npath.split(os.sep)[2])
    if os.path.exists(path0):        memo = "Root path found >> " + path0
    else:        memo = 'Error finding root path in working dir:' + npath        
    return path0
    print(memo)
def download_mmCIF_gz2(pdb_id, destination_path):   # Version 2 of download mmCIF gz (without exception handling)
    urlA = 'https://files.rcsb.org/download/'
    urlB = '.cif.gz'
    url = urlA + pdb_id.upper() + urlB
    file_path = destination_path + '\\' + pdb_id + urlB
    if not os.path.isfile(file_path):
        wget.download(url, destination_path)
        print('Downloading: ', pdb_id + urlB)
        return file_path
    else:        return file_path
def download_mmCIF_lig(lig_id, destination_path):   # Download mmCIF for ligands (without exception handling)
    urlA = 'https://files.rcsb.org/ligands/view/'
    urlB = '.cif'
    url = urlA + lig_id.upper() + urlB
    file_path = destination_path + '\\' + lig_id + urlB
    if not os.path.isfile(file_path):
        wget.download(url, destination_path)
        print('Downloading: ', lig_id + urlB)
        return file_path
    else:        return file_path
def add_log(msg, log_file):     # Create error log
    msg = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + '\t' + msg
    with open(path_root + '\\' + log_file, 'a') as file:
        file.write(msg + '\n')
script_name = os.path.basename(__file__)    #log_file = script_name[:-3] + '_rejected_res_' + infile1[:-4] + '.log'
log_file_dnld = script_name + '_downloadErrors.log' #log_file_dnld = job_id + '_' + script_name + '_downloadErrors' + '.log'
def next_path(path_pattern):    # Create incrementing directory name for each job
    i = 1
    # First do an exponential search
    while os.path.exists(path_pattern % i):
        i = i * 2
    # Result lies somewhere in the interval (i/2..i]
    # We call this interval (a..b] and narrow it down until a + 1 = b
    a, b = (i // 2, i)
    while a + 1 < b:
        c = (a + b) // 2 # interval midpoint
        a, b = (c, b) if os.path.exists(path_pattern % c) else (a, c)
    return path_pattern % b
def search_query_history(new_query_name, past_queries_filename):    # Find past job under the same query name, if found, return that job id
    dict_q = dict()
    try:
        with open (pathQRS + '\\' + past_queries_filename, 'r') as in_q:
            for line in in_q:
                dict_q[line.split('-')[0]] = line.split('-')[1][:-1]
        if new_query_name in dict_q.keys():            return dict_q[new_query_name]
        else:            return 0
    except:        return 0

    
    
    

# Set directories
path_root = root_path() + r'\Documents\Bioinfo_local\Ions\datasets_local\APO_candidates\webserver'
#path_root = r'C:\Users\TopOffice\Documents\GitHub\workDir\apoholo_web'
pathSIFTS = path_root + r'\SIFTS'           # Pre compiled files with UniProt PDB mapping
pathSTRUCTS = path_root + r'\structures'    # Directory with ALL pdb structures (used for fetch/download)
pathLIGS = path_root + r'\ligands'    # Directory with ALL pdb ligands (used for fetch/download)
pathRSLTS = next_path(path_root + r'\results' + '\\job_%s')     #pathRSLTS = path_root + r'\results' + '\\' + 'job_' + str(job_id)
pathQRS = path_root + r'\queries'             # Directory/index with parameters of previously run jobs

job_id = os.path.basename(os.path.normpath(pathRSLTS))

# Create directories if they don't exist
print('Setting up directories')
if os.path.isdir(pathSTRUCTS):    print('Structure directory:\t', pathSTRUCTS)
else:
    print('Creating structure directory:\t', pathSTRUCTS)
    os.makedirs(pathSTRUCTS)
if os.path.isdir(pathLIGS):    print('Ligands directory:\t', pathLIGS)
else:
    print('Creating ligands directory:\t', pathLIGS)
    os.makedirs(pathLIGS)
if save_separate == 1 or multisave == 1 or save_session == 1:
    if os.path.isdir(pathRSLTS):    print('Results directory:\t', pathRSLTS)
    else:
        print('Results directory:\t', pathRSLTS)
    os.makedirs(pathRSLTS)
if os.path.isdir(pathQRS):    print('Queries directory:\t', pathQRS)
else:
    print('Creating queries directory:\t', pathQRS)
    os.makedirs(pathQRS)
print('Done\n')

# Declare and load SIFTS input file(s)
#pathSIFTS = path0 + r'\ownCloud\Bioinfo_ownCloud\Projects\Ions\Uniprot_PDBchain\autodownload'
fileSIFTSdict = pathSIFTS + '\\' + "pdb_chain_uniprot_dict.txt" # downloaded SIFTS file unedited
fileRSIFTS = pathSIFTS + '\\' + "pdb_chain_uniprot_REVERSE_SPnum.txt" # pre-compiled rSIFTS file (reverse_SIFTS_SPnum.py)

print('Loading SIFTS dictionary')   # Load normal SIFTS dictionary as dict_SIFTS
with open (fileSIFTSdict, 'r') as input1:
    data = input1.read()
dict_SIFTS = ast.literal_eval(data)

print('Loading reverse SIFTS dictionary')   # Load reverse_SIFTS (SPnum) dictionary as dict_rSIFTS
with open (fileRSIFTS, 'r') as input2:
    data = input2.read()
dict_rSIFTS = ast.literal_eval(data)
print('Done\n')


## Parse single line input (line by line mode, 1 holo structure per line)
# if no chains specified, consider all chains #TODO
print('Parsing input')
struct = single_line_input.split()[0].lower()       # adjust case, struct = lower
user_chains = single_line_input.split()[1].upper()  # adjust case, chains = upper
ligand_names = single_line_input.split()[2].upper() # adjust case, ligands = upper
#user_position = single_line_input.split()[3] #TODO

# Parse chains
if not user_chains == 'ALL':
    user_chains = ''.join(user_chains)
    user_chains = user_chains.split(',')
    user_chains_bundle = '+'.join(user_chains)
    # Convert chains to structchain combos
    user_structchains = list()
    for user_chain in user_chains:
        user_structchain = struct.lower() + user_chain.upper()
        user_structchains.append(user_structchain)

# Parse ligands
ligand_names = ''.join(ligand_names)
ligand_names = ligand_names.split(',')
ligand_names_bundle = '+'.join(ligand_names)

# Print input info
print('Input structure:\t', struct)
if user_chains == 'ALL':
    print('Input chains:\t\t', user_chains)
else:
    print('Input chains:\t\t', user_chains)#, '\t', user_chains_bundle)
    print('Input structchains:\t', user_structchains)
print('Input ligands:\t\t', ligand_names)#, '\t', ligand_names_bundle)
#print('PyMOL version: ', cmd.get_version())
print('Done\n')


## Download or get path to query structure & ligand
try:
    struct_path = download_mmCIF_gz2(struct, pathSTRUCTS)
    print('Loading structure:\t', struct_path, '\n')
except:
    print('Error downloading structure:\t', struct, '\n')
print('Verifying ligands:\t', ligand_names)
for lig_id in ligand_names:
    try:
        lig_path = download_mmCIF_lig(lig_id, pathLIGS)
        with open (lig_path, 'r') as in_lig:
            for line in in_lig:
                if line.startswith('_chem_comp.name'):
                    lig_name = line.split()[1:]
                if line.startswith('_chem_comp.pdbx_synonyms'):
                    lig_syn = line.split()[1:]
                    print(lig_id, ' '.join(lig_name), ' '.join(lig_syn))
                    break
    except:        print('Error verifying ligand:\t', lig_id)



## Find Apo candidates (rSIFTS)

# Find & VERIFY input chains by UniProt ID (if they don't exist in uniprot, we cannot process them)
print('\nFinding & verifying query chains "', user_chains, '" by UniProt ID')
discarded_chains = list()   # Discarded chains (format: structchain  discard_msg)

if user_chains == 'ALL':
    user_chains = list()
    user_structchains = list()
    #print('Considering all chains in query structure, finding chains..')
    for key in dict_SIFTS:
        if key[:4] == struct:
            user_chains.append(key[4:])
            user_structchains.append(key)
            print(key, dict_SIFTS[key])
else:
    for user_structchain in user_structchains:
        try:
            print(user_structchain, dict_SIFTS[user_structchain])
        except:
            #print('User specified chain does not exist in UniProt, removing it from input:\t', user_structchain)
            user_chains.remove(user_structchain[4:])
            user_structchains.remove(user_structchain)
            discarded_chains.append(user_structchain + '\t' + 'No assigned UniProt ID\n')
user_chains_bundle = '+'.join(user_chains)
print('Input chains verified:\t', user_structchains, user_chains)
if len(discarded_chains) > 0:    print('Input chains rejected:\t', discarded_chains)






## Look up query in history, if found, return job path and end script
user_input_parameters = struct + '_' + ','.join(user_chains) + '_' + ','.join(ligand_names)
query_full = user_input_parameters + '_' + settings_param + '-' + job_id
#print(query_full)
if look_in_archive == 1:
    print('\nLooking for the same job in history')
    old_same_job = search_query_history(query_full.split('-')[0], 'queries.txt')
    # print('Result of lookup:', old_same_job)
    if old_same_job == 0:
        print('Same job not found, continuing process\n')
    else:
        if os.path.isdir(os.path.dirname(pathRSLTS) + '\\' + old_same_job):
            print('Same job found in history, printing path of old job results: ', os.path.dirname(pathRSLTS) + '\\' + old_same_job, '\n')
            print('Printing alignments of old job')
            with open (os.path.dirname(pathRSLTS) + '\\' + old_same_job + '\\' + 'results.csv', 'r') as old_in:
                for line in old_in:
                    print(line[:-1])
            print('\nDeleting new job folder', job_id)
            if os.path.isdir(pathRSLTS):
                os.rmdir(pathRSLTS)
            print('Done\nExiting')
            sys.exit(0)
        print('Old job directory not found, continuing process\n')
    
    
    



# Get apo candidates from rSIFTS dict
print('\nLooking for Apo candidates')
dictApoCandidates = dict()
#positive_overlap = dict()
#negative_overlap = dict()
uniprot_overlap = dict()

for user_structchain in user_structchains:
    for key, values in dict_rSIFTS.items(): # iterate over reverse SIFTS chains/uniprot IDs
        for i in values: # iterate over the values in each key, i = struct/chain combo
            if i.split()[0] == user_structchain:
                #uniprot_id = key # not needed
                structchain = i.split()[0] # pass structchain to variable
                x1 = i.split()[1]    # holo SP BEG
                x2 = i.split()[2]    # holo SP END
                
                # Find apo candidates by filtering by coverage
                for candidate in values:

                    # Discard same structure and/or chain with input (holo)
                    if candidate.split()[0][:4] != structchain[:4]: # discard any structure same with query (stricter)
                    #if candidate.split(' ')[0] != structchain: # only discard same chain
                    
                        y1 = candidate.split()[1]
                        y2 = candidate.split()[2]                        
                        
                        # Calculate overlap and % of overlap on query (x1,x2)
                        result =  min(int(x2), int(y2)) - max(int(x1), int(y1))
                        percent = result / (int(x2) - int(x1)) * 100
                        percent = round(percent, 1)     # round the float
                                                
                        # Only consider positive overlap (negative overlap may occur cause of numbering)
                        #uniprot_overlap.setdefault(i.split()[0], []).append(candidate.split()[0]+' '+str(percent))
                        uniprot_overlap.setdefault(candidate.split()[0], []).append(i.split()[0] + ' ' + str(percent))
                        if overlap_threshold != 0 and percent >= overlap_threshold or overlap_threshold == 0:
                            dictApoCandidates.setdefault(i, []).append(candidate+' '+str(result)+' '+str(percent))
                        

#print('Candidates with positive overlap: ', positive_overlap)
#print('Candidates with negative overlap: ', negative_overlap)
print('Candidate chains over user-specified overlap threshold [', overlap_threshold, '%]: ', sum([len(dictApoCandidates[x]) for x in dictApoCandidates if isinstance(dictApoCandidates[x], list)]))
print('')



## Apo candidates evaluation

# Put all structures for downloading into set
apo_candidate_structs = set()
for key, values in dictApoCandidates.items():
    for i in values:    # iterate over the values in each key, i = struct/chain combo & SP_BEG SP_END etc
        struct = i.split()[0][:4] # split value strings to get structure only
        apo_candidate_structs.add(struct)
print('Total structures to download for parsing: ', len(apo_candidate_structs), '\n')

# Download/fetch/load the Apo candidate structures to specified directory [this should be replaced by fetch later #TODO]
for apo_candidate_structure in apo_candidate_structs:
    try:
        structPath = download_mmCIF_gz2(apo_candidate_structure, pathSTRUCTS)
    except Exception as ex1:
        template = "Exception {0} occurred. \n\t\t\t\t\tArguments:{1!r}"
        message = template.format(type(ex1).__name__, ex1.args) + apo_candidate_structure
        add_log(message, log_file_dnld)
        print(f'*apo file {apo_candidate_structure} not found')
        
'''# Add extra structures for testing [NMR struct: 1hko | cryo-em struct: 6nt5]
apo_candidate_structs.add('6nt5') #EM
apo_candidate_structs.add('1hko') #NMR'''

# Parse (mmCIF) structures to get resolution & method. Apply cut-offs
print('Checking resolution and experimental method of Apo candidate structures')
for apo_candidate_struct in apo_candidate_structs:
    apo_candidate_structPath = download_mmCIF_gz2(apo_candidate_struct, pathSTRUCTS)
    resolution = '?'
    with gzip.open (apo_candidate_structPath, 'rt') as mmCIFin:
        for line in mmCIFin:
            try:
                if line.split()[0] == '_exptl.method':
                    method = line.split("'")[1] # capture experimental method #method = ' '.join(line.split()[1:])
                    if method == 'SOLUTION NMR': # break fast if method is 'NMR'
                        break
                elif line.split()[0] == '_refine.ls_d_res_high' and float(line.split()[1]):
                    resolution = float(line.split()[1]) # X-ray highest resolution
                    break
                elif line.split()[0] == '_em_3d_reconstruction.resolution' and float(line.split()[1]):
                    resolution = float(line.split()[1]) # EM resolution
                    break
            except:# Exception as ex: # getting weird but harmless exceptions
                print('Problem parsing structure: ', apo_candidate_struct)#, ex)
        try:
            if NMR == 1 and method == 'SOLUTION NMR' or resolution <= res_threshold:
                print(apo_candidate_struct, ' resolution:\t', resolution, '\t', method, '\tPASS') # Xray/EM
            else:
                discarded_chains.append(apo_candidate_struct + '\t' + 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']\n')
                print(apo_candidate_struct, ' resolution:\t', resolution, '\t', method, '\tFAIL') # Xray/EM
        except:
            discarded_chains.append(apo_candidate_struct + '\t' + 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']\n')
            print(apo_candidate_struct, ' resolution:\t', resolution, '\t\t', method, '\t\tFAIL') # NMR
print('Done\n')


# Compile set of structures to be discarded (non verbose)
discard_structs = set()
for i in discarded_chains:
    if len(i.split()[0]) == 4:
        discard_structs.add(i.split()[0])

# Discard apo entries below threshold(s). Put remainder into new dict, discard UniProt numbering
print('Structures to discard:\t', len(discard_structs))
if len(discard_structs) > 0:    print('Discarding structures\t', discard_structs)
dictApoCandidates_1 = dict()
for key, values in dictApoCandidates.items():
    for i in values[:apo_chain_limit]:            
        if i.split()[0][:4] in discard_structs:
            print('removing apo', i.split()[0], 'from holo', key.split()[0])
            #pass
        else:
            dictApoCandidates_1.setdefault(key.split()[0], []).append(i.split()[0])
print('\nApo candidate chains satisfying user requirements (method/resolution) [', res_threshold, 'Å ]: ', sum([len(dictApoCandidates_1[x]) for x in dictApoCandidates_1 if isinstance(dictApoCandidates_1[x], list)]))


# Open apo winner structures, align to holo, and check if the superimposed (ligand) sites are ligand-free
if beyond_hetatm == 1:    search_name = 'resn '
else:    search_name = 'hetatm and resn '

holo_lig_positions = dict()
apo_holo_dict = dict()
for holo_structchain, apo_structchains in dictApoCandidates_1.items():
    print('')
    holo_struct = holo_structchain[:4]
    holo_chain = holo_structchain[4:]
    holo_struct_path = download_mmCIF_gz2(holo_struct, pathSTRUCTS)
    
    # Initialize PyMOL but don't reload Holo if present
    if holo_struct in cmd.get_object_list('all'): # object names
    #if holo_struct in cmd.get_names('all'): # object and selection names
        cmd.delete('not ' + holo_struct)
    else:
        cmd.reinitialize('everything')
        cmd.load(holo_struct_path)
    
    cmd.select(holo_struct + holo_chain, holo_struct + '& chain ' + holo_chain) 
    if len(dictApoCandidates_1) > 0 and save_separate == 1 and not os.path.isfile(pathRSLTS + '\\holo_' + holo_struct + '.cif.gz'):
        cmd.save(pathRSLTS + '\\holo_' + holo_struct + '.cif.gz', holo_struct)
    
    # Find & name specified ligands
    #ligands_selection = cmd.select('query_ligands', 'hetatm and resn ' + ligand_names_bundle + ' and chain ' + holo_chain) # resn<->name
    ligands_selection = cmd.select('query_ligands', search_name + ligand_names_bundle + ' and chain ' + holo_chain) # resn<->name
    if ligands_selection == 0:
        print('No ligands found in author chain, trying PDB chain')
        ligands_selection = cmd.select('query_ligands', search_name + ligand_names_bundle + ' and segi ' + holo_chain)
        if ligands_selection == 0:
            print('No ligands found in PDB chain, skipping: ', holo_structchain)
            continue
    
    ligands_atoms = cmd.identify('query_ligands', mode=0)
    print('Query ligand selection atoms:\t', ligands_atoms, holo_chain)    
    
    # Get positions of specified ligands (better than atom ids to specify during alignment)
    myspace = {'positions': []} # temporary dict with fixed key name
    for atom in ligands_atoms: #this is a bulk of atoms for all ligands, many atoms can belong to a single ligand
        cmd.iterate('id ' + str(atom), 'positions.append(resi +" "+ chain +" "+ resn)', space = myspace)

    # Transfer temporary list with positions to dict
    for key, values in myspace.items():
      for i in values:
          holo_lig_positions.setdefault(holo_structchain , []).append(i)
    print('Holo ligands positions for chain: ', holo_structchain,  holo_lig_positions.get(holo_structchain))
    
    # Name holo ligands as PyMOL selections. Put real (detected) ligand names into set
    holo_lig_names = set()
    #holo_lig_names.update(ligand_names)
    for ligand in holo_lig_positions[holo_structchain]:
        resi = ligand.split()[0]
        chain = ligand.split()[1]
        resn = ligand.split()[2]
        ligand_ = ligand.replace(' ', '_')
        holo_lig_names.add(resn)
        s1 = cmd.select('holo_' + ligand_, 'model ' + holo_struct + '& resi ' + resi + '& chain ' + chain + '& resn ' + resn)

    
    # Start Apo chain loop. Align and mark atom selections around holo ligand
    for apo_structchain in apo_structchains:
        apo_struct = apo_structchain[:4]
        apo_chain = apo_structchain[4:]
        apo_struct_path = download_mmCIF_gz2(apo_struct, pathSTRUCTS)
        
        if apo_struct in cmd.get_object_list('all'):
            pass
        else:    
            cmd.load(apo_struct_path)
        cmd.select(apo_struct + apo_chain, apo_struct + '& chain ' + apo_chain)
        
        # Align apo-holo
        try:
            aln_rms = cmd.align(apo_struct + '& chain ' + apo_chain, holo_struct + '& chain ' + holo_chain, cutoff=10.0, cycles=1)
            aln_tm = psico.fitting.tmalign(apo_struct + '& chain ' + apo_chain, holo_struct + '& chain ' + holo_chain, quiet=1, transform=0)
            print('Alignment RMSD/TM score:', apo_structchain, holo_structchain, round(aln_rms[0], 3), aln_tm)
        except:
            discarded_chains.append(apo_structchain + '\t' + 'Alignment error\n')
            print('Alignment RMSD/TM score: ERROR')
            print('*poor alignment, discarding chain ', apo_structchain)
            continue
        
        # Discard bad alignments
        if aln_tm < min_tmscore:
            discarded_chains.append(apo_structchain + '\t' + 'Poor alignment [RMSD/TM]: ' + str(round(aln_rms[0], 3)) +'/'+ str(aln_tm) + '\n')
            print('*poor alignment, discarding chain ', apo_structchain)
            continue

        found_ligands = set()
        found_ligands_xtra = set()
        
        # Start ligand look-up loop for apo chain
        for ligand in holo_lig_positions[holo_structchain]:

            # Around selection [this is looking for ligands in every (valid) chain alignment, not just the standard locus of holo ligand(s)]
            ligand_ = ligand.replace(' ', '_') #remove spaces for selection name
            #s2 = cmd.select(apo_structchain + '_arnd_' + ligand_, 'model ' + apo_struct + '& chain ' + apo_chain + ' near_to ' + ligand_scan_radius + ' of holo_' + ligand_)
            s2 = cmd.select(apo_structchain + '_arnd_' + ligand_, 'model ' + apo_struct + '& hetatm & not solvent' + ' near_to ' + ligand_scan_radius + ' of holo_' + ligand_)
        
            # Put selected atoms in a list, check their name identifiers to see if holo ligand name is present
            myspace_a = {'a_positions': []}
            for a_atom in cmd.identify(apo_structchain + '_arnd_' + ligand_):
                cmd.iterate('id ' + str(a_atom), 'a_positions.append(resi +" "+ chain +" "+ resn)', space = myspace_a)
            
            # Transfer previous dict[key] values (just resn) into set (for easier handling)
            apo_lig_names = set()
            for a_position in myspace_a['a_positions']:
                a_atom_lig_name = a_position.split()[2]
                apo_lig_names.add(a_atom_lig_name)
            
            # Assess apo ligands
            for i in apo_lig_names:
                if i in holo_lig_names or i in ligand_names:    # check in both lists (detected holo ligs + query ligs)
                    found_ligands.add(i)    #break #print('Holo ligand found in Apo: ', apo_structchain, i)
                elif i not in nolig_resn:
                    found_ligands_xtra.add(i)
                    
        
        # Print verdict for chain & save it as a separate object/file
        print(f'*query ligands: {ligand_names}\tdetected ligands: {holo_lig_names}\t detected apo ligands: {apo_lig_names}\tfound query ligands: {found_ligands}\tfound non-query ligands: {found_ligands_xtra}')
        if lig_free_sites == 1 and len(found_ligands_xtra) == 0 and len(found_ligands) == 0 or lig_free_sites == 0 and len(found_ligands) == 0:
            apo_holo_dict.setdefault(holo_structchain , []).append(apo_structchain + ' ' + uniprot_overlap[apo_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)))
            print('PASS')   #print('*===> Apo chain', apo_structchain, ' clean of query ligands ', holo_lig_names)
            if save_separate ==1:                cmd.save(pathRSLTS + '\\' + apo_structchain + '_aln_to_' + holo_structchain + '.cif.gz', apo_structchain)
        else:
            print('FAIL')   #print('*apo chain', apo_structchain, ' includes query ligands ', found_ligands)
            
    
    # Clean objects/selections
    apo_win_structs = set()
    for key, values in apo_holo_dict.items():
        for value in values:
            apo_win_structs.add(value.split()[0][:4])
    add_sele = ['quer', 'holo']     # add first 4 chars of selection names we want to keep
    good_selections = apo_win_structs.copy()    # shallow copy
    good_selections.update(add_sele)
    all_selections = cmd.get_names('all')
    for i in all_selections:
        if i[:4] == holo_struct or i[:4] in good_selections:
            pass
        else:
            cmd.delete(i)
    all_selections = cmd.get_names('all')
    for i in all_selections:    # Delete 0 atom selections
        if cmd.count_atoms(i) == 0:
            cmd.delete(i)
    
    cmd.disable('all') # toggles off the display of all objects
    cmd.enable(holo_struct)
    cmd.deselect()
    cmd.reset()
    cmd.center('query_ligands')
    
    # Save results as session (.pse.gz) or multisave (.pdb)
    #filename_body = pathRSLTS + '\\' + 'aln_' + holo_structchain + '_to_' + '_'.join(cmd.get_object_list('all and not ' + holo_struct))
    filename_body = pathRSLTS + '\\' + 'aln_' + holo_struct
    filename_pse = filename_body + '.pse.gz'
    filename_multi = filename_body + '_multi.cif'    
    #apo_win_structs_filename = '_'.join(list(apo_win_structs))
    #filename_pse = pathRSLTS + '\\' + 'aln_' + holo_structchain + '_to_' + apo_win_structs_filename + '.pse.gz'

    #if len(cmd.get_object_list('all')) > 1:
    if len(dictApoCandidates_1) > 0:
        #if save_separate ==1:            
            #cmd.save(pathRSLTS + '\\holo_' + holo_structchain + '.pdb.gz', holo_structchain)
            #if not os.path.isfile(pathRSLTS + '\\holo_' + holo_struct + '.pdb.gz'):                cmd.save(pathRSLTS + '\\holo_' + holo_struct + '.pdb.gz', holo_struct)
        if save_session == 1:            cmd.save(filename_pse)
        if multisave == 1:            cmd.multisave(filename_multi, append=1)
        
    
print('')
if len(apo_holo_dict) > 0:
    if save_separate == 1 or multisave == 1 or save_session == 1:
        
        # Write dictionary to file
        header = "#HEADER: {holo_chain: [apo_chain %UniProt_overlap RMSD TM_score]\n"
        filename_aln = pathRSLTS + '\\aln_' + '_'.join(list(apo_holo_dict.keys()))
        with open (filename_aln + '.txt', 'wt') as out1:
            out1.write(header)
            out1.write(str(apo_holo_dict))
            
        # Write CSV file
        filename_csv = pathRSLTS + '\\results.csv'
        header = "#holo_chain,apo_chain,%UniProt_overlap,RMSD,TM_score\n"
        with open (filename_csv, 'w') as csv_out:
            csv_out.write(header)
            for key, values in apo_holo_dict.items():
                for value in values:
                    #csv_out.write("%s,%s,%s,%s,%s,%s\n"%(key,value,apo_holo_dict[key][0].split()[0],apo_holo_dict[key][0].split()[1],apo_holo_dict[key][0].split()[2],apo_holo_dict[key][0].split()[3]))
                    csv_out.write("%s,%s\n"%(key,','.join(value.split())))
    # Print dict
    print('Apo holo results: ')
    for key in apo_holo_dict: print(key, apo_holo_dict.get(key))
else:    print('No apo forms found')
print('\nDone')


# Append the name of the query and the job_id in the queries.txt
if job_id:
    with open (pathQRS + '\\' + 'queries.txt', 'a') as out_q:
        out_q.write(query_full + '\n')
    


    



