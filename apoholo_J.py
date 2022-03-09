# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:24:57 2021

@author: ChrisX
"""
# Apo-Holo Juxtaposition - AHoJ
from common import get_workdir, load_dict_binary, tmalign2

import __main__
__main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
# import pymol.cmd as cmd
import psico.fitting
import psico.fullinit
import pymol2

import ast
import gzip
import os
import wget
import time
import argparse
import sys
from dataclasses import dataclass

from concurrent.futures import ThreadPoolExecutor as PoolExecutor; import threading           # multi-threading
# from concurrent.futures import ProcessPoolExecutor as PoolExecutor; import multiprocessing  # multi-processing (doesn't work atm)


_global_lock = threading.Lock()                      # multi-threading
# global_lock = multiprocessing.Manager().Lock()     # multi-processing (must be moved to main)

'''
Given an experimental protein structure (PDB code), with optionally specified chain(s) and ligand(s), find its equivalent apo and holo forms.
The program will look for both apo and holo forms of the query structure. Structures are processed chain by chain.

The user can specify the following input arguments depending on the mode of search
i) When looking for apo from holo:
-Min arguments: PDB code
-Max arguments: PDB code, chain(s), ligand(s)
ii) When looking for holo from apo:
-Min arguments: PDB code
-Max arguments: PDB code, chain(s)
'''

# TODO add smart synching with PDB
# TODO adjust search radius according to mol. weight of ligand


##########################################################################################################
# Define functions
##########################################################################################################

def download_mmCIF_gz2(pdb_id, destination_path):   # Version 2 of download mmCIF gz (without exception handling)
    urlA = 'https://files.rcsb.org/download/'
    urlB = '.cif.gz'
    url = urlA + pdb_id.upper() + urlB
    file_path = destination_path + '/' + pdb_id + urlB
    if not os.path.isfile(file_path):
        wget.download(url, destination_path)
        print('Downloading: ', pdb_id + urlB)
        return file_path
    else:
        return file_path


def download_mmCIF_lig(lig_id, destination_path):   # Download mmCIF for ligands (without exception handling)
    urlA = 'https://files.rcsb.org/ligands/view/'
    urlB = '.cif'
    url = urlA + lig_id.upper() + urlB
    file_path = destination_path + '/' + lig_id + urlB
    if not os.path.isfile(file_path):
        wget.download(url, destination_path)
        print('Downloading: ', lig_id + urlB)
        return file_path
    else:
        return file_path


def add_log(msg, log_file):     # Create error log
    msg = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + '\t' + msg
    with open(log_file, 'a') as file:
        file.write(msg + '\n')


def next_job(path_pattern):    # Create incrementing directory name for each job
    i = 1
    while os.path.exists(path_pattern % i):     # First do an exponential search
        i = i * 2
    a, b = (i // 2, i)  # Result lies somewhere in the interval (i/2..i]    # We call this interval (a..b] and narrow it down until a + 1 = b
    while a + 1 < b:
        c = (a + b) // 2  # interval midpoint
        a, b = (c, b) if os.path.exists(path_pattern % c) else (a, c)
    return path_pattern % b


def search_query_history(pathQRS, new_query_name, past_queries_filename):    # Find past job under the same query name, if found, return the job id
    dict_q = dict()
    try:
        with open(pathQRS + '/' + past_queries_filename, 'r') as in_q:
            for line in in_q:
                dict_q[line.split('-')[0]] = line.split('-')[1][:-1]
        if new_query_name in dict_q.keys():
            return dict_q[new_query_name]
        else:
            return 0
    except:
        return 0


def wrong_input_error():  # arg_job_id, arg_pathRSLTS):
    print('ERROR: Wrong input format\nPlease use a whitespace character to separate input arguments')
    print('Input format: <pdb_id> <chains> <ligands> or <pdb_id> <chains> or <pdb_id> <ligands> or <pdb_id>')
    print('Input examples: "3fav A,B ZN" or "3fav ZN" or "3fav ALL ZN" or "3fav"')
    #print('Exiting & deleting new results folder:', path_job_results)
    #if os.path.isdir(path_job_results):
    #    os.rmdir(path_job_results)
    sys.exit(1)  # exit with error

    
# Join list of ligand codes to a single string, sorted
# separator: '-'
# if the list is empty return '-'
def join_ligands(ligands):
    if not ligands:
        return '-'
    return '-'.join(sorted(ligands))


##########################################################################################################

@dataclass
class Query:
    struct: str
    chains: str           # maybe even change to list
    ligands: str          # maybe even change to list
    position: str
    autodetect_lig: bool
    water_as_ligand: bool


@dataclass
class QueryResult:
    """ Result of a single query """
    result_dir: str
    num_apo_chains: int = 0
    num_holo_chains: int = 0
    error: str = None

    def __str__(self):
        return f"QueryResult[ apo: {self.num_apo_chains}, holo: {self.num_holo_chains}, error: {self.error}, dir: {self.result_dir})]"


@dataclass
class PrecompiledData:
    """
    Pre-compiled data needed for computation
    UniProt PDB mapping
    """
    dict_SIFTS: dict   # regular SIFTS dictionary
    dict_rSIFTS: dict  # reverse SIFTS (SPnum) dictionary


def verify_ligands(ligand_names, pathLIGS):
    #if autodetect_lig == 0 or ligand_names is not None:
    #print('Verifying ligands:\t', ligand_names)
    for lig_id in ligand_names:
        #try:
        lig_path = download_mmCIF_lig(lig_id, pathLIGS)
        with open(lig_path, 'r') as in_lig:
            for line in in_lig:
                if line.startswith('_chem_comp.name'):
                    lig_name = line.split()[1:]
                if line.startswith('_chem_comp.pdbx_synonyms'):
                    lig_syn = line.split()[1:]
                    #print('>', lig_id, ' '.join(lig_name), ' '.join(lig_syn))
                    print('Verifying ligands:\t', ligand_names, ' > ', lig_id, ' '.join(lig_name), ' '.join(lig_syn))
                    break
        #except:
            #print('Error verifying ligand:\t', lig_id)


def parse_query(query: str, autodetect_lig: bool = False, water_as_ligand: bool = False) -> Query:

    # Parse single line input (line by line mode, 1 holo structure per line)
    # if no chains specified, consider all chains
    print(f"Parsing query '{query}'")
    query = query.strip()
    parts = query.split()

    struct = parts[0].lower()
    chains = 'ALL'
    ligands = None
    position = None

    # Define non-ligands (3-letter names of amino acids and h2o)
    std_rsds = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    nonstd_rsds = "SEP TPO PSU MSE MSO 1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG PYG PYL SEC PHA HOH".split()
    d_aminoacids = "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA".split()
    nolig_resn = list()
    nolig_resn.extend(std_rsds + nonstd_rsds + d_aminoacids)


    if len(struct) != 4:
        raise ValueError(f"Invalid query '{query}': '{struct}' is not a valid PDB structure code")

    if len(parts) == 1:
        autodetect_lig = 1               # overrides cmd line param
    elif len(parts) == 2:  # and autodetect_lig == 1:
        chains = parts[1].upper()
        autodetect_lig = 1
    #elif len(parts) == 2 and autodetect_lig == 0:  # this triggers "ALL" chains mode
        #ligands = parts[1].upper()
    elif len(parts) == 3:
        chains = parts[1].upper()        # adjust case, chains = upper
        ligands = parts[2].upper()       # adjust case, ligands = upper
    
    # When position is specified, there has to be a single ligand/residue specified
    elif len(parts) == 4 and len(parts[2]) < 4 and len(parts[2].split(',')) == 1:# and int(parts[3]):
        try:
            chains = parts[1].upper()
            ligands = parts[2].upper()
            position = str(int(parts[3]))   # test if int
            if ligands in std_rsds:
                autodetect_lig = 1
                #print('\nLigand is standard residue')
        except:
            raise ValueError(f"Invalid query '{query}': wrong number of parts")
    else:
        raise ValueError(f"Invalid query '{query}': wrong number of parts")

    if chains == '*' or chains == '?':
        chains = 'ALL'
    elif not all(chain.isalnum() for chain in chains.split(',')):
        raise ValueError(f"Invalid query '{query}': only alphanumeric characters allowed as chains")
    if ligands == '*' or ligands == '?':
        ligands = None
        autodetect_lig = 1

    # Remove star from ligands str
    if ligands is not None:
        if ',*,' in ligands:
            autodetect_lig = 1
            ligands = ligands.replace(",*,", ",")
        elif ',*' in ligands:
            autodetect_lig = 1
            ligands = ligands.replace(",*", "")
        elif '*,' in ligands:
            autodetect_lig = 1
            ligands = ligands.replace("*,", "")
        elif '*' in ligands and len(ligands) > 1:
            autodetect_lig = 1
            ligands = ligands.replace("*", "")

    # If ligand is HOH or std residue or non-std residue,  make sure that:
    # i) there is just one specified (handled earlier)
    # ii) there is a fourth argument (position)
    if ligands == 'HOH' and position is not None:
        water_as_ligand = 1

    for i in nolig_resn:
        if ligands == i and position is None:
            raise ValueError(f"Invalid query '{query}': specify index position of HOH or residue") 
    
    return Query(struct=struct, chains=chains, ligands=ligands, position=position, autodetect_lig=autodetect_lig, water_as_ligand=water_as_ligand)


def load_precompiled_data_txt(workdir) -> PrecompiledData:
    pathSIFTS = workdir + '/SIFTS'
    fileSIFTSdict = pathSIFTS + '/' + "pdb_chain_uniprot_dict.txt"        # regular SIFTS file stripped
    fileRSIFTS = pathSIFTS + '/' + "pdb_chain_uniprot_REVERSE_SPnum.txt"  # pre-compiled rSIFTS file (reverse_SIFTS_SPnum.py)

    print('Loading SIFTS dictionary')   # Load normal SIFTS dictionary as dict_SIFTS
    with open(fileSIFTSdict, 'r') as input1:
        data = input1.read()
    dict_SIFTS = ast.literal_eval(data)

    print('Loading reverse SIFTS dictionary')   # Load reverse_SIFTS (SPnum) dictionary as dict_rSIFTS
    with open(fileRSIFTS, 'r') as input2:
        data = input2.read()
    dict_rSIFTS = ast.literal_eval(data)

    return PrecompiledData(dict_SIFTS=dict_SIFTS, dict_rSIFTS=dict_rSIFTS)


def load_precompiled_data_bin(workdir) -> PrecompiledData:
    pathSIFTS = workdir + '/SIFTS'
    fileSIFTSdict = pathSIFTS + '/pdb_chain_uniprot_dict.bin'
    fileRSIFTS = pathSIFTS + '/pdb_chain_uniprot_REVERSE_SPnum.bin'

    print('Loading SIFTS dictionary')
    dict_SIFTS = load_dict_binary(fileSIFTSdict)

    print('Loading reverse SIFTS dictionary')
    dict_rSIFTS = load_dict_binary(fileRSIFTS)

    return PrecompiledData(dict_SIFTS=dict_SIFTS, dict_rSIFTS=dict_rSIFTS)


def load_precompiled_data(workdir) -> PrecompiledData:
    """Load pre-compiled data generated by prepare.py"""
    res = load_precompiled_data_bin(workdir)
    print('Done loading pre-compiled data\n')
    return res


def process_query(query, workdir, args, data: PrecompiledData = None) -> QueryResult:
    """
    Process single line query
    :param query: single line query with format "<pdb_id> <chains> <ligands>" (see README.md)
    :param workdir: global work directory
    :param args: all parsed cmd line args
    :param data: if not provided will be loaded on the fly
    :return:
    """
    
    ''' Test input (overrides argparse) '''
    #multiline_input = '3fav all zn\n1a73 a zn,MG,HEM\n5ok3 all tpo'
    #query = '1a0u' #hem, big search
    #query = '1a73 a zn'#',MG,HEM'
    #query = '5ok3 all tpo' #phosphothreonine, no apos
    #query = '2ZB1 all gk4'
    #query = '7l1f all F86' # too long
    #query = '1SI4 cyn'
    #query = '2v7c a'
    #query = '5gss all gsh' # slow
    #query = '1jq8 so4'
    #query = '1l5h b CLF'
    #query = '1DB1 vdx' #vitamin D3 study
    #query = '3IXJ all 586' # beta-secretase 1 with inhibitor (cryptic?) # too long
    #query = '2jds all L20' # cAMP-dependent protein kinase w inhibitor #202 chains 145 structs, long
    #query = '1pzo all cbt' # TEM-1 Beta-Lactamase with Core-Disrupting Inhibitor #115 chains, 58 structs, longish
    #query = '1qsh d heg' # long
    
    #query = '3N3I all roc' # HIV mutations
    #query = '2whh all ppn' # related to upper example, multiple identical structchains in values of dict
    
    # Fast examples
    #query = '2v0v' # Fully apo structure
    #query = '3CQV all hem'#,coh'# hem,f86,mg,tpo,act,jkl,ue7,909' # apohaemoglobin study [OK]
    #query = '3fav all zn' # [OK]
    #query = '1py2 d frh' # 228 chains, 180 valid, long - run only on one chain [OK*]
    #query = '2hka all c3s' # bovine NPC2 complex with cholesterol sulfate [OK]
    #query = '2v57 a,c prl' # apo-holo SS changes in TetR-like transcriptional regulator LfrR in complex with proflavine [OK]


    # Init independent pymol instance
    pm = pymol2.PyMOL()
    pm.start()
    cmd = pm.cmd

    # Basic
    res_threshold = args.res_threshold
    NMR = args.NMR
    xray_only = args.xray_only
    lig_free_sites = args.lig_free_sites
    autodetect_lig = args.autodetect_lig
    reverse_search = args.reverse_search

    # Advanced
    overlap_threshold = args.overlap_threshold
    lig_scan_radius = args.lig_scan_radius
    min_tmscore = args.min_tmscore

    # Experimental
    water_as_ligand = args.water_as_ligand
    nonstd_rsds_as_lig = args.nonstd_rsds_as_lig
    d_aa_as_lig = args.d_aa_as_lig
    beyond_hetatm = args.beyond_hetatm
    look_in_archive = args.look_in_archive

    # Internal
    apo_chain_limit = args.apo_chain_limit

    # Saving
    save_oppst = args.save_oppst
    save_separate = args.save_separate
    save_session = args.save_session
    multisave = args.multisave

    # Adjust input, resolve conflicts
    if reverse_search == 1:
        autodetect_lig = 1
    lig_scan_radius = str(lig_scan_radius)  # needs to be str
    reverse_mode = False

    # Pass settings to a string
    settings_str = 'res' + str(res_threshold) + '_NMR' + str(NMR) + '_ligfree' + str(lig_free_sites) + '_h2olig' + str(water_as_ligand) + '_overlap' + str(overlap_threshold) + '_ligrad' + str(lig_scan_radius) + '_tmscore' + str(min_tmscore) + '_beyondhet' + str(beyond_hetatm) + '_nonstdrsds' + str(nonstd_rsds_as_lig) + '_drsds' + str(d_aa_as_lig)


    # Define non-ligands (3-letter names of amino acids and h2o)
    nolig_resn = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    std_rsds = list(nolig_resn)
    if water_as_ligand == 0:
        nolig_resn.append('HOH')
    # Non-standard residues [SEP TPO PSU MSE MSO][1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG]
    nonstd_rsds = "SEP TPO PSU MSE MSO 1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG PYG PYL SEC PHA".split()
    if nonstd_rsds_as_lig == 0:
        nolig_resn.extend(nonstd_rsds)
    # D-amino acids
    d_aminoacids = "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA".split()
    if d_aa_as_lig == 0:
        nolig_resn.extend(d_aminoacids)

    # Set directories, create job_id
    path_root = workdir
    pathSTRUCTS = path_root + '/structures'    # Directory with ALL pdb structures (used for fetch/download)
    pathLIGS = path_root + '/ligands'          # Directory with ALL pdb ligands (used for fetch/download)
    pathQRS = path_root + '/queries'           # Directory/index with parameters of previously run jobs


    # TODO make next job_id generation less clumsy
    global _global_lock
    if _global_lock is None:
        _global_lock = threading.Lock()   # to allow unit tests
    with _global_lock:
        generated_path_job_results = next_job(path_root + '/results/job_%s')     #pathRSLTS = path_root + r'/results' + '/' + 'job_' + str(job_id)
        if args.out_dir is not None:
            path_results = args.out_dir  # user defined
        else:
            path_results = generated_path_job_results   # generated
        if not os.path.isdir(path_results):
            os.makedirs(path_results)                     # must be created inside lock to ensure each next_job is unique
    print('Results directory:\t', path_results)
    job_id = os.path.basename(os.path.normpath(path_results))


    if data is None:
        data = load_precompiled_data(workdir)
    dict_SIFTS = data.dict_SIFTS
    dict_rSIFTS = data.dict_rSIFTS

    # Get additional info
    # script_name = os.path.basename(__file__)    #log_file = script_name[:-3] + '_rejected_res_' + infile1[:-4] + '.log'
    # log_file_dnld = script_name + '_downloadErrors.log' #log_file_dnld = job_id + '_' + script_name + '_downloadErrors' + '.log'
    #log_file_dnld = path_root + '/download_errors.log'

    print('PyMOL version: ', cmd.get_version())

    # Create directories if they don't exist
    print('Setting up directories')

    if os.path.isdir(pathSTRUCTS):
        print('Structure directory:\t', pathSTRUCTS)
    else:
        print('Creating structure directory:\t', pathSTRUCTS)
        os.makedirs(pathSTRUCTS)
    if os.path.isdir(pathLIGS):
        print('Ligands directory:\t', pathLIGS)
    else:
        print('Creating ligands directory:\t', pathLIGS)
        os.makedirs(pathLIGS)
    #if save_separate == 1 or multisave == 1 or save_session == 1: # bypass

    if os.path.isdir(pathQRS):
        print('Queries directory:\t', pathQRS)
    else:
        print('Creating queries directory:\t', pathQRS)
        os.makedirs(pathQRS)
    print('Done\n')



    # Parse single line input (line by line mode, 1 holo structure per line)
    # if no chains specified, consider all chains
    print('Parsing input')


    try:
        q = parse_query(query, autodetect_lig, water_as_ligand)
    except ValueError as e:
        print(e)
        wrong_input_error()

    user_chains = q.chains
    struct = q.struct
    ligand_names = q.ligands
    position = q.position
    autodetect_lig = q.autodetect_lig
    water_as_ligand = q.water_as_ligand

    # Parse chains
    if not user_chains == 'ALL':
        user_chains = ''.join(user_chains)
        user_chains = user_chains.split(',')
        #user_chains_bundle = '+'.join(user_chains)
        # Convert chains to structchain combos
        user_structchains = list()
        for user_chain in user_chains:
            user_structchain = struct.lower() + user_chain.upper()
            user_structchains.append(user_structchain)

    #if ligand_names is None:  # This should be safe to remove as well
    #    print("Input ligands were not defined!")
        #ligand_names = 'autodetect'
        # sys.exit(1) ?
    
    # Verify input structure here
    try:
        print('Verifying structure:', struct)
        download_mmCIF_gz2(struct, pathSTRUCTS)
    except:
        raise ValueError(f"Invalid PDB ID '{query}': use a valid 4-letter PDB code") 
        
    # Verify ligands here
    if autodetect_lig == 0 or ligand_names is not None:
        try:
            verify_ligands(ligand_names.split(','), pathLIGS)
        except:
            raise ValueError(f"Invalid ligands in query '{query}': use PDB ligand names")

    # Parse ligands
    if autodetect_lig == 1 and ligand_names is not None or autodetect_lig == 0:
        ligand_names = ''.join(ligand_names)
        ligand_names = ligand_names.split(',')
        ligand_names_bundle = '+'.join(ligand_names)
    

    # Print input info
    print('')
    print('Input structure:\t', struct)
    print('Input chains:\t\t', user_chains)
    if not user_chains == 'ALL':
        print('Input structchains:\t', user_structchains)
    if autodetect_lig == 1 and ligand_names is None:
        print('Input ligands:\t\tauto-detect')
    elif autodetect_lig == 1 and ligand_names is not None:
        print('Input ligands:\t\t', ligand_names, '+ auto-detect')
    else:
        print('Input ligands:\t\t', ligand_names) #, '\t', ligand_names_bundle)
    if position is not None:
        print('Input position:\t\t', position)
    print('Done\n')


    # Download or get path to query structure & ligand
    try:
        struct_path = download_mmCIF_gz2(struct, pathSTRUCTS)
        print('Loading structure:\t', struct_path, '\n')
    except:
        print('Error downloading structure:\t', struct, '\n')
        sys.exit(1)
        # TODO move into parse_query to fail fast
    
    # Verify ligands (moved)
    '''
    if autodetect_lig == 0 or ligand_names is not None:
        print('Verifying ligands:\t', ligand_names)
        for lig_id in ligand_names:
            try:
                lig_path = download_mmCIF_lig(lig_id, pathLIGS)
                with open(lig_path, 'r') as in_lig:
                    for line in in_lig:
                        if line.startswith('_chem_comp.name'):
                            lig_name = line.split()[1:]
                        if line.startswith('_chem_comp.pdbx_synonyms'):
                            lig_syn = line.split()[1:]
                            print(lig_id, ' '.join(lig_name), ' '.join(lig_syn))
                            break
            except:
                print('Error verifying ligand:\t', lig_id)
    '''


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
    #user_chains_bundle = '+'.join(user_chains)
    print('Input chains verified:\t', user_structchains, user_chains)
    if len(discarded_chains) > 0:
        print('Input chains rejected:\t', discarded_chains)



    ## Look up query in history, if found, return job path and end script
    if autodetect_lig == 0:
        user_input_parameters = struct + '_' + ','.join(user_chains) + '_' + ','.join(ligand_names)
    else:
        user_input_parameters = struct + '_' + ','.join(user_chains) + '_autodetect_lig'
    query_full = user_input_parameters + '_' + settings_str + '-' + job_id
    #print(query_full)
    if look_in_archive == 1:
        print('\nLooking for the same job in history')
        old_same_job = search_query_history(pathQRS, query_full.split('-')[0], 'queries.txt')
        # print('Result of lookup:', old_same_job)
        if old_same_job == 0:
            print('Same job not found, continuing process\n')
        else:
            if os.path.isdir(os.path.dirname(path_results) + '/' + old_same_job):
                print('Same job found in history, printing path of old job results: ', os.path.dirname(path_results) + '/' + old_same_job, '\n')
                print('Printing alignments of old job')
                with open(os.path.dirname(path_results) + '/' + old_same_job + '/' + 'results.csv', 'r') as old_in:
                    for line in old_in:
                        print(line[:-1])
                print('\nDeleting new job folder', job_id)
                if os.path.isdir(path_results):
                    os.rmdir(path_results)
                print('Done\nExiting')
                sys.exit(0)
            print('Old job directory not found, continuing process\n')
    print('')
    
    
    # Get apo candidates from rSIFTS dict, calculate sequence overlap
    # Look for longest UniProt mapping of query chains
    dictApoCandidates = dict()
    uniprot_overlap = dict()
    
    for user_structchain in user_structchains:
        print('\nLooking for longest UniProt mapping for query chain', user_structchain)
        query_uniprot_lengths = list()
        query_uniprot_lengths_dict = dict()
        
        # Iterate through UniProt IDs and their structchains
        for key, values in dict_rSIFTS.items():  # key = UniProt accession, values = structchains + SP_BEG + SP_END
            # Iterate through each structchain
            for i in values:  # i = structchain SP_BEG SP_END
                # Find the match(es) for the query structchain
                if i.split()[0] == user_structchain:
                    
                    uniprot_id = key  # UniProt ID where query belongs
                    structchain_s = i.split()[0] # pass structchain to variable
                    x1s = i.split()[1]      # query SP BEG
                    x2s = i.split()[2]      # query SP END
                    xlength = i.split()[3]  # query SP length

                    query_uniprot_lengths_dict[xlength] = uniprot_id + ' ' + structchain_s + ' ' + x1s + ' ' + x2s
                    query_uniprot_lengths.append(int(xlength))
        
        # Find longest mapping
        top_xlength = max(query_uniprot_lengths)
        top_mapping = query_uniprot_lengths_dict[str(top_xlength)]
        uniprot_id = top_mapping.split()[0]
        user_structchain = top_mapping.split()[1]
        x1 = top_mapping.split()[2]
        x2 = top_mapping.split()[3]
        
        print(query_uniprot_lengths_dict)
        print('->', query_uniprot_lengths_dict[str(top_xlength)], top_xlength)

        # Find candidates overlap for longest stretch
        print('Calculating apo/holo candidate overlap for mapping:', uniprot_id, user_structchain, x1, x2, top_xlength)
        #for values in dict_rSIFTS[uniprot_id]:
        for candidate in dict_rSIFTS[uniprot_id]:
            if candidate.split()[0][:4] != user_structchain[:4]: # discard any structure same with query (stricter)
            
                y1 = candidate.split()[1]    # candidate SP BEG
                y2 = candidate.split()[2]    # candidate SP END
    
                # Calculate overlap and % of overlap on query (x1,x2)
                result =  min(int(x2), int(y2)) - max(int(x1), int(y1))
                percent = result / (int(x2) - int(x1)) * 100
                percent = round(percent, 1)     # round the float
    
                # Build dict with calculated overlap
                #uniprot_overlap.setdefault(i.split()[0], []).append(candidate.split()[0]+' '+str(percent))
                uniprot_overlap.setdefault(candidate.split()[0], []).append(user_structchain + ' ' + str(percent))
                
                # Only consider positive overlap (negative overlap may occur cause of wrong numbering)
                if overlap_threshold != 0 and percent >= overlap_threshold or overlap_threshold == 0 and percent > 0:
                    dict_key = user_structchain+' '+x1+' '+x2
                    dictApoCandidates.setdefault(dict_key, []).append(candidate+' '+str(result)+' '+str(percent))

        print(f'Total chains for {uniprot_id}, {len(dict_rSIFTS[uniprot_id])}')
        print(f'Candidate chains over user-specified overlap threshold [{overlap_threshold}%]:\t{len(dictApoCandidates[dict_key])}') # - {dictApoCandidates[dict_key]}') #[dict_key][0].split()[0]}]')
    
    total_chains = sum([len(dictApoCandidates[x]) for x in dictApoCandidates if isinstance(dictApoCandidates[x], list)])
    print(f'\nTotal candidate chains over user-specified overlap threshold [{overlap_threshold}%]:\t{total_chains}\n')



    ## Apo/holo candidate evaluation

    # Put all structures for downloading into set
    apo_candidate_structs = set()
    for key, values in dictApoCandidates.items():
        for i in values:    # iterate over the values in each key, i = struct/chain combo & SP_BEG SP_END etc
            struct = i.split()[0][:4] # split value strings to get structure only
            apo_candidate_structs.add(struct)
    print('Total structures to be parsed: ', len(apo_candidate_structs), '\n')

    # Download/load the Apo candidate structures to specified directory [TODO this should be replaced by load later]
    for apo_candidate_structure in apo_candidate_structs:
        try:
            download_mmCIF_gz2(apo_candidate_structure, pathSTRUCTS)
        except Exception as ex1:
            #template = "Exception {0} occurred. \n\t\t\t\t\tArguments:{1!r}"
            #message = template.format(type(ex1).__name__, ex1.args) + apo_candidate_structure
            #add_log(message, log_file_dnld)
            print(f'*apo file {apo_candidate_structure} not found, removing from candindates list')

            # Instead of fail, remove structure from queue
            discarded_chains.append(apo_candidate_structure + '\t' + 'PDB structure not found\n')
            apo_candidate_structs.remove(apo_candidate_structure)


    # Parse (mmCIF) structures to get resolution & method. Apply cut-offs
    print('Checking resolution and experimental method of Apo candidate structures')
    for apo_candidate_struct in apo_candidate_structs:
        apo_candidate_structPath = download_mmCIF_gz2(apo_candidate_struct, pathSTRUCTS)
        resolution = '?\t'
        with gzip.open (apo_candidate_structPath, 'rt') as mmCIFin:
            for line in mmCIFin:
                try:
                    if line.split()[0] == '_exptl.method':
                        method = line.split("'")[1]  # capture experimental method #method = ' '.join(line.split()[1:])
                        if method == 'SOLUTION NMR':  # break fast if method is 'NMR'
                            break
                    elif line.split()[0] == '_refine.ls_d_res_high' and float(line.split()[1]):
                        resolution = round(float(line.split()[1]), 3)  # X-ray highest resolution
                        break
                    elif line.split()[0] == '_em_3d_reconstruction.resolution' and float(line.split()[1]):
                        resolution = round(float(line.split()[1]), 3)  # EM resolution
                        break
                except:  # Exception as ex: # getting weird but harmless exceptions
                    #print('Problem parsing structure: ', apo_candidate_struct)
                    pass  # ignore and hide exceptions from stdout
            try:
                if NMR == 1 and method == 'SOLUTION NMR' and xray_only == 0 or xray_only == 1 and method == 'X-RAY DIFFRACTION' and resolution <= res_threshold or xray_only == 0 and resolution <= res_threshold:
                    print(apo_candidate_struct, ' resolution:\t', resolution, '\t', method, '\tPASS')  # Xray/EM
                else:
                    discarded_chains.append(apo_candidate_struct + '\t' + 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']\n')
                    print(apo_candidate_struct, ' resolution:\t', resolution, '\t', method, '\tFAIL')  # Xray/EM
            except:
                discarded_chains.append(apo_candidate_struct + '\t' + 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']\n')
                print('*Exception', apo_candidate_struct, ' resolution:\t', resolution, '\t', method, '\t\tFAIL')  # NMR

    print('Done\n')


    # Compile set of structures to be discarded (non verbose)
    discard_structs = set()
    for i in discarded_chains:
        if len(i.split()[0]) == 4:
            discard_structs.add(i.split()[0])

    # Discard apo entries below threshold(s). Put remainder into new dict, discard UniProt numbering
    print('Structures to discard:\t', len(discard_structs))
    if len(discard_structs) > 0:
        print('Discarding structures\t', discard_structs)
    dictApoCandidates_1 = dict()
    for key, values in dictApoCandidates.items():
        for i in values[:apo_chain_limit]:
            if i.split()[0][:4] in discard_structs:
                print('removing apo', i.split()[0], 'from holo', key.split()[0])
                #pass
            else:
                dictApoCandidates_1.setdefault(key.split()[0], []).append(i.split()[0])
    
    eligible_chains = sum([len(dictApoCandidates_1[x]) for x in dictApoCandidates_1 if isinstance(dictApoCandidates_1[x], list)])
    print(f'\nCandidate chains satisfying user requirements (method/resolution) [{res_threshold} Å ]:\t{eligible_chains}')
    print(dictApoCandidates_1) # helper print, delete later


    # Define ligand search query
    if position is not None: # (4 args) assumes that everything is specified (chains(taken care of), 1 ligand, position) ignore_autodetect_lig
        
        if ligand_names[0] == 'HOH': # mark selection & change lig scan radius
            lig_scan_radius = '3'
            search_term = 'resi ' + position + ' and resn ' + ligand_names_bundle

        elif ligand_names[0] in nonstd_rsds: # find ligands
            if nonstd_rsds_as_lig == 1: # mark selection
                search_term = 'resi ' + position + ' and resn ' + ligand_names_bundle
            elif nonstd_rsds_as_lig == 0: # treat as residue, find ligand

                search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'

        elif ligand_names[0] in d_aminoacids:
            if d_aa_as_lig == 1: # mark selection
                search_term = 'resi ' + position + ' and resn ' + ligand_names_bundle
            elif d_aa_as_lig == 0: # treat as residue, find ligands
                if water_as_ligand == 1:
                    search_term = 'hetatm near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                else:
                    search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'

        elif ligand_names[0] in std_rsds: # find ligands
            if water_as_ligand == 1:
                search_term = 'hetatm near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                print('\n*Search term = ', search_term)
            else:
                search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                print('\n*Search term = ', search_term)
        
        else: # find ligands
            print('\nUnaccounted-for query selection case, using default search\n') # TODO quit?
            print(ligand_names)
            sys.exit(1)
            search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
            print('\n*Search term = ', search_term)

    else: # position == None (3 or less args)
        
        # Configure autodetect lig search expression according to variable
        if autodetect_lig == 1:
            autodetect_lig_expression = ' or (hetatm and not solvent and not polymer)'
        else:
            autodetect_lig_expression = ''
    
        if ligand_names is not None: # ligands specified
            
            if ligand_names[0] in nonstd_rsds:
                print('\n====== Ligands specified + auto-detecting ligands ======\n')
                search_term = 'resn ' + ligand_names_bundle + autodetect_lig_expression
                print('\n*Search term = ', search_term)
        
        elif ligand_names is None and autodetect_lig == 1: 
            
            search_term = 'hetatm and not solvent and not polymer' # water as ligand should be ignored here without a position
        
        #elif ligand_names is None and autodetect_lig == 0: # TODO Here is the chance to force reverse search
            # Here the user has not specified ligands, and they have not turned autodetect_lig ON

                
             
        
        else: # ligands not specified (*/?) 
            print('\n====== No ligands specified: auto-detecting all ligands ======\n')
            search_term = 'hetatm and not solvent and not polymer'
            print('\n*Search term = ', search_term)
            #print('\n No ligands were specified and user has auto-detect OFF, turning it ON to continue')
            autodetect_lig = 1 # Force autodetect ON

        '''
        if autodetect_lig == 1 and ligand_names is not None: # 
            print('\n====== Ligands specified + auto-detecting ligands ======\n')
            search_term = 'resn ' + ligand_names_bundle + ' or (hetatm and not solvent and not polymer)'
            print('\n*Search term = ', search_term)
        elif autodetect_lig == 1 and ligand_names is None:
            print('\n====== No ligands specified: auto-detecting all ligands ======\n')
            search_term = 'hetatm and not solvent and not polymer'
            print('\n*Search term = ', search_term)
        else:
            search_term = 'hetatm and resn ' + ligand_names_bundle
            print('\n*Search term = ', search_term)
        
        #cmd.select(apo_structchain + '_arnd_' + ligand_, 'model ' + apo_struct + '& hetatm & not solvent' + ' near_to ' + lig_scan_radius + ' of holo_' + ligand_)
        '''
    '''
    elif beyond_hetatm == 1: # To be tested
        search_term = 'resn ' + ligand_names_bundle
        print('\n*Search term = ', search_term)
    
    '''

    ##################  Query ligand detection  ##################

    holo_lig_positions = dict()
    apo_holo_dict = dict()
    apo_holo_dict_H = dict()
    for holo_structchain, apo_structchains in dictApoCandidates_1.items():
        print('')
        print(f'=== Processing query chain {holo_structchain} ===')
        holo_struct = holo_structchain[:4]
        holo_chain = holo_structchain[4:]
        holo_struct_path = download_mmCIF_gz2(holo_struct, pathSTRUCTS)

        # Initialize PyMOL but don't reload Holo if present
        if holo_struct in cmd.get_object_list('all'):  # object names
            #if holo_struct in cmd.get_names('all'): # object and selection names
            cmd.delete('not ' + holo_struct)
        else:
            cmd.reinitialize('everything')
            cmd.load(holo_struct_path)

        cmd.select(holo_struct + holo_chain, holo_struct + '& chain ' + holo_chain)


        # Find & name specified ligands in query structure
        ligands_selection = cmd.select('query_ligands', holo_struct +' and '+ search_term + ' and chain ' + holo_chain)  # resn<->name
        if ligands_selection == 0:
            print('No ligands found in author chain, trying PDB chain')
            ligands_selection = cmd.select('query_ligands', holo_struct +' and '+ search_term + ' and segi ' + holo_chain)
            if ligands_selection == 0 and reverse_search == 0:
                print('No ligands found in PDB chain, skipping: ', holo_structchain)
                continue
            elif ligands_selection == 0 and reverse_search == 1:
                print('No ligands found in PDB chain\n====== Reverse mode active, considering input as apo ======\n')
                reverse_mode = True
                cmd.delete('query_ligands')


        if not reverse_mode:  # If query is not APO

            # Identify atom IDs of selected ligand atoms
            ligands_atoms = cmd.identify('query_ligands', mode=0)
            
            # Get positions of specified ligands (better than atom IDs to specify during alignment)
            myspace = {'positions': []}  # temporary dict with fixed key name
            for atom in ligands_atoms:  # this is a bulk of atoms for all ligands, many atoms can belong to a single ligand
                cmd.iterate('id ' + str(atom), 'positions.append(resi +" "+ chain +" "+ resn)', space = myspace)

            # Transfer temporary list with positions to dict
            for key, values in myspace.items():
                for i in values:
                    holo_lig_positions.setdefault(holo_structchain, []).append(i)

            # Remove duplicate values from holo_lig_positions
            for key,value in holo_lig_positions.items():
                holo_lig_positions[key] = list(holo_lig_positions.fromkeys(value))   # preserves the order of values

            print('Ligand information')
            print('Atom IDs: ', ligands_atoms)
            print('Total atoms: ', len(ligands_atoms))
            print('Atom positions/chains/names: ', holo_lig_positions.get(holo_structchain))  #, '/', len(holo_lig_positions.get(holo_structchain)))

            # Name holo ligands as PyMOL selections. Put real (detected) ligand names into set
            holo_lig_names = set()
            for ligand in holo_lig_positions[holo_structchain]:
                resi = ligand.split()[0]
                chain = ligand.split()[1]
                resn = ligand.split()[2]
                ligand_ = ligand.replace(' ', '_')
                holo_lig_names.add(resn)
                cmd.select('holo_' + ligand_, 'model ' + holo_struct + '& resi ' + resi + '& chain ' + chain + '& resn ' + resn) # s1
            if autodetect_lig == 1:
                ligand_names = holo_lig_names.copy()



        # Start candidate chain loop
        # Align chain to query and mark atom selections around superimposed holo ligand binding sites
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
            print('')
            try:
                aln_rms = cmd.align(apo_struct + '& chain ' + apo_chain, holo_struct + '& chain ' + holo_chain, cutoff=10.0, cycles=1)
                aln_tm = tmalign2(cmd, apo_struct + '& chain ' + apo_chain, holo_struct + '& chain ' + holo_chain, quiet=1, transform=0)
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

            # Start ligand look-up loop for apo chain
            if not reverse_mode:
                found_ligands = set()
                found_ligands_xtra = set()
                for ligand in holo_lig_positions[holo_structchain]:

                    # Around selection [looks for ligands in every (valid) chain alignment, not just the standard locus of holo ligand(s)]
                    ligand_ = ligand.replace(' ', '_') # remove spaces for selection name
                    #s2 = cmd.select(apo_structchain + '_arnd_' + ligand_, 'model ' + apo_struct + '& chain ' + apo_chain + ' near_to ' + lig_scan_radius + ' of holo_' + ligand_)
                    cmd.select(apo_structchain + '_arnd_' + ligand_, 'model ' + apo_struct + '& hetatm & not solvent' + ' near_to ' + lig_scan_radius + ' of holo_' + ligand_)
                    #cmd.select(apo_structchain + '_arnd_' + ligand_, apo_struct + '& hetatm & not (polymer or solvent)' + ' near_to ' + lig_scan_radius + ' of holo_' + ligand_)

                    # Put selected atoms in a list, check their name identifiers to see if holo ligand name is present
                    myspace_a = {'a_positions': []}
                    for a_atom in cmd.identify(apo_structchain + '_arnd_' + ligand_):
                        cmd.iterate('id ' + str(a_atom), 'a_positions.append(resi +" "+ chain +" "+ resn)', space = myspace_a)

                    # Transfer dict[key] values (just resn) into set for easier handling
                    apo_lig_names = set()
                    for a_position in myspace_a['a_positions']:
                        a_atom_lig_name = a_position.split()[2]
                        apo_lig_names.add(a_atom_lig_name)

                    # Assess apo ligands
                    for i in apo_lig_names:
                        if i in holo_lig_names or i in ligand_names:    # check in both lists (detected holo ligs + query ligs)   # TODO holo_lig_names may be undefined
                            found_ligands.add(i)    #print('Holo ligand found in Apo: ', apo_structchain, i)
                        elif i not in nolig_resn:
                            found_ligands_xtra.add(i)

            # Start reverse mode, if query = APO (no ligands). Find identical structures with/wo ligands
            else:
                found_ligands_r = set()
                found_ligands_xtra = set()
                # Find ligands in holo candidate
                cmd.select('holo_ligands_' + apo_structchain, apo_struct + '& chain ' + apo_chain + '& hetatm & not (polymer or solvent)')
                myspace_r = {'r_positions': []}
                for r_atom in cmd.identify('holo_ligands_' + apo_structchain, mode=0):
                    cmd.iterate('id ' + str(r_atom), 'r_positions.append(resi +" "+ chain +" "+ resn)', space = myspace_r)
                # Transfer dict[key] values (just resn) into set for easier handling
                for r_position in myspace_r['r_positions']:
                    r_atom_lig_name = r_position.split()[2]
                    if r_atom_lig_name not in nolig_resn:  # exclude non ligands
                        found_ligands_r.add(r_atom_lig_name)


            # Print verdict for chain & save it as ".cif.gz"
            if not reverse_mode:

                print(f'*query ligand(s)/position: {ligand_names}/[{position}]\tdetected ligands: {holo_lig_names}\t detected apo ligands: {apo_lig_names}\tfound query ligands: {found_ligands}\tfound non-query ligands: {found_ligands_xtra}')
                
                # Apo
                if lig_free_sites == 1 and len(found_ligands_xtra) == 0 and len(found_ligands) == 0 or lig_free_sites == 0 and len(found_ligands) == 0:
                    ligands_str = join_ligands(found_ligands.union(found_ligands_xtra))
                    apo_holo_dict.setdefault(holo_structchain, []).append(apo_structchain + ' ' + uniprot_overlap[apo_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    if len(found_ligands_xtra) > 0:
                        print('APO*')
                    else:
                        print('APO')
                    if save_separate == 1:
                        if not os.path.isfile(path_results + '/holo_' + holo_struct + '.cif.gz'):
                            cmd.save(path_results + '/holo_' + holo_struct + '.cif.gz', holo_struct) # save query structure
                        cmd.save(path_results + '/a_' + apo_structchain + '_aln_to_' + holo_structchain + '.cif.gz', apo_structchain) # save apo chain
                # Holo
                else:
                    ligands_str = join_ligands(found_ligands.union(found_ligands_xtra))
                    apo_holo_dict_H.setdefault(holo_structchain, []).append(apo_structchain + ' ' + uniprot_overlap[apo_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    if len(found_ligands) > 0:
                        print('HOLO')
                    else:
                        print('HOLO*')
                    if save_separate == 1 and save_oppst == 1:
                        if not os.path.isfile(path_results + '/holo_' + holo_struct + '.cif.gz'):
                            cmd.save(path_results + '/holo_' + holo_struct + '.cif.gz', holo_struct) # save query structure
                        cmd.save(path_results + '/h_' + apo_structchain + '_aln_to_' + holo_structchain + '.cif.gz', apo_structchain) # save holo chain

            else:  # reverse mode
                
                print('Found ligands: ', found_ligands_r)  # TODO found_ligands_r may be undefined
                
                # Holo
                if len(found_ligands_r) > 0:
                    ligands_str = join_ligands(found_ligands_r)
                    apo_holo_dict_H.setdefault(holo_structchain, []).append(apo_structchain + ' ' + uniprot_overlap[apo_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    print('HOLO')
                    if save_separate == 1:
                        if not os.path.isfile(path_results + '/holo_' + holo_struct + '.cif.gz'):
                            cmd.save(path_results + '/holo_' + holo_struct + '.cif.gz', holo_struct) # save query structure
                        cmd.save(path_results + '/h_' + apo_structchain + '_aln_to_' + holo_structchain + '.cif.gz', apo_structchain) # save apo chain
                # Apo
                else:
                    ligands_str = join_ligands(found_ligands_r.union(found_ligands_xtra))
                    apo_holo_dict.setdefault(holo_structchain, []).append(apo_structchain + ' ' + uniprot_overlap[apo_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    if len(found_ligands_xtra) > 0:
                        print('APO*')
                    else:
                        print('APO')
                    if save_separate == 1 and save_oppst == 1:
                        if not os.path.isfile(path_results + '/holo_' + holo_struct + '.cif.gz'):
                            cmd.save(path_results + '/holo_' + holo_struct + '.cif.gz', holo_struct)  # save query structure
                        cmd.save(path_results + '/a_' + apo_structchain + '_aln_to_' + holo_structchain + '.cif.gz', apo_structchain)  # save holo chain


        # Clean objects/selections in PyMOL session
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
        cmd.disable('all')  # toggle off the display of all objects
        cmd.enable(holo_struct)
        cmd.deselect()
        cmd.reset()
        try:
            cmd.center('query_ligands')
        except:
            cmd.center(holo_structchain)

        # Save results as session (.pse.gz) or multisave (.cif)
        filename_body = path_results + '/' + 'aln_' + holo_struct     #filename_body = pathRSLTS + '/' + 'aln_' + holo_structchain + '_to_' + '_'.join(cmd.get_object_list('all and not ' + holo_struct))
        filename_pse = filename_body + '.pse.gz'
        filename_multi = filename_body + '_multi.cif'
        if len(apo_holo_dict) > 0:    #len(dictApoCandidates_1) > 0:    #if len(cmd.get_object_list('all')) > 1:
            if save_session == 1:
                cmd.save(filename_pse)
            if multisave == 1:
                cmd.multisave(filename_multi, append=1)

    # end for

    print('')


    ## Save results in text output
    
    #if reverse_mode:    query_chain = 'apo_chain'
    #else:   query_chain = holo_chain

    # Apo results
    num_apo_chains = sum([len(apo_holo_dict[x]) for x in apo_holo_dict if isinstance(apo_holo_dict[x], list)])  # number of found APO chains
    if len(apo_holo_dict) > 0:  #if save_separate == 1 or multisave == 1 or save_session == 1:
        '''
        # Write dictionary to file
        filename_aln = pathRSLTS + '/apo_aln_' + '_'.join(list(apo_holo_dict.keys()))
        if reverse_mode:    header = "#HEADER: {apo_chain: [apo_chain %UniProt_overlap RMSD TM_score]\n"
        else:        header = "#HEADER: {holo_chain: [apo_chain %UniProt_overlap RMSD TM_score]\n"
        with open (filename_aln + '.txt', 'wt') as out1:
            out1.write(header)
            out1.write(str(apo_holo_dict))
        '''
        # Write CSV file
        filename_csv = path_results + '/results_apo.csv'
        if reverse_mode:
            header = "#apo_chain,apo_chain,%UniProt_overlap,RMSD,TM_score,ligands\n"
        else:
            header = "#holo_chain,apo_chain,%UniProt_overlap,RMSD,TM_score,ligands\n"
        #header = "#holo_chain,apo_chain,%UniProt_overlap,RMSD,TM_score\n"
        with open(filename_csv, 'w') as csv_out:
            csv_out.write(header)
            for key, values in apo_holo_dict.items():
                for value in values:
                    csv_out.write("%s,%s\n" % (key, ','.join(value.split())))

        # Print apo dict
        print('Apo results: ')
        for key in apo_holo_dict:
            print(key, apo_holo_dict.get(key))
    else:
        print('No apo forms found')


    # Holo results
    num_holo_chains = sum([len(apo_holo_dict_H[x]) for x in apo_holo_dict_H if isinstance(apo_holo_dict_H[x], list)])  # number of found HOLO chains
    if len(apo_holo_dict_H) > 0:
        '''
        # Write dictionary to file
        filename_aln = pathRSLTS + '/holo_aln_' + '_'.join(list(apo_holo_dict_H.keys()))
        if reverse_mode:    header = "#HEADER: {apo_chain: [holo_chain %UniProt_overlap RMSD TM_score ligands]\n"
        else:   header = "#HEADER: {holo_chain: [holo_chain %UniProt_overlap RMSD TM_score ligands]\n"
        with open (filename_aln + '.txt', 'wt') as out1:
            out1.write(header)
            out1.write(str(apo_holo_dict_H))
        '''
        # Write CSV file
        filename_csv = path_results + '/results_holo.csv'
        if reverse_mode:
            header = "#apo_chain,holo_chain,%UniProt_overlap,RMSD,TM_score,ligands\n"
        else:
            header = "#holo_chain,holo_chain,%UniProt_overlap,RMSD,TM_score,ligands\n"
        with open(filename_csv, 'w') as csv_out:
            csv_out.write(header)
            for key, values in apo_holo_dict_H.items():
                for value in values:
                    csv_out.write("%s,%s\n" % (key, ','.join(value.split())))

        # Print holo dict
        print('\nHolo results: ')
        for key in apo_holo_dict_H:
            print(key, apo_holo_dict_H.get(key))
    else:
        print('\nNo holo forms found')

    if len(apo_holo_dict) == 0 and len(apo_holo_dict_H) == 0:
        print('\nConsider reversing the search or revising the input query')
        # Note: we don't want to delete empty job folder but keep it for potential further processing by the webserver
        # Delete empty job folder
        # print('\nDeleting empty job folder')    
        # try:
        #     os.rmdir(pathRSLTS)
        # except OSError as error:
        #     print('Job folder not empty, job ID: ', job_id, error)

    # Append the name of the query and the job_id in the queries.txt
    if job_id:
        print('\nSaving query:', query_full)
        with open(pathQRS + '/' + 'queries.txt', 'a') as out_q:
            out_q.write(query_full + '\n')

    # Print argparse arguments
    #print('\n', args)
    #print(vars(args))

    pm.stop()

    print('\nResults saved to directory:\t', path_results)
    print('\nDone processing query: ', query)

    return QueryResult(
        result_dir=path_results,
        num_apo_chains=num_apo_chains,
        num_holo_chains=num_holo_chains)

# end process_query()


##########################################################################################################

def try_process_query(query, workdir, args, data: PrecompiledData = None) -> QueryResult:
    """ Returns QueryResult with error if calculation fails """
    try:
        return process_query(query, workdir, args, data)
    except Exception as ex:
        print(ex)
        return QueryResult(result_dir=None, error=ex)


def process_queries(query_lines: list, workdir, args, data: PrecompiledData = None) -> list:

    queries = [parse_query(query_str) for query_str in query_lines]  # parse all here to check for possible invalid queries

    if data is None:
        data = load_precompiled_data(workdir)

    def process_q(query):
        return try_process_query(query, workdir, args, data)

    with PoolExecutor(max_workers=args.threads) as pool:
        res = pool.map(process_q, query_lines)

    # TODO better way to print summary results and report errors
    print('\n--------------------')
    print('Results:\n')
    print('\n'.join([str(r) for r in res]))
    print('--------------------\n')

    return res


def parse_args(argv):
    parser = argparse.ArgumentParser()

    # Main user query
    parser.add_argument('--query', type=str,   default='1a73 a zn', help='main input query')
    #parser.add_argument('--query', type=str,   default='1a73 a', help='main input query')
    #parser.add_argument('--query', type=str,   default='1a73 a ser 97', help='main input query')
    
    #parser.add_argument('--query', type=str,   default='1a73 b hoh 509', help='main input query')
    #parser.add_argument('--query', type=str,   default='6h3c b,g zn', help='main input query')
    #parser.add_argument('--query', type=str,   default='2v0v', help='main input query')
    #parser.add_argument('--query', type=str,   default='2hka all c3s', help='main input query')
    #parser.add_argument('--query', type=str,   default='2v57 a,c prl', help='main input query')
    
    


    # Basic
    parser.add_argument('--res_threshold',     type=float, default=3.8,  help='resolution cut-off for apo chains (angstrom), condition is <=')
    parser.add_argument('--NMR',               type=int,   default=1,    help='0/1: discard/include NMR structures')
    parser.add_argument('--xray_only',         type=int,   default=0,    help='0/1: only consider X-ray structures')
    parser.add_argument('--lig_free_sites',    type=int,   default=1,    help='0/1: resulting apo sites will be free of any other known ligands in addition to specified ligands')
    parser.add_argument('--autodetect_lig',    type=int,   default=0,    help='0/1: if the user does not know the ligand, auto detection will consider non-protein heteroatoms as ligands')
    parser.add_argument('--reverse_search',    type=int,   default=0,    help='0/1: look for holo structures from apo')

    # Advanced
    parser.add_argument('--overlap_threshold', type=float, default=0,    help='% of overlap between apo and holo chain (w UniProt numbering), condition is ">=", "0" will not allow (erroneously) negative overlap')
    parser.add_argument('--lig_scan_radius',   type=float, default=5,    help='angstrom radius to look around holo ligand(s) superposition (needs to be converted to str)')
    parser.add_argument('--min_tmscore',       type=float, default=0.5,  help='minimum acceptable TM score for apo-holo alignments (condition is "<" than)')

    # Experimental
    parser.add_argument('--water_as_ligand',   type=int,   default=0,    help='0/1: consider HOH atoms as ligands (can be used in combination with lig_free_sites)(strict)')
    parser.add_argument('--nonstd_rsds_as_lig',type=int,   default=0,    help='0/1: ignore/consider non-standard residues as ligands')
    parser.add_argument('--d_aa_as_lig',       type=int,   default=0,    help='0/1: ignore/consider D-amino acids as ligands')
    parser.add_argument('--beyond_hetatm',     type=int,   default=0,    help='0/1: when enabled, does not limit holo ligand detection to HETATM records for specified ligand/residue')  # [might need to apply this to apo search too #TODO]
    parser.add_argument('--look_in_archive',   type=int,   default=0,    help='0/1: search if the same query has been processed in the past (can give very fast results)')

    # Internal
    parser.add_argument('--apo_chain_limit',   type=int,   default=999,  help='limit number of apo chains to consider when aligning (for fast test runs)')
    parser.add_argument('--work_dir',          type=str,   default=None, help='global root working directory for pre-computed and intermediary data')
    parser.add_argument('--out_dir',           type=str,   default=None, help='explicitly specified output directory')
    parser.add_argument('--threads',           type=int,   default=4,    help='number of concurrent threads for processing multiple queries')
    # Saving
    parser.add_argument('--save_oppst',        type=int,   default=1,    help='0/1: also save chains same with query (holo chains when looking for apo, and apo chains when looking for holo)')
    parser.add_argument('--save_separate',     type=int,   default=1,    help='0/1: save each chain object in a separate file (default save)')
    parser.add_argument('--save_session',      type=int,   default=0,    help='0/1: save each result as a PyMOL ".pse" session (zipped, includes annotations -less recommended)')
    parser.add_argument('--multisave',         type=int,   default=0,    help='0/1: save each result in a .pdb file (unzipped, no annotations -least recommended)')
    '''
    # print help if there are no arguments
    if len(argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    '''

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    workdir = get_workdir(args)
    query = args.query

    # TODO read multi line queries from file
    query_lines = [q.strip() for q in query.splitlines()]  # TODO ignore blank lines
    if len(query_lines) > 1:
        process_queries(query_lines, workdir, args)
    else:
        process_query(query.strip(), workdir, args)

    print('All done')
    return 0




if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
