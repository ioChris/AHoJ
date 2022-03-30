# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:24:57 2021

@author: ChrisX
"""
# Apo-Holo Juxtaposition - AHoJ
from common import get_workdir, load_dict_binary, tmalign2, write_file

import __main__
__main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
# import pymol.cmd as cmd
#import psico.fitting
#import psico.fullinit
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

#import rich.traceback

#rich.traceback.install(show_locals=True, extra_lines=4, max_frames=1)

VERSION = '0.2.0'


_global_lock = threading.Lock()                      # multi-threading
# global_lock = multiprocessing.Manager().Lock()     # multi-processing (must be moved to main)

'''
Given an experimental protein structure (PDB code), with optionally specified chain(s), ligand(s) and position, find its equivalent apo and holo forms.
The program will look for both apo and holo forms of the query structure. Structures are processed chain by chain.

The user can specify the following input arguments depending on the mode of search
i) When looking for apo from holo:
-Min arguments: PDB code
-Max arguments: PDB code, chain(s), ligand(s), position
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
    print('Input format: <pdb_id> <chain> <ligand> <position> or <pdb_id> <chains> <ligands> or <pdb_id> <chains> or <pdb_id> <ligands> or <pdb_id>')
    print('Input examples: "3fav A,B ZN" or "3fav * ZN" or "3fav ALL ZN" or "3fav"')
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


def is_not_blank(s):
    return bool(s and not s.isspace())

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
                    print('Verifying ligand:\t', ligand_names, '> ', lig_id, ' '.join(lig_name), ' '.join(lig_syn))
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
    #nonstd_rsds = "SEP TPO PSU MSE MSO 1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG PYG PYL SEC PHA ".split() # don't use as no-lig rsds
    d_rsds_hoh = "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA HOH".split()
    nolig_resn = list()
    nolig_resn.extend(std_rsds + d_rsds_hoh)


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

    if chains == '*':
        chains = 'ALL'
    elif not all(chain.isalnum() for chain in chains.split(',')):  # check that chains are alphanumeric characters
        raise ValueError(f"Invalid query '{query}': only alphanumeric characters allowed as chains")
    if ligands == '*':
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


def write_results_apo_csv(apo_holo_dict, path_results, broad_search_mode):
    '''
    # Write dictionary to file
    filename_aln = pathRSLTS + '/apo_aln_' + '_'.join(list(apo_holo_dict.keys()))
    if broad_search_mode:    header = "#HEADER: {apo_chain: [apo_chain %UniProt_overlap RMSD TM_score]\n"
    else:        header = "#HEADER: {query_chain: [apo_chain %UniProt_overlap RMSD TM_score]\n"
    with open (filename_aln + '.txt', 'wt') as out1:
        out1.write(header)
        out1.write(str(apo_holo_dict))
    '''
    # Write CSV file
    filename_csv = path_results + '/results_apo.csv'
    if broad_search_mode:
        header = "#query_apo_chain, apo_chain, %UniProt_overlap, RMSD, TM_score, ligands\n"
    else:
        header = "#query_holo_chain, apo_chain, %UniProt_overlap, RMSD, TM_score, ligands\n"
    #header = "#query_chain,apo_chain,%UniProt_overlap,RMSD,TM_score\n"
    with open(filename_csv, 'w') as csv_out:
        csv_out.write(header)
        for key, values in apo_holo_dict.items():
            for value in values:
                csv_out.write("%s,%s\n" % (key, ','.join(value.split())))


def write_results_holo_csv(apo_holo_dict_H, path_results, broad_search_mode):
    '''
    # Write dictionary to file
    filename_aln = pathRSLTS + '/holo_aln_' + '_'.join(list(apo_holo_dict_H.keys()))
    if broad_search_mode:    header = "#HEADER: {apo_chain: [query_chain %UniProt_overlap RMSD TM_score ligands]\n"
    else:   header = "#HEADER: {query_chain: [query_chain %UniProt_overlap RMSD TM_score ligands]\n"
    with open (filename_aln + '.txt', 'wt') as out1:
        out1.write(header)
        out1.write(str(apo_holo_dict_H))
    '''
    # Write CSV file
    filename_csv = path_results + '/results_holo.csv'
    if broad_search_mode:
        header = "#query_apo_chain, query_chain, %UniProt_overlap, RMSD, TM_score, ligands\n"    # I'm confused here, shouldn't 2. column be a result "holo_chain"?
    else:
        header = "#query_holo_chain, holo_chain, %UniProt_overlap, RMSD, TM_score, ligands\n"
    with open(filename_csv, 'w') as csv_out:
        csv_out.write(header)
        for key, values in apo_holo_dict_H.items():
            for value in values:
                csv_out.write("%s,%s\n" % (key, ','.join(value.split())))


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
    #autodetect_lig = args.autodetect_lig
    #reverse_search = args.reverse_search # should be renamed to "start with apo" or "broad search"
    water_as_ligand = args.water_as_ligand
    
    # Advanced
    overlap_threshold = args.overlap_threshold
    lig_scan_radius = args.lig_scan_radius
    min_tmscore = args.min_tmscore
    nonstd_rsds_as_lig = args.nonstd_rsds_as_lig
    d_aa_as_lig = args.d_aa_as_lig

    # Experimental
    #beyond_hetatm = args.beyond_hetatm
    look_in_archive = args.look_in_archive

    # Internal
    apo_chain_limit = args.apo_chain_limit

    # Saving
    save_oppst = args.save_oppst
    save_separate = args.save_separate
    save_session = args.save_session
    multisave = args.multisave

    # Adjust input, resolve conflicts
    autodetect_lig = 0 # default OFF
    #if reverse_search == 1:
     #   autodetect_lig = 1
    lig_scan_radius = str(lig_scan_radius)  # needs to be str
    broad_search_mode = False

    # Pass settings to a string
    #settings_str = 'res' + str(res_threshold) + '_NMR' + str(NMR) + '_xrayonly' + str(xray_only) + '_ligfree' + str(lig_free_sites) + '_autodtctlig' + str(autodetect_lig) + '_reverse' + str(reverse_search) + '_h2olig' + str(water_as_ligand) + '_overlap' + str(overlap_threshold) + '_ligrad' + str(lig_scan_radius) + '_tmscore' + str(min_tmscore) + '_nonstdaas' + str(nonstd_rsds_as_lig) + '_daas' + str(d_aa_as_lig)
    settings_str = 'res' + str(res_threshold) + '_NMR' + str(NMR) + '_xrayonly' + str(xray_only) + '_ligfree' + str(lig_free_sites) + '_autodtctlig' + str(autodetect_lig) + '_h2olig' + str(water_as_ligand) + '_overlap' + str(overlap_threshold) + '_ligrad' + str(lig_scan_radius) + '_tmscore' + str(min_tmscore) + '_nonstdaas' + str(nonstd_rsds_as_lig) + '_daas' + str(d_aa_as_lig)


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

    # Verify input structure here  # TODO move verifications into parse_query to fail fast
    try:
        struct_path = download_mmCIF_gz2(struct, pathSTRUCTS)
        print('Verifying structure:', struct, '\t> ', struct_path.split('/')[-1]) #'/'.join(struct_path.split('/')[-3:]))
    except:
        raise ValueError(f"Invalid PDB ID in query '{query}': use a valid 4-letter PDB code")

    # Verify ligand/residue names here
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
        print('Input ligands:\t\t', 'auto-detect')
    elif autodetect_lig == 1 and ligand_names is not None: # this scenario should not occur anymore and thus be removed
        print('Input ligands:\t\t', ligand_names, '+ auto-detect')
    else:
        print('Input ligands:\t\t', ligand_names) #, '\t', ligand_names_bundle)
    if position is not None:
        print('Input position:\t\t', position)
    print('Done\n')

    '''
    # Download or get path to query structure & ligand
    try:
        struct_path = download_mmCIF_gz2(struct, pathSTRUCTS)
        print('Loading structure:\t', struct_path, '\n')
    except:
        print('Error downloading structure:\t', struct, '\n')
        sys.exit(1)
    '''



    ## Find Apo candidates (rSIFTS)

    # Find & VERIFY input chains by UniProt ID (if they don't exist in uniprot, we cannot process them)
    # TODO: allow non-UniProt chains, because ligands can be assigned to non-polymer chains
    print(f'\nFinding & verifying query chains {user_chains} by UniProt ID')
    discarded_chains = list()   # Discarded chains (format: structchain + '\t' + discard_msg)
    usr_structchains_unverified = list()
    
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
                #user_chains.remove(user_structchain[4:])
                user_structchains.remove(user_structchain)
                discarded_chains.append(user_structchain + '\t' + 'No assigned UniProt ID\n')
                usr_structchains_unverified.append(user_structchain)
                
    #user_chains_bundle = '+'.join(user_chains)
    print('Input chains verified:\t\t', user_structchains)#, user_chains)
    print('Input chains unverified:\t', usr_structchains_unverified)
    if len(discarded_chains) > 0:
        print('Input chains rejected:\t', discarded_chains)



    # Form the full query expression - to use for indexing the query and searching past jobs
    if autodetect_lig == 0:
        user_query_parameters = struct + '_' + ','.join(user_chains) + '_' + ','.join(ligand_names)
    elif autodetect_lig == 1 and ligand_names is not None:
        user_query_parameters = struct + '_' + ','.join(user_chains) + '_' + ','.join(ligand_names) + '_autodtctligloc'
    else:
        user_query_parameters = struct + '_' + ','.join(user_chains) + '_autodtctlig'
    query_full = user_query_parameters + '_' + settings_str + '-' + job_id
    #print(query_full)


    ## Look up query in history, if found, return job path and end script
    if look_in_archive == 1:
        print('\nLooking for the same job in job history')
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
    uniprot_overlap_all = dict() # mostly for reference/monitoring
    uniprot_overlap = dict() # overlap of successful hits

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

        # Find candidates overlap for longest mapping
        print('Calculating apo/holo candidate overlap for mapping:', uniprot_id, user_structchain, x1, x2, top_xlength)
        #for values in dict_rSIFTS[uniprot_id]:
        own_chains = list() # the structchains that belong to the same UniProt ID and same structure with query
        for candidate in dict_rSIFTS[uniprot_id]:
            if candidate.split()[0][:4] != user_structchain[:4]: # filter out any structure same with query (stricter)

                y1 = candidate.split()[1]    # candidate SP BEG
                y2 = candidate.split()[2]    # candidate SP END

                # Calculate overlap and % of overlap on query (x1,x2)
                result =  min(int(x2), int(y2)) - max(int(x1), int(y1))
                percent = result / (int(x2) - int(x1)) * 100
                percent = round(percent, 1)     # round the float

                # Build dict with calculated overlap
                #uniprot_overlap.setdefault(i.split()[0], []).append(candidate.split()[0]+' '+str(percent))
                uniprot_overlap_all.setdefault(candidate.split()[0], []).append(user_structchain + ' ' + str(percent)) # this keeps all calculated overlaps (from candidate to all query chains)
                #uniprot_overlap[candidate.split()[0]] = user_structchain + ' ' + str(percent)

                # Only consider positive overlap (negative overlap may occur cause of wrong numbering)
                if overlap_threshold != 0 and percent >= overlap_threshold or overlap_threshold == 0 and percent > 0:
                    dict_key = user_structchain + ' ' + x1 + ' ' + x2
                    dictApoCandidates.setdefault(dict_key, []).append(candidate+' '+str(result)+' '+str(percent))
                    uniprot_overlap.setdefault(candidate.split()[0], []).append(user_structchain + ' ' + str(percent)) # this keeps only successful overlaps
            else:
                own_chains.append(candidate)

        print(f'Total chains for {uniprot_id}: {len(dict_rSIFTS[uniprot_id])}')
        print(f'Total chains for {uniprot_id} excluding query structure: {len(dict_rSIFTS[uniprot_id]) - len(own_chains)}')
        #print(dict_rSIFTS[uniprot_id], own_chains)

        if len(dict_rSIFTS[uniprot_id]) > len(own_chains):
            print(f'Candidate chains over user-specified overlap threshold [{overlap_threshold}%]:\t{len(dictApoCandidates[dict_key])}') # - {dictApoCandidates[dict_key]}') #[dict_key][0].split()[0]}]')
        else:
            print('No other UniProt chains found')
            #sys.exit(2) # TODO don't exit script, continue for other UniProt chains within same struct or multi-query

    # Handle unverified chains [first find all valid UniProt chains of the structure, and then their candidates]
    #if len(usr_structchains_unverified) > 0:
    #for unverified_structchain in usr_structchains_unverified:
        #TODO

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
        except: # Exception as ex1:
            #template = "Exception {0} occurred. \n\t\t\t\t\tArguments:{1!r}"
            #message = template.format(type(ex1).__name__, ex1.args) + apo_candidate_structure
            #add_log(message, log_file_dnld)
            print(f'*apo file {apo_candidate_structure} not found, removing from candindates list')

            # Don't fail, instead remove candidate structure from queue
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


    # Compile set of structures to be discarded (discard comments)
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
    print(f'\nCandidate chains satisfying user requirements (method/resolution) [{res_threshold} Å]:\t{eligible_chains}')
    print(dictApoCandidates_1) # helper print, delete later


    # Configure autodetect lig search expression according to variable
    if autodetect_lig == 1:
        autodetect_lig_expression = ' or (hetatm and not solvent and not polymer)'
    else:
        autodetect_lig_expression = ''


    # Define ligand search query
    if position is not None: # (4 args) assumes that everything is specified (chains(taken care of), 1 ligand, position) ignore_autodetect_lig
        
        if ligand_names[0] == 'HOH': # mark selection & change lig scan radius (unless user has set it to lower than 3)
            if float(lig_scan_radius) > 3:
                lig_scan_radius = '3'
            search_term = 'resi ' + position + ' and resn ' + ligand_names_bundle
            print('\n*Search term = ', search_term)

        elif ligand_names[0] in nonstd_rsds: # find ligands
            if nonstd_rsds_as_lig == 1: # mark selection
                search_term = 'resi ' + position + ' and resn ' + ligand_names_bundle
            elif nonstd_rsds_as_lig == 0: # treat as residue, find ligand
                search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                print('\n*Search term = ', search_term)

        elif ligand_names[0] in d_aminoacids:
            if d_aa_as_lig == 1: # mark selection
                search_term = 'resi ' + position + ' and resn ' + ligand_names_bundle
                print('\n*Search term = ', search_term)
            elif d_aa_as_lig == 0: # treat as residue, find ligands
                if water_as_ligand == 1:
                    search_term = 'hetatm near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                    print('\n*Search term = ', search_term)
                else:
                    search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                    print('\n*Search term = ', search_term)

        elif ligand_names[0] in std_rsds: # find ligands
            if water_as_ligand == 1:
                search_term = 'hetatm near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                print('\n*Search term = ', search_term)
            else:
                search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'
                print('\n*Search term = ', search_term)

        else: # find ligands # normal ligand? Error: ligand has different chain
            #print('\nUnaccounted-for query selection case, using default search\n') # TODO examine several examples
            print('\n', ligand_names, query)
            search_term = '(resn ' + ligand_names_bundle + ' and resi ' + position + autodetect_lig_expression + ')'
            print('\n*Search term = ', search_term)

    else: # position == None (3 or less args)

        if ligand_names is not None: # ligands specified

            if ligand_names[0] in nonstd_rsds:
                search_term = '(resn ' + ligand_names_bundle + autodetect_lig_expression + ')'
                print('\n*Search term = ', search_term)

            else: # assume ligand is a real ligand
                search_term = '(hetatm and resn ' + ligand_names_bundle + autodetect_lig_expression + ')'
                print('\n*Search term = ', search_term)

        elif ligand_names is None and autodetect_lig == 1: 
            search_term = 'hetatm and not solvent and not polymer' # water as ligand should be ignored here without a position
            print('\n*Search term = ', search_term)

        #elif ligand_names is None and autodetect_lig == 0: # TODO Here is the chance to force reverse search
            # Here the user has not specified ligands, and they have not turned autodetect_lig ON
            #reverse_search = 1 # TEST doesn't work because autodetect has been forced ON already

        else: # ligands not specified (*/?) 
            print('\n====== No ligands specified: auto-detecting all ligands ======\n')
            search_term = 'hetatm and not solvent and not polymer'
            print('\n*Search term = ', search_term)
            #print('\n No ligands were specified and user has auto-detect OFF, turning it ON to continue')
            autodetect_lig = 1 # Force autodetect ON



    ##################  Query ligand detection  ##################

    query_lig_positions = dict()
    apo_holo_dict = dict()
    apo_holo_dict_H = dict()

    progress_total_candidates = sum([len(lst) for lst in dictApoCandidates_1.values()])
    progress_processed_candidates = 0


    def track_progress(write_results: bool = False):
        if args.track_progress:
            write_file(path_results + '/.progress', f"{progress_processed_candidates}/{progress_total_candidates}")
            if write_results:
                write_results_apo_csv(apo_holo_dict, path_results, broad_search_mode)
                write_results_holo_csv(apo_holo_dict_H, path_results, broad_search_mode)


    for query_structchain, candidates_structchains in dictApoCandidates_1.items():
        print('')
        print(f'=== Processing query chain {query_structchain} ===')

        query_struct = query_structchain[:4]
        query_chain = query_structchain[4:]
        query_struct_path = download_mmCIF_gz2(query_struct, pathSTRUCTS)

        # Initialize PyMOL but don't reload Holo if present
        if query_struct in cmd.get_object_list('all'):  # object names
            #if query_struct in cmd.get_names('all'): # object and selection names
            cmd.delete('not ' + query_struct)
        else:
            cmd.reinitialize('everything')
            cmd.load(query_struct_path)

        # Name selection for query structchain
        cmd.select(query_structchain, query_struct + '& chain ' + query_chain)

        # Find & name selection for user-specified (or autodetect) ligands in query structure
        ligands_selection = cmd.select('query_ligands', query_struct +' and '+ search_term + ' and chain ' + query_chain)
        if ligands_selection == 0:
            print('No ligands found in author chain, trying PDB chain')
            ligands_selection = cmd.select('query_ligands', query_struct +' and '+ search_term + ' and segi ' + query_chain)
            #if ligands_selection == 0 and reverse_search == 0:
            if ligands_selection == 0 and autodetect_lig == 0:
                print('No ligands found in PDB chain, skipping: ', query_structchain)
                continue
            #elif ligands_selection == 0 and reverse_search == 1:
            elif ligands_selection == 0 and autodetect_lig == 1:
                print('No ligands found in PDB chain\n====== Reverse mode active, considering input as apo ======\n')
                broad_search_mode = True
                cmd.delete('query_ligands')


        if not broad_search_mode:  # When query is not fully APO

            # Identify atom IDs of selected ligand atoms
            ligands_atoms = cmd.identify('query_ligands', mode=0)  # bulk atoms of all query ligands

            # Get positions of specified ligands (iterate through atoms)
            myspace = {'positions': []}
            cmd.iterate('query_ligands', 'positions.append(resi +" "+ chain +" "+ resn)', space = myspace) # iterates through atoms of selection

            # Transfer temporary list with positions to dict
            for key, values in myspace.items():
                for i in values:
                    query_lig_positions.setdefault(query_structchain, []).append(i)

            # Remove duplicate values from query_lig_positions
            for key, value in query_lig_positions.items():
                query_lig_positions[key] = list(query_lig_positions.fromkeys(value))   # preserves the order of values

            print('Query ligand information')
            print('Total atoms: ', len(ligands_atoms))
            print('Atom IDs: ', ligands_atoms)
            print('Position/chain/name: ', query_lig_positions.get(query_structchain))  #, '/', len(query_lig_positions.get(query_structchain)))

            # Name query ligands as seperate selections per "residue"/position. Put real (detected) ligand names into set
            query_lig_names = set()
            for ligand in query_lig_positions[query_structchain]:
                resi = ligand.split()[0]
                chain = ligand.split()[1]
                resn = ligand.split()[2]
                ligand_ = ligand.replace(' ', '_')
                query_lig_names.add(resn)
                cmd.select('holo_' + ligand_, query_struct + '& resi ' + resi + '& chain ' + chain + '& resn ' + resn) # s1
            if autodetect_lig == 1 and ligand_names is None:
                ligand_names = query_lig_names.copy() # ligand_names = user-specified ligands, when no ligands specified, they might be undefined


        # Start candidate chain loop
        # Align candidate chain to query chain and mark atom selections around superimposed query ligand binding sites
        for candidate_structchain in candidates_structchains:
            progress_processed_candidates += 1
            track_progress(write_results=True)

            candidate_struct = candidate_structchain[:4]
            candidate_chain = candidate_structchain[4:]
            candidate_struct_path = download_mmCIF_gz2(candidate_struct, pathSTRUCTS)

            if candidate_struct in cmd.get_object_list('all'):
                pass
            else:
                cmd.load(candidate_struct_path)
            cmd.select(candidate_structchain, candidate_struct + '& chain ' + candidate_chain)

            # Align candidate to query
            print('')
            try:
                aln_rms = cmd.align(candidate_struct + '& chain ' + candidate_chain, query_struct + '& chain ' + query_chain, cutoff=2.0, cycles=1)
                aln_tm = tmalign2(cmd, candidate_struct + '& chain ' + candidate_chain, query_struct + '& chain ' + query_chain, quiet=1, transform=0)
                print('Alignment RMSD/TM score:', candidate_structchain, query_structchain, round(aln_rms[0], 3), aln_tm)
            except:
                discarded_chains.append(candidate_structchain + '\t' + 'Alignment error\n')
                print('Alignment RMSD/TM score: ERROR')
                print('*poor alignment (error), discarding chain ', candidate_structchain)
                continue

            # Discard poor alignments
            if aln_tm < min_tmscore:
                discarded_chains.append(candidate_structchain + '\t' + 'Poor alignment [RMSD/TM]: ' + str(round(aln_rms[0], 3)) +'/'+ str(aln_tm) + '\n')
                print('*poor alignment (below threshold), discarding chain ', candidate_structchain)
                continue


            # Look for ligands in candidate chain
            if not broad_search_mode:
                found_ligands = set()
                found_ligands_xtra = set()
                #print('All query ligands:', query_lig_positions[query_structchain])
                #print(f'Assessing detected ligands in {query_lig_names} and {ligand_names}')
                for ligand in query_lig_positions[query_structchain]:
                    #print('scanning ligand:', ligand)
                    resi = ligand.split()[0]
                    chain = ligand.split()[1]
                    resn = ligand.split()[2]
                    ligand_ = ligand.replace(' ', '_') # remove spaces for selection name

                    # Around selection: look for candidate ligands in the superimposed sites of the aligned query ligands # TODO possible to extract binding site residues here
                    cndt_sele_expression = '(' + candidate_struct + ' and hetatm)' # only limit search to candidate structure (not chain) #cndt_sele_expression = candidate_struct + ' and chain ' + candidate_chain + ' and hetatm'
                    qr_lig_sele_expression = '(' + query_struct + ' and chain ' + chain + ' and resi ' + resi + ' and resn ' + resn + ')'

                    #full_sele =  cndt_sele_expression + ' near_to ' + lig_scan_radius + ' of ' + qr_lig_sele_expression
                    #print('full sele:', full_sele)
                    cmd.select(candidate_structchain + '_arnd_' + ligand_, cndt_sele_expression + ' near_to ' + lig_scan_radius + ' of ' + qr_lig_sele_expression)
                    # If cndt ligand belongs to different PDB chain than candidate structchain, it will not be saved, we need to merge the two selections here
                    cmd.select(candidate_structchain, candidate_structchain + '_arnd_' + ligand_, merge=1)

                    # Put selected atoms in a list, check their name identifiers to see if query ligand name is present
                    myspace_cndt = {'cndt_positions': []}
                    cmd.iterate(candidate_structchain + '_arnd_' + ligand_, 'cndt_positions.append(resi +" "+ chain +" "+ resn)', space = myspace_cndt)

                    # Remove duplicate values from myspace_cndt
                    for key, value in myspace_cndt.items():
                        myspace_cndt[key] = list(myspace_cndt.fromkeys(value))   # preserves the order of values - better than set()

                    print(f'-candidate ligands in query ligand binding site [{ligand}]: {myspace_cndt["cndt_positions"]}')

                    # Transfer dict[key] values (just resn) into set for easier handling
                    cndt_lig_names = set()
                    for cndt_position in myspace_cndt['cndt_positions']: # use set to remove redundant positions
                        cndt_atom_lig_name = cndt_position.split()[2]
                        cndt_lig_names.add(cndt_atom_lig_name)
                    #print('Found candidate ligands:', cndt_lig_names)


                    # Assess apo ligands
                    for i in cndt_lig_names:
                        if i in query_lig_names or i in ligand_names:  # check in both lists (detected holo ligs + query ligs)  # TODO query_lig_names may be undefined
                            found_ligands.add(i)
                        elif i not in nolig_resn:
                            found_ligands_xtra.add(i)

            # Start reverse mode, if query = fully APO (no ligands). Find identical structures with/wo ligands
            else:
                found_ligands_r = set()
                found_ligands_xtra = set()
                # Find ligands in holo candidate
                cmd.select('holo_ligands_' + candidate_structchain, candidate_struct + '& chain ' + candidate_chain + '& hetatm & not (polymer or solvent)')
                myspace_r = {'r_positions': []}
                cmd.iterate('holo_ligands_' + candidate_structchain, 'r_positions.append(resi +" "+ chain +" "+ resn)', space = myspace_r)
                # Remove duplicate values from myspace_r
                for key, value in myspace_r.items():
                    myspace_r[key] = list(myspace_r.fromkeys(value))   # preserves the order of values
                print(f'-candidate ligands found: {myspace_r["r_positions"]}') # use set to remove redundant positions
                # Transfer dict[key] values (just resn) into set for easier handling
                for r_position in myspace_r['r_positions']: # use set to remove redundant positions
                    r_atom_lig_name = r_position.split()[2]
                    if r_atom_lig_name not in nolig_resn:  # exclude non ligands
                        found_ligands_r.add(r_atom_lig_name)


            # Print verdict for chain & save it as ".cif.gz"
            if not broad_search_mode:

                print(f'*query ligand(s)/position: {ligand_names}/[{position}]\t detected ligands: {query_lig_names}\t detected candidate ligands: {cndt_lig_names}\t found query ligands: {found_ligands}\t found non-query ligands: {found_ligands_xtra}')

                # Apo
                if lig_free_sites == 1 and len(found_ligands_xtra) == 0 and len(found_ligands) == 0 or lig_free_sites == 0 and len(found_ligands) == 0:
                    ligands_str = join_ligands(found_ligands.union(found_ligands_xtra))
                    apo_holo_dict.setdefault(query_structchain, []).append(candidate_structchain + ' ' + uniprot_overlap[candidate_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    if len(found_ligands_xtra) > 0:
                        print('APO*')
                    else:
                        print('APO')
                    if save_separate == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results + '/query_' + query_struct + '.cif.gz', query_struct) # save query structure
                        cmd.save(path_results + '/apo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain) # save apo chain

                # Holo
                else:
                    ligands_str = join_ligands(found_ligands.union(found_ligands_xtra))
                    apo_holo_dict_H.setdefault(query_structchain, []).append(candidate_structchain + ' ' + uniprot_overlap[candidate_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    if len(found_ligands) > 0:
                        print('HOLO')
                    else:
                        print('HOLO*')
                    if save_separate == 1 and save_oppst == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results + '/query_' + query_struct + '.cif.gz', query_struct) # save query structure
                        #save_sele = (candidate_structchain + '_arnd_' + ligand_) for ligand in query_lig_positions[query_structchain])
                        cmd.save(path_results + '/holo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain) # save holo chain


            else:  # reverse mode

                print('Found ligands: ', found_ligands_r)  # TODO found_ligands_r may be undefined

                # Holo
                if len(found_ligands_r) > 0:
                    ligands_str = join_ligands(found_ligands_r)
                    apo_holo_dict_H.setdefault(query_structchain, []).append(candidate_structchain + ' ' + uniprot_overlap[candidate_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    print('HOLO')
                    if save_separate == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results + '/query_' + query_struct + '.cif.gz', query_struct) # save query structure
                        cmd.save(path_results + '/holo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain) # save apo chain
                # Apo
                else:
                    ligands_str = join_ligands(found_ligands_r.union(found_ligands_xtra))
                    apo_holo_dict.setdefault(query_structchain, []).append(candidate_structchain + ' ' + uniprot_overlap[candidate_structchain][0].split()[1] + ' ' + str(round(aln_rms[0], 3)) + ' ' + str(round(aln_tm, 3)) + ' ' + ligands_str)
                    if len(found_ligands_xtra) > 0:
                        print('APO*')
                    else:
                        print('APO')
                    if save_separate == 1 and save_oppst == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results + '/query_' + query_struct + '.cif.gz', query_struct)  # save query structure
                        cmd.save(path_results + '/apo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain)  # save holo chain



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
            if i[:4] == query_struct or i[:4] in good_selections:
                pass
            else:
                cmd.delete(i)
        all_selections = cmd.get_names('all')
        for i in all_selections:    # Delete 0 atom selections
            if cmd.count_atoms(i) == 0:
                cmd.delete(i)
        cmd.disable('all')  # toggle off the display of all objects
        cmd.enable(query_struct)
        cmd.deselect()
        cmd.reset()
        try:
            cmd.center('query_ligands')
        except:
            cmd.center(query_structchain)


        # Save results as session (.pse.gz) or multisave (.cif)
        filename_body = path_results + '/' + 'aln_' + query_struct     #filename_body = pathRSLTS + '/' + 'aln_' + query_structchain + '_to_' + '_'.join(cmd.get_object_list('all and not ' + query_struct))
        filename_pse = filename_body + '.pse.gz'
        filename_multi = filename_body + '_multi.cif'
        if len(apo_holo_dict) > 0:    #len(dictApoCandidates_1) > 0:    #if len(cmd.get_object_list('all')) > 1:
            if save_session == 1:
                cmd.save(filename_pse)
            if multisave == 1:
                cmd.multisave(filename_multi, append=1)

    track_progress()
    # end for

    print('')

    # Print chains that were discarded
    if len(discarded_chains) > 0:
        print('Discarded candidate chains: ', len(discarded_chains))
        #print(f"{' '.join(map(str, discarded_chains))}\n")
    
    # Print calculated UniProt overlap
    #print('Uniprot overlap dictionary\n', uniprot_overlap)

    ## Save results in text output

    #if broad_search_mode:    query_chain = 'apo_chain'
    #else:   query_chain = query_chain


    # Apo results
    write_results_apo_csv(apo_holo_dict, path_results, broad_search_mode)
    num_apo_chains = sum([len(apo_holo_dict[x]) for x in apo_holo_dict if isinstance(apo_holo_dict[x], list)])  # number of found APO chains
    if len(apo_holo_dict) > 0:  #if save_separate == 1 or multisave == 1 or save_session == 1:
        # Print apo dict
        print('\nApo chains: ', num_apo_chains)
        for key in apo_holo_dict:
            print(key, apo_holo_dict.get(key))
    else:
        print('No apo forms found')


    # Holo results
    write_results_holo_csv(apo_holo_dict_H, path_results, broad_search_mode)
    num_holo_chains = sum([len(apo_holo_dict_H[x]) for x in apo_holo_dict_H if isinstance(apo_holo_dict_H[x], list)])  # number of found HOLO chains
    if len(apo_holo_dict_H) > 0:
        # Print holo dict
        print('\nHolo chains: ', num_holo_chains)
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
        print(f'\nSaving query search: {[user_query_parameters]} for query {[query]}') #query_full)  # Don't show full query string
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

    # TODO validate all queries before processing

    if data is None:
        data = load_precompiled_data(workdir)

    def process_q(query):
        return try_process_query(query, workdir, args, data)

    with PoolExecutor(max_workers=args.threads) as pool:
        results = list(pool.map(process_q, query_lines))

    # TODO better way to print summary results and report errors
    print('\n--------------------')
    print('Results:\n')
    print('\n'.join([str(r) for r in results]))
    print('--------------------\n')

    return results


def parse_args(argv):
    parser = argparse.ArgumentParser()

    # Main user query
    # Ligand
    #parser.add_argument('--query', type=str,   default='1a73 a zn', help='main input query') # OK apo 0, holo 16
    #parser.add_argument('--query', type=str,   default='1a73 a,b zn', help='main input query') # OK apo 0, holo 32
    #parser.add_argument('--query', type=str,   default='1a73 * zn', help='main input query') # OK apo 0, holo 32
    #parser.add_argument('--query', type=str,   default='1a73 a zn', help='main input query') # reverse_search=1, OK apo 0, holo 16
    #parser.add_argument('--query', type=str,   default='1a73 * zn', help='main input query') # reverse_search=1, OK apo 0, holo 32
    #parser.add_argument('--query', type=str,   default='1a73 e mg 205', help='main input query') # fail, ligand assigned non-polymer chain
    #parser.add_argument('--query', type=str,   default='1a73 a', help='main input query') # apo 0, holo 16
    #parser.add_argument('--query', type=str,   default='1a73 * *', help='main input query') # OK apo 0, holo 32
    #parser.add_argument('--query', type=str,   default='5j72 a na 703', help='main input query') # apo 0, holo 0 (no UniProt chains)
    #parser.add_argument('--query', type=str,   default='1a73 b mg 206', help='main input query') # OK, apo 4, holo 12
    #parser.add_argument('--query', type=str,   default='1a73 b mg 206', help='main input query') # water_as_ligand=1 OK, apo 4, holo 12
    #parser.add_argument('--query', type=str,   default='1a73 b mg 206', help='main input query') # autodetect_lig=1 OK, apo 0, holo 16
    #parser.add_argument('--query', type=str,   default='7s4z a *', help='main input query') # apo 103, holo 104 *many irrelevant ligands
    #parser.add_argument('--query', type=str,   default='3fav all zn', help='main input query') 
    #parser.add_argument('--query', type=str,   default='3fav all', help='main input query') 
    #parser.add_argument('--query', type=str,   default='1y57 a mpz', help='main input query') 
    
    #parser.add_argument('--query', type=str,   default='6h3c b,g zn', help='main input query') # OK apo 0, holo 4
    #parser.add_argument('--query', type=str,   default='2v0v', help='main input query') # OK apo 0, holo 0
    #parser.add_argument('--query', type=str,   default='2v0v', help='main input query') # reverse_search=1, apo 8, holo 24
    #parser.add_argument('--query', type=str,   default='2hka all c3s', help='main input query') # OK apo 2, holo 0
    #parser.add_argument('--query', type=str,   default='2v57 a,c prl', help='main input query') # OK apo 4, holo 0
    #parser.add_argument('--query', type=str,   default='3CQV all hem', help='main input query') # OK apo 6, holo 5
    #parser.add_argument('--query', type=str,   default='2npq a bog', help='main input query') # long, apo 149, holo 114, p38 MAP kinase cryptic sites
    #parser.add_argument('--query', type=str,   default='1ksw a NBS', help='main input query') # apo 4, holo 28 Human c-Src Tyrosine Kinase (Thr338Gly Mutant) in Complex with N6-benzyl ADP
    parser.add_argument('--query', type=str,   default='1ai5', help='main input query') # negative uniprot overlap (fixed)
    
    # Residue
    #parser.add_argument('--query', type=str,   default='1a73 a ser', help='main input query') # expected parsing fail
    #parser.add_argument('--query', type=str,   default='1a73 a ser 97', help='main input query') # OK apo 4, holo 12
    
    # SARS CoV 2
    #parser.add_argument('--query', type=str,   default='7aeh a R8H', help='main input query') # apo 116, holo 431, job 314 H PC, CoV2 Mpro caspase-1 inhibitor SDZ 224015
    #parser.add_argument('--query', type=str,   default='2amq a his 164', help='main input query') # apo 41, holo 66, job 318 H, Crystal Structure Of SARS_CoV Mpro in Complex with an Inhibitor N3
    #parser.add_argument('--query', type=str,   default='6lu7 a R8H', help='main input query')

    # Water
    #parser.add_argument('--query', type=str,   default='1a73 b hoh 509', help='main input query') # OK apo 9, holo 7
    #parser.add_argument('--query', type=str,   default='1a73 * hoh 509', help='main input query') # OK apo 9, holo 7
    #parser.add_argument('--query', type=str,   default='1a73 * hoh', help='main input query') # expected parsing fail
    #parser.add_argument('--query', type=str,   default='3i34 x hoh 311', help='main input query') # apo 113, holo 94 *many irrelevant ligands show up
    
    #parser.add_argument('--query', type=str,   default='1pkz a tyr 9', help='main input query') # apo 7, holo 93, water as lig, marian, allosteric effect of hoh
    
    # Non standard residues
    #parser.add_argument('--query', type=str,   default='6sut a tpo', help='main input query') # OK apo 0, holo 3
    #parser.add_argument('--query', type=str,   default='6sut a tpo 285', help='main input query') # OK apo 0, holo 3
    #parser.add_argument('--query', type=str,   default='6sut a tpo,*', help='main input query') # OK apo 0, holo 3

    #parser.add_argument('--query', type=str,   default='1a73 a zn 201', help='main input query') # OK apo 0, holo 16
    

    # Basic
    parser.add_argument('--res_threshold',     type=float, default=3.8,   help='Lowest allowed resolution for result structures (applies to highest resolution value for scattering methods, expressed in angstroms), condition is <=')
    parser.add_argument('--NMR',               type=int,   default=1,     help='0/1: Discard/include NMR structures')
    parser.add_argument('--xray_only',         type=int,   default=0,     help='0/1: Only consider X-ray structures')
    parser.add_argument('--lig_free_sites',    type=int,   default=1,     help='0/1: Ligand-free binding sites. When on, resulting apo sites will be free of any other known ligands in addition to specified ligands')
    #parser.add_argument('--autodetect_lig',    type=int,   default=0,     help='0/1: This will find and consider any non-protein and non-solvent heteroatoms as ligands and mark their binding sites, in addition to any specified ligands (useful when the user does not know the ligand)')
    #parser.add_argument('--reverse_search',    type=int,   default=0,     help='0/1: Start the search with an apo structure that does not bind any ligands')
    parser.add_argument('--water_as_ligand',   type=int,   default=0,     help='0/1: Consider HOH atoms as ligands when examining the superimposed candidate binding sites (can be used in combination with lig_free_sites - strict condition)')

    # Advanced
    parser.add_argument('--overlap_threshold', type=float, default=0,     help='Minimum % of sequence overlap between query and result chains (using the SIFTS residue-level mapping with UniProt), condition is ">="')
    parser.add_argument('--lig_scan_radius',   type=float, default=4.5,   help='Angstrom radius to look around the query ligand(s) superposition (needs to be converted to str)')
    parser.add_argument('--min_tmscore',       type=float, default=0.5,   help='Minimum acceptable TM score for apo-holo alignments (condition is "<" than)')
    parser.add_argument('--nonstd_rsds_as_lig',type=int,   default=0,     help='0/1: Ignore/consider non-standard residues as ligands')
    parser.add_argument('--d_aa_as_lig',       type=int,   default=0,     help='0/1: Ignore/consider D-amino acids as ligands')

    # Experimental
    #parser.add_argument('--beyond_hetatm',     type=int,   default=0,     help='0/1: when enabled, does not limit holo ligand detection to HETATM records for specified ligand/residue')  # [might need to apply this to apo search too #TODO remove?]
    parser.add_argument('--look_in_archive',   type=int,   default=0,     help='0/1: Search if the same query has been processed in the past (can give very fast results)')

    # Internal
    parser.add_argument('--apo_chain_limit',   type=int,   default=999,   help='limit number of apo chains to consider when aligning (for fast test runs)')
    parser.add_argument('--work_dir',          type=str,   default=None,  help='global root working directory for pre-computed and intermediary data')
    parser.add_argument('--out_dir',           type=str,   default=None,  help='explicitly specified output directory')
    parser.add_argument('--threads',           type=int,   default=4,     help='number of concurrent threads for processing multiple queries')
    parser.add_argument('--track_progress',    type=bool,  default=False, help='track the progress of long queries in .progress file, update result csv files continually (not just at the end)')
    # Saving
    parser.add_argument('--save_oppst',        type=int,   default=1,     help='0/1: also save chains same with query (holo chains when looking for apo, and apo chains when looking for holo)')
    parser.add_argument('--save_separate',     type=int,   default=1,     help='0/1: save each chain object in a separate file (default save)')
    parser.add_argument('--save_session',      type=int,   default=0,     help='0/1: save each result as a PyMOL ".pse" session (zipped, includes annotations -less recommended)')
    parser.add_argument('--multisave',         type=int,   default=0,     help='0/1: save each result in a .pdb file (unzipped, no annotations -least recommended)')
    '''
    # print help if there are no arguments
    if len(argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    '''

    return parser.parse_args(argv)


def main(argv):
    print(f'APO-HOLO JUXTAPOSITION (v{VERSION})')
    args = parse_args(argv)

    workdir = get_workdir(args)
    query = args.query

    # TODO read multi line queries from file
    query_lines = [q.strip() for q in query.splitlines() if is_not_blank(q)]


    if len(query_lines) > 1:
        results = process_queries(query_lines, workdir, args)

        # Report errors
        error_results = [r for r in results if r.error is not None]
        if error_results:
            print(f'ERRORS: {len(error_results)} of {len(results)} queries finished with errors:')
            for qr in error_results:
                print(f' >  {qr.error}')
            print('Some queries finished with ERROR.')
            sys.exit(1)

    else:
        process_query(query.strip(), workdir, args)

    print('All done')
    return 0




if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
