# -*- coding: utf-8 -*-
"""
@author: ChrisX
"""
# Apo-Holo Juxtaposition - AHoJ
import shutil

'''
Given a protein structure (PDB code), with optionally specified chain(s), ligand(s) and their position,
find its equivalent apo and holo forms, for each chain separately.
When no ligand(s) are specified, the program will detect and consider all ligands in the query.

The user can specify the following input arguments
i) When looking for apo from holo:
-Min arguments: PDB code
-Max arguments: PDB code, chain(s), ligand(s), position

*When position is specified, only one ligand can be defined.
'''

VERSION = '0.4.8'

import copy
import pathlib

from common import get_workdir, load_dict_binary, tmalign2, write_file
from residue_mapping import map_pdb_resnum_to_uniprot, group_mapped_res_by_chain, examine_cndt_mapped_bs_res, remove_negative_duplicate_cndt_bs_res_pos, evaluate_candidate_bs_rsds, print_dict_readable, good_candidates_from_residue_mapping, bad_candidates_from_residue_mapping, get_scores_from_residue_mapping

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

from concurrent.futures import ThreadPoolExecutor; import threading           # multi-threading
from concurrent.futures import ProcessPoolExecutor; import multiprocessing    # multi-processing (doesn't work atm)

#import rich.traceback

#rich.traceback.install(show_locals=True, extra_lines=4, max_frames=1)


_global_lock = threading.Lock()                      # multi-threading
# global_lock = multiprocessing.Manager().Lock()     # multi-processing (must be moved to main)

# TODO add sync with PDB

##########################################################################################################
# Define functions
##########################################################################################################

def download_mmCIF_gz2(pdb_id, pdb_dir):   # Version 2 of download mmCIF gz (without exception handling)
    # in pdb_dir mimic directory structure of FTP/rsynced whole PDB
    # e.g.: 4ZZW is in {pdb_dir}/zz/4zzw.cif.gz

    urlA = 'https://files.rcsb.org/download/'
    ext = '.cif.gz'
    url = urlA + pdb_id.upper() + ext
    pdb_id = pdb_id.lower()
    middle_bit = pdb_id[1:3]
    subdir = f'{pdb_dir}/{middle_bit}'
    file_path = f'{subdir}/{pdb_id}{ext}'

    if not os.path.isfile(file_path):
        pathlib.Path(subdir).mkdir(exist_ok=True)
        print(f'Downloading: {pdb_id + ext}')
        wget.download(url, subdir)
        return file_path
    else:
        return file_path
    # TODO(rdk): solve problem where structure is downloaded at the same time by multiple processes and saved as '3hku.cif (1).gz'


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


def download_sifts_xml_gz(pdb_id, sifts_dir):
    urlA = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/'
    ext = '.xml.gz'
    middle_bit = pdb_id[1:3]
    pdb_id = pdb_id.lower()
    url = urlA + f'{middle_bit}/' + pdb_id + ext
    subdir = f'{sifts_dir}/{middle_bit}'
    file_path = f'{subdir}/{pdb_id}{ext}'
    if not os.path.isfile(file_path):
        pathlib.Path(subdir).mkdir(exist_ok=True)
        print(f'Downloading: {pdb_id + ext}')
        wget.download(url, subdir)
        return file_path
    else:
        return file_path


def merge_fragmented_unp_overlaps(fragmented_overlap_dict):
    merged_overlap_dict = dict()
    for key, values in fragmented_overlap_dict.items(): # key = candidate structchain, values = <query_structchain %UNP_overlap>
        if len(values) > 1:
            #total_overlap = 0
            temp_dict = dict()
            for value in values:
                chain = value.split()[0]
                overlap = value.split()[1]
                if chain in temp_dict.keys(): # structchain has already a unp fragment
                    new_overlap = float(temp_dict[chain]) + float(overlap)
                    temp_dict[chain] = str(round(new_overlap, 1))
                else:
                    temp_dict[chain] = overlap
                #total_overlap += float(overlap)
            #merged_overlap_dict.setdefault(key, []).append(chain + ' ' + str(round(total_overlap, 1)))
            temp_list = list()
            for x, y in temp_dict.items():
                temp_list.append(x + ' ' + y)
            merged_overlap_dict[key] = temp_list
        else:
            merged_overlap_dict[key] = values
    return merged_overlap_dict


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


def wrong_input_error():
    print('\n=== ERROR: Wrong input format ===')
    print('-use a whitespace character to separate input arguments\n-chains are case-sensitive')
    print('\nInput format structure:\n<pdb_id> <chain> <ligand/residue> <position> or\n<pdb_id> <chains> <ligands> or\n<pdb_id> <chains> or\n<pdb_id> <ligands> or\n<pdb_id>')
    print('\nInput examples:\n"3fav A ZN 101"\n"3fav A,B ZN"\n"3fav ! ZN"\n"3fav * ZN"\n"3fav ALL ZN"\n"3fav"\n')
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
    water_as_ligand_auto: bool
    nonstd_rsds_as_lig_auto: bool
    d_aa_as_lig_auto: bool


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
class CandidateChainResult:
    """ Result of candidate chain evaluation """
    # TODO remodel, include/remove attributes
    query_struct: str
    query_chain: str
    candidate_struct: str
    candidate_chain: str
    #query_lig_positions: dict
    passed: bool = False
    discard_reason: str = None
    tm_score: float = None
    tm_score_i: float = None
    rmsd: float = None
    #bndg_rsd_ratio: str
    #bndg_rsd_percent: str

    apo_holo_dict_instance: dict = None
    apo_holo_dict_H_instance: dict = None
    cndt_lig_positions_instance: dict = None


@dataclass
class PrecompiledData:
    """
    Pre-compiled data needed for computation
    UniProt PDB mapping
    """
    dict_SIFTS: dict   # regular SIFTS dictionary
    dict_rSIFTS: dict  # reverse SIFTS (SPnum) dictionary


def verify_ligands(ligand_names, pathLIGS):

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
        #except Exception as ex:
            #print('Error verifying ligand:\t', lig_id, ex)


def parse_query(query: str, autodetect_lig: bool = False, water_as_ligand_auto: bool = False, nonstd_rsds_as_lig_auto: bool = False, d_aa_as_lig_auto: bool = False) -> Query: #

    # Parse single line input (line by line mode, 1 structure per line)
    # If only first argument specified (structure but no chains), consider all chains
    print(f"Parsing query '{query}'")
    query = query.strip()
    parts = query.split()

    struct = parts[0].lower()
    chains = 'ALL'
    ligands = None
    position = None

    # Define non-ligands (3-letter names of amino acids and h2o)
    std_rsds = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    nonstd_rsds = "SEP TPO PSU MSE MSO 1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG PYG PYL SEC PHA".split()
    d_rsds = "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA".split()
    #nolig_resn = list()    #nolig_resn.extend(std_rsds + nonstd_rsds + d_rsds_hoh)


    if len(struct) != 4:
        raise ValueError(f"Invalid query '{query}': '{struct}' is not a valid PDB structure code")

    if len(parts) == 1:
        autodetect_lig = 1  # overrides cmd line param
    elif len(parts) == 2:
        chains = parts[1]
        autodetect_lig = 1
    elif len(parts) == 3:
        chains = parts[1]
        ligands = parts[2].upper()       # adjust case, ligands = upper

    # When position is specified, there has to be a single ligand/residue specified
    elif len(parts) == 4 and len(parts[2]) < 4 and len(parts[2].split(',')) == 1:# and int(parts[3]):
        try:
            chains = parts[1]
            ligands = parts[2].upper()
            position = str(int(parts[3]))   # test if int
            #if ligands in std_rsds:
                #autodetect_lig = 1 # not needed
                #print('\nLigand is standard residue')
        except Exception:
            raise ValueError(f"Invalid query '{query}': wrong number of parts")
    else:
        raise ValueError(f"Invalid query '{query}': wrong number of parts")

    if chains == '*' or chains == 'all':
        chains = 'ALL'
    elif chains == '!' and ligands is not None: # this will only process chain(s) that include query ligand(s)
        pass
    elif not all(chain.isalnum() for chain in chains.split(',')):  # check that chains are alphanumeric characters
        raise ValueError(f"Invalid query '{query}': only alphanumeric characters allowed as chains or '!/*'")


    # If ligand is HOH or std residue or non-std residue,  make sure that:
    # i) there is just one specified (handled earlier)
    # ii) there is a fourth argument (position)

    if ligands is not None:
        for ligand in ligands.split(','):
            ligand = ligand.upper()
            if ligand in std_rsds and position is None:
                raise ValueError(f"Invalid query '{query}': specify index position of residue")
            elif ligand == 'HOH' and position is not None:
                water_as_ligand_auto = 1
            elif ligand == 'HOH' and position is None:
                raise ValueError(f"Invalid query '{query}': specify index position of HOH")
            elif ligand in nonstd_rsds:
                nonstd_rsds_as_lig_auto = 1
            elif ligand in d_rsds:
                d_aa_as_lig_auto = 1
            elif ligand == '*':  # overide any other ligands in input and replace with autodetect all ligands
                ligands = None
                autodetect_lig = 1
                break

    return Query(struct=struct, chains=chains, ligands=ligands, position=position, autodetect_lig=autodetect_lig, water_as_ligand_auto=water_as_ligand_auto, nonstd_rsds_as_lig_auto=nonstd_rsds_as_lig_auto, d_aa_as_lig_auto=d_aa_as_lig_auto)


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
    #fileSIFTSdict = pathSIFTS + '/pdb_chain_uniprot_dict.bin'
    #fileRSIFTS = pathSIFTS + '/pdb_chain_uniprot_REVERSE_SPnum.bin'

    # Observed residues only
    fileSIFTSdict = pathSIFTS + '/uniprot_segments_observed_dict.bin'
    fileRSIFTS = pathSIFTS + '/uniprot_segments_observed_REVERSE_SPnum.bin'

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


def write_ligands_csv(query_lig_positions, cndt_lig_positions, path_results):  # Write dict(s) to csv
    # we don't want to edit original dicts while computation is still running
    # doing deepcopy because they are dicts of lists
    cndt_lig_positions = copy.deepcopy(cndt_lig_positions)
    query_lig_positions = copy.deepcopy(query_lig_positions)

    filename_csv = path_results + '/ligands.csv'
    header = "#chain, ligand_positions\n"

    with open(filename_csv, 'w') as csv_out:
        csv_out.write(header)

        # Write query ligands
        for key, values in query_lig_positions.items():
            for idx, item in enumerate(values):  # replace " " with "_"
                values[idx] = item.replace(" ", "_")
            csv_out.write("%s,%s\n" % (key, '-'.join(values)))

        # Write (holo or apo) candidate ligands
        for key, values in cndt_lig_positions.items():
            for idx, item in enumerate(values):  # replace " " with "_"
                values[idx] = item.replace(" ", "_")
            csv_out.write("%s,%s\n" % (key, '-'.join(values)))


def write_results_apo_csv(apo_holo_dict, path_results):

    # Write CSV file
    filename_csv = path_results + '/results_apo.csv'
    header = "#query_chain, apo_chain, Resolution, R-free, %UniProt_overlap, Mapped_bndg_rsds, %Mapped_bndg_rsds, RMSD, TM_score, iTM_score, ligands\n"
    with open(filename_csv, 'w') as csv_out:
        csv_out.write(header)
        for key, values in apo_holo_dict.items():
            for value in values:
                csv_out.write("%s,%s\n" % (key, ','.join(value.split())))


def write_results_holo_csv(apo_holo_dict_H, path_results):

    # Write CSV file
    filename_csv = path_results + '/results_holo.csv'
    header = "#query_chain, apo_chain, Resolution, R-free, %UniProt_overlap, Mapped_bndg_rsds, %Mapped_bndg_rsds, RMSD, TM_score, iTM_score, ligands\n"
    with open(filename_csv, 'w') as csv_out:
        csv_out.write(header)
        for key, values in apo_holo_dict_H.items():
            for value in values:
                csv_out.write("%s,%s\n" % (key, ','.join(value.split())))


def write_discarded_chains(discarded_chains_dict, path_results):

    # Write txt file
    filename_txt = path_results + '/discarded_chains.txt'
    #header = 'Discarded_chain, Discard_reason'
    with open(filename_txt, 'w') as txt_out:
        #txt_out.write(header)
        for key, values in discarded_chains_dict.items():
            txt_out.write("%s\t%s\n" % (key, ' | '.join(values)))


def write_pymol_script(path_results):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    pymol_file = "load_results_into_PyMOL.pml"
    # copyfile() should be fast and does not copy permissions and metadata which is what we want
    shutil.copyfile(f"{script_dir}/{pymol_file}", f"{path_results}/{pymol_file}")


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
    include_nmr = args.include_nmr
    xray_only = args.xray_only
    lig_free_sites = args.lig_free_sites
    #autodetect_lig = args.autodetect_lig
    #reverse_search = args.reverse_search # should be renamed to "start with apo" or "broad search"
    water_as_ligand_usr = args.water_as_ligand

    # Advanced
    overlap_threshold = args.overlap_threshold
    bndgrsds_threshold = args.bndgrsds_threshold
    lig_scan_radius = args.lig_scan_radius
    min_tmscore = args.min_tmscore
    nonstd_rsds_as_lig_usr = args.nonstd_rsds_as_lig
    d_aa_as_lig_usr = args.d_aa_as_lig

    # Experimental
    #beyond_hetatm = args.beyond_hetatm
    look_in_archive = args.look_in_archive

    # Internal
    apo_chain_limit = args.apo_chain_limit
    intrfc_lig_radius = args.intrfc_lig_radius
    hoh_scan_radius = args.hoh_scan_radius # TODO replace with dynamic function for scan radius

    # Saving
    save_oppst = args.save_oppst
    save_separate = args.save_separate
    #save_session = args.save_session
    #multisave = args.multisave

    # Adjust input, resolve conflicts
    autodetect_lig = 0 # default OFF
    #if reverse_search == 1:
    #    autodetect_lig = 1
    lig_scan_radius = str(lig_scan_radius)      # needs to be str
    intrfc_lig_radius = str(intrfc_lig_radius)  # needs to be str
    hoh_scan_radius = str(hoh_scan_radius)      # needs to be str
    #cndtlig_scan_radius = lig_scan_radius      # why is this "not used", since it is required later on? -local vrbl
    #broad_search_mode = False # previously called "reverse_mode"


    # Set directories, create job_id
    path_root = workdir
    pathSTRUCTS = path_root + '/structures'    # Directory with ALL pdb structures (used for fetch/download)
    pathLIGS = path_root + '/ligands'          # Directory with ALL pdb ligands (used for fetch/download)
    pathQRS = path_root + '/queries'           # Directory/index with parameters of previously run jobs
    pathXML = path_root + '/rsd_mappings'      # Directory with xml files for residue mappings

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
    path_results_structs = path_results + '/structure_files'  # Create subdirectory to store mmCIF files
    if not os.path.isdir(path_results_structs):
        os.makedirs(path_results_structs) 
    job_id = os.path.basename(os.path.normpath(path_results))


    if data is None:
        data = load_precompiled_data(workdir)
    dict_SIFTS = data.dict_SIFTS
    dict_rSIFTS = data.dict_rSIFTS

    # Get additional info
    # script_name = os.path.basename(__file__)    #log_file = script_name[:-3] + '_rejected_res_' + infile1[:-4] + '.log'
    # log_file_dnld = script_name + '_downloadErrors.log' #log_file_dnld = job_id + '_' + script_name + '_downloadErrors' + '.log'
    # log_file_dnld = path_root + '/download_errors.log'


    print('PyMOL version: ', cmd.get_version()[0:3])

    # Create directories if they don't exist
    print('Setting up directories')

    if os.path.isdir(pathSTRUCTS):
        print('Structure directory:', pathSTRUCTS)
    else:
        print('Creating structure directory:\t', pathSTRUCTS)
        os.makedirs(pathSTRUCTS)
    if os.path.isdir(pathLIGS):
        print('Ligands directory:\t', pathLIGS)
    else:
        print('Creating ligands directory:\t', pathLIGS)
        os.makedirs(pathLIGS)
    if os.path.isdir(pathQRS):
        print('Queries directory:\t', pathQRS)
    else:
        print('Creating queries directory:\t', pathQRS)
        os.makedirs(pathQRS)
    if os.path.isdir(pathXML):
        print('XML file directory:\t', pathXML)
    else:
        print('Creating XML file directory:\t', pathXML)
        os.makedirs(pathXML)
    print('Done\n')



    # Parse input
    print('Parsing input')

    try:
        #q = parse_query(query, autodetect_lig, water_as_ligand_auto, nonstd_rsds_as_lig_auto, d_aa_as_lig_auto)
        q = parse_query(query, autodetect_lig) # don't define, use function defaults
    except ValueError as e:
        print(e)
        wrong_input_error()

    user_chains = q.chains
    struct = q.struct
    ligand_names = q.ligands
    position = q.position
    autodetect_lig = q.autodetect_lig
    water_as_ligand_auto = q.water_as_ligand_auto
    nonstd_rsds_as_lig_auto = q.nonstd_rsds_as_lig_auto
    d_aa_as_lig_auto = q.d_aa_as_lig_auto



    # Parse chains
    if user_chains != 'ALL' and user_chains != '!':
        user_chains = ''.join(user_chains)
        user_chains = user_chains.split(',')
        #user_chains_bundle = '+'.join(user_chains)

        # Convert chains to structchain combos
        user_structchains = list()
        for user_chain in user_chains:
            user_structchain = struct.lower() + user_chain#.upper()
            user_structchains.append(user_structchain)


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
        except Exception as ex:
            print(ex)
            raise ValueError(f"Invalid ligands in query '{query}': use PDB ligand names")

    # Parse ligands
    if autodetect_lig == 1 and ligand_names is not None or autodetect_lig == 0:
        ligand_names = ''.join(ligand_names)
        ligand_names = ligand_names.split(',')
        ligand_names_bundle = '+'.join(ligand_names)


    # Print input info
    print('')
    print('Input structure:\t', struct)
    if user_chains == "!":
        print('Input chains:\t\t', 'LIG-BINDING ONLY')
    else:
        print('Input chains:\t\t', user_chains)
    if user_chains != 'ALL' and user_chains != '!':
        print('Input structchains:\t', user_structchains)
    if autodetect_lig == 1 and ligand_names is None:
        print('Input ligands:\t\t', 'auto-detect')
    elif autodetect_lig == 1 and ligand_names is not None: # this scenario should not occur anymore and thus be removed
        print('Input ligands:\t\t', ligand_names, '+ auto-detect')
    else:
        print('Input ligands:\t\t', ligand_names) #, '\t', ligand_names_bundle)
    if position is not None:
        print('Input position:\t\t', position)
    print('Done')

    # Toggle "!" chain switch on = examine only chains with defined ligands
    only_lig_chains = 0
    if user_chains == "!":
        only_lig_chains = 1


    # Finished parsing query and configuring settings, store final settings and define non-ligands here

    # Merge user settings and auto settings from parse query here
    if water_as_ligand_auto == 1 or water_as_ligand_usr == 1:
        water_as_ligand = 1
    else:
        water_as_ligand = 0
    if nonstd_rsds_as_lig_auto == 1 or nonstd_rsds_as_lig_usr == 1:
        nonstd_rsds_as_lig = 1
    else:
        nonstd_rsds_as_lig = 0
    if d_aa_as_lig_auto == 1 or d_aa_as_lig_usr == 1:
        d_aa_as_lig = 1
    else:
        d_aa_as_lig = 0

    # Define non-ligands (3-letter names of amino acids and h2o)
    nolig_resn = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    std_rsds = list(nolig_resn)
    if water_as_ligand == 0:
        nolig_resn.append('HOH')
    # Non-standard residues
    nonstd_rsds = "SEP TPO PSU MSE MSO 1MA 2MG 5MC 5MU 7MG H2U M2G OMC OMG PSU YG PYG PYL SEC PHA".split()
    if nonstd_rsds_as_lig == 0:
        nolig_resn.extend(nonstd_rsds)
    # D-amino acids
    d_aminoacids = "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA".split()
    if d_aa_as_lig == 0:
        nolig_resn.extend(d_aminoacids)

    # Pass final settings to a string
    settings_str = 'res' + str(res_threshold) + '_NMR' + str(include_nmr) + '_xrayonly' + str(xray_only) + '_ligfree' + str(lig_free_sites) + '_autodtctlig' + str(autodetect_lig) + '_h2olig' + str(water_as_ligand) + '_overlap' + str(overlap_threshold) + '_ligrad' + str(lig_scan_radius) + '_tmscore' + str(min_tmscore) + '_nonstdaas' + str(nonstd_rsds_as_lig) + '_daas' + str(d_aa_as_lig)

    '''# Test print
    print('\nuser chains, struct, ligand names, position, autodetect lig, water as ligand, nonstd rsds, d_aa_as_lig:\n', 
          user_chains, struct, ligand_names, position, autodetect_lig, water_as_ligand, nonstd_rsds_as_lig, d_aa_as_lig)
    print('\nnolig resis:', nolig_resn)
    #sys.exit(1)'''


    ## Find Apo candidates (rSIFTS)

    # Find & VERIFY input chains by UniProt ID (if they don't exist in uniprot, we cannot process them)
    # allow non-UniProt chains, because ligands can be assigned to non-protein chains
    print(f'\nFinding & verifying query chains ["{user_chains}"] by UniProt ID')
    discarded_chains = dict()   # Discarded chains (format: structchain + '\t' + discard_msg)
    usr_structchains_unverified = list()
    qr_uniprot_ids = dict()


    # Traceback "ALL" / "*" / "!" chains from SIFTS file
    if user_chains == 'ALL' or user_chains == "!":
        user_chains = list()
        user_structchains = list()

        #print('Considering all chains in query structure, finding chains..')
        for key in dict_SIFTS:
            if key[:4] == struct:
                user_chains.append(key[4:])
                user_structchains.append(key)
                print(key, dict_SIFTS[key])
                qr_uniprot_ids[key] = dict_SIFTS[key]
    else:
        for user_structchain in user_structchains[:]:
            try:
                print(user_structchain, dict_SIFTS[user_structchain])
                qr_uniprot_ids[user_structchain] = dict_SIFTS[user_structchain]
            except:
                user_structchains.remove(user_structchain)
                #discarded_chains.append(user_structchain + '\t' + 'Query chain not assigned UniProt ID\n')
                usr_structchains_unverified.append(user_structchain)

    # Put the UniProt found chains into a list
    query_unp_chains = list()
    for user_structchain in user_structchains:
        usr_chain = user_structchain[4:]
        query_unp_chains.append(usr_chain)


    # Map ligands with non-protein chains to protein chains [not allowed in broad search]
    non_protein_lig_chains = dict()
    if len(usr_structchains_unverified) > 0 and ligand_names is not None:
        print(f'-ligand assigned non-protein chain(s): {usr_structchains_unverified}, attempting to map to protein chains')

        for unverified_structchain in usr_structchains_unverified[:]:

            if position is not None:
                non_protein_lig_expression = struct + ' and chain ' + unverified_structchain[4:] + ' and resn ' + ligand_names_bundle + ' and resi ' + str(position)
            else:
                non_protein_lig_expression = struct + ' and chain ' + unverified_structchain[4:] + ' and resn ' + ligand_names_bundle

            cmd.reinitialize('everything')
            cmd.load(struct_path)
            s1n = cmd.select('non-protein_lig_' + ligand_names_bundle, non_protein_lig_expression)
            non_protein_lig_atoms = cmd.identify('non-protein_lig_' + ligand_names_bundle)
            non_protein_lig_atoms = '+'.join(str(i) for i in non_protein_lig_atoms)
            s2n = 'ID ' + non_protein_lig_atoms

            try: # Catch wrong chain exception (wrong case or non-existing chains)
                s3n = cmd.select('around_non-protein_lig' + unverified_structchain, 'polymer.protein near_to ' + intrfc_lig_radius + ' of ' + s2n)
            except Exception:
                wrong_input_error()

            if s1n != 0 and s3n != 0:
                non_protein_lig_chains[unverified_structchain] = cmd.get_chains('around_non-protein_lig' + unverified_structchain)
                print('found protein binding chains:', non_protein_lig_chains[unverified_structchain])
                print(f'Replacing non-protein chain [{unverified_structchain[4:]}] with binding chain {non_protein_lig_chains[unverified_structchain]}')
                for found_chain in non_protein_lig_chains[unverified_structchain]:
                    new_structchain = struct + found_chain  # Convert detected chains into structchains

                    # Verify new structchain with SIFTS before saving it
                    try:
                        print(new_structchain, dict_SIFTS[new_structchain])
                        user_structchains.append(new_structchain)
                    except Exception:
                        print('-remapped chain not found in SIFTS', found_chain)

                #user_structchains.remove(unverified_structchain) # already removed
            else:
                print('User specified chain does not exist in UniProt, removing it from input:\t', unverified_structchain)
                usr_structchains_unverified.remove(unverified_structchain)
                #user_chains.remove(unverified_structchain[4:])
                #user_structchains.remove(unverified_structchain)
                #discarded_chains.append(unverified_structchain + '\t' + 'No assigned UniProt ID\n')
                discard_msg = 'Query chain: No assigned UniProt ID (unverified)'
                discarded_chains.setdefault(unverified_structchain, []).append(discard_msg)
                usr_structchains_unverified.append(unverified_structchain)
            cmd.delete('not ' + struct)


    ## Handle"!" chain operator: Detect ligand chains and discard non-ligand query chains
    # If "!" operator used and ligand is defined but not detected, exit
    # If "!" and no ligands defined, detect and consider only chains with ligands #TODO
    if only_lig_chains == 1 and ligand_names is not None:# and len(user_structchains) > 1: 
        # Run this even for 1 chain and exit if defined ligand(s) not present
        print('\nFinding & removing non ligand-binding query chains')

        #if len(user_structchains) == 1:    print('Only UniProt chain is available, skipping process..\nDone')
        lig_bndng_prot_chains = dict()        
        qr_lig_positions = list()
        good_qr_chains = list()     # lig binding chains (that are also uniprot chains) priority is givent to assigned chain if it's a UNP chain, otherwise, all other binding chains are kept
        usr_structchains_non_ligbndng = list()

        if struct in cmd.get_object_list('all'):  # object names
            cmd.delete('not ' + struct)
        else:
            cmd.reinitialize('everything')
            cmd.load(struct_path)

        # Select defined ligands in whole structure
        all_struct_ligs_expr = struct + ' and resn ' + ligand_names_bundle
        s1lc_as = cmd.select('all_ligs_' + ligand_names_bundle, all_struct_ligs_expr)

        if s1lc_as > 0:
            myspace_qr_lig_positions = {'qr_ligs': []}
            cmd.iterate('all_ligs_' + ligand_names_bundle, 'qr_ligs.append(chain+"_"+resn+"_"+resi)', space=myspace_qr_lig_positions)

            # Remove duplicate values from query_lig_positions & Transfer binding sites to global dictionary
            for key, values in myspace_qr_lig_positions.items():
                qr_lig_positions = list(myspace_qr_lig_positions.fromkeys(values))

            # Build dict with ligands' assigned chains
            for lig_position in qr_lig_positions:
                qr_lig_chain = lig_position.split('_')[0]
                qr_lig_name = lig_position.split('_')[1]
                qr_lig_pos = lig_position.split('_')[2]

                # Find binding chains and build dict
                lig_sele_expr = struct + ' and chain ' + qr_lig_chain + ' and resn ' + qr_lig_name + ' and resi ' + qr_lig_pos
                s1lc_b = cmd.select('around_lig_' + lig_position, 'polymer.protein near_to ' + intrfc_lig_radius + ' of ' + lig_sele_expr)
                if s1lc_b > 0:
                    lig_bndng_prot_chains[lig_position] = cmd.get_chains('around_lig_' + lig_position)
                    #binding_chains.update(cmd.get_chains('around_lig_' + lig_position))

            # If none of the query ligands have (protein) binding residues, exit
            if len(lig_bndng_prot_chains) == 0:
                print(f'\nProtein binding residues not detected for query ligand(s) {ligand_names} on any query chain. Exiting..')
                sys.exit(1)

            #print(f'All query lig positions\n{qr_lig_positions}')
            #print_dict_readable(lig_bndng_prot_chains, 'Ligand assigned (left, one per ligand) and binding protein chains (right, may be multiple)')

            # Build set with nessessary UNP chains (per ligand approach)
            #print('\nBuilding set of good binding chains')
            for ligand, binding_chains in lig_bndng_prot_chains.items():
                assigned_chain = ligand.split('_')[0]
                if assigned_chain in (binding_chains and query_unp_chains):
                    #print(f'-assigned chain [{assigned_chain}] for ligand [{ligand}] present in {binding_chains} AND in {query_unp_chains}\n*adding to good chains')
                    good_qr_chains.append(assigned_chain)
                else:
                    #print(f'-assigned chain [{assigned_chain}] for ligand [{ligand}] NOT present in {binding_chains} AND in {query_unp_chains}\n*looking for other binding UNP chains')
                    if len(binding_chains) > 0:
                        #print(f'*adding binding chain(s) {binding_chains} to good chains')
                        good_qr_chains.extend(binding_chains)

            # Remove verified query structchains that are not in good chains (binding the ligands)
            # Note: if ligand binds multiple protein chains, and assigned chain is among them, keep only assigned chain -else keep any valid UNP chain
            for user_structchain in user_structchains[:]:
                usr_chain = user_structchain[4:]                    
                if usr_chain not in good_qr_chains:
                    #print(f'\nRemoving chain {usr_chain}')
                    usr_structchains_non_ligbndng.append(user_structchain)
                    user_structchains.remove(user_structchain)
                    #discarded_chains.append(user_structchain + '\t' + 'Query chain: Non ligand-binding chain ("!")\n')
                    discard_msg = 'Query chain: Non ligand-binding chain ("!")'
                    discarded_chains.setdefault(user_structchain, []).append(discard_msg)
                    #print(f'[{usr_chain}] not in assigned chains[{assigned_chains}] or binding chains[{binding_chains}]. Removing..')

            cmd.delete('not ' + struct)
            print(f'Approved binding chains:\t{good_qr_chains}')

        else: # Query ligand not detected (exit?)
            print(f'\nQuery ligand(s) {ligand_names} not detected on any query chain. Exiting..')
            sys.exit(1)


    # Print report of chain verification
    print('\nInput chains verified:\t\t', user_structchains)#, user_chains)
    print('Input chains unverified:\t', usr_structchains_unverified)
    if len(non_protein_lig_chains) > 0:
        print('*Remapped chains:\t\t\t', non_protein_lig_chains)  # TODO unverified chains are not discarded
    if only_lig_chains == 1: #len(usr_structchains_non_ligbndng) > 0:
        try:
            print('Removed non-binding chains:\t', usr_structchains_non_ligbndng)
        except Exception:
            pass

    #sys.exit(1)
    # Capture the full query expression - for query indexing and searching previous jobs
    if autodetect_lig == 0:
        user_query_parameters = struct + '_' + ','.join(user_chains) + '_' + ','.join(ligand_names)
    elif autodetect_lig == 1 and ligand_names is not None:  # this probably does not occur anymore
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


    # Get apo candidates from rSIFTS dict, calculate sequence overlap (observed residues)
    print(f'Calculating overall UniProt mappings for query chains\n{user_structchains}')

    dictApoCandidates = dict()      # candidate structchains over user-threshold UNP overlap
    dictApoCandidates_b = dict()    # candidate structchains with non-zero UNP overlap (more)
    uniprot_overlap_all = dict()    # overlap of all UniProt chains
    uniprot_overlap_invrs = dict()  # same but keys are query chains for easier manipulation
    #uniprot_overlap = dict()       # overlap of successful hits over user threshold
    uniprot_segments_dict = dict()  # the segments of the observed residues in chains of interest only
    user_structchains_unp = dict()
    own_chains = dict()

    for user_structchain in user_structchains:

        query_chain_total_unp_length = 0 # reset total length
        query_chain_segments = list()
        cndt_chains_segments = dict()

        # Get UniProt ID of query chain to build rsd mappings dict
        uniprot_id = dict_SIFTS[user_structchain]  # Get primary UniProt AC
        user_structchains_unp[user_structchain] = uniprot_id  # Pass it as value to query structchain
        #print(user_structchain, uniprot_id)

        # Get all chains with their segments for given UNP AC (i.e., list of values)
        all_segments_for_UNP_AC = dict_rSIFTS[uniprot_id]  # list
        #print('\nall_segments_for_UNP_AC\n', all_segments_for_UNP_AC)

        # Put all query segments together and all candidate segments together
        for segment_overlap in all_segments_for_UNP_AC:
            #current_struct = segment_overlap.split()[0][:4]
            #print(current_struct, user_structchain[:4])
            if segment_overlap.split()[0][:4] == user_structchain[:4]: # same struct
                if segment_overlap.split()[0] == user_structchain: # same structchain
                    #query_chain_segments.append(segment_overlap.split()[1:3])  # append query chain segment start/end
                    query_chain_segments.append(segment_overlap.split()[1] + ' ' + segment_overlap.split()[2])  # append query chain segment start/end
                else: # same struct, different chain ("own_chain") #if segment_overlap.split()[0][:4] == user_structchain[:4] and segment_overlap.split()[0] != user_structchain:
                    #own_chains.append(segment_overlap.split()[0]) # for list
                    own_chains.setdefault(user_structchain, []).append(segment_overlap.split()[0])
            else:
                cndt_chains_segments.setdefault(segment_overlap.split()[0], []).append(segment_overlap.split()[1] + ' ' + segment_overlap.split()[2])

        #print('\nquery_chain_segments\n', query_chain_segments)
        #print_dict_readable(cndt_chains_segments, '\ncndt_chains_segments dict')


        # Calculate total length of query chain segments
        for query_chain_segment in query_chain_segments:
            seg_start = query_chain_segment.split()[0]
            seg_end = query_chain_segment.split()[1]
            seg_l = int(seg_end) - int(seg_start)
            query_chain_total_unp_length += seg_l
        #print('\ntotal query chain segment length:', query_chain_total_unp_length, '\n')


        # Calculate percentage of coverage of cndt chains (segments) on total query unp length
        for cndt_chain, cndt_segments in cndt_chains_segments.items():
            #cndt_total_unp_legnth = 0
            cndt_seg_prcnt = 0
            for cndt_segment in cndt_segments:
                cndt_unp_legnth = 0
                cndt_seg_start = cndt_segment.split()[0]
                cndt_seg_end = cndt_segment.split()[1]
                cndt_unp_legnth = int(cndt_seg_end) - int(cndt_seg_start)

                for query_chain_segment in query_chain_segments:
                    seg_start = query_chain_segment.split()[0]
                    seg_end = query_chain_segment.split()[1]

                    # Calculate overlap and % of overlap on query (x1,x2)
                    result =  min(int(seg_end), int(cndt_seg_end)) - max(int(seg_start), int(cndt_seg_start))

                    if result > 0: # ignore neg overlap
                        cndt_seg_prcnt += result / query_chain_total_unp_length * 100
                        percent = result / query_chain_total_unp_length * 100
                        percent = round(percent, 1)
                    else:
                        percent = 0
                    #print(cndt_chain, cndt_segment, '->', user_structchain, query_chain_segment, '[overlap length, seg_percent, total_prcnt]:', result, round(percent), round(cndt_seg_prcnt, 2))
                    #dictApoCandidates.setdefault(user_structchain+' '+query_chain_segment, []).append(cndt_chain+' '+cndt_segment+' '+str(cndt_unp_legnth)+' '+str(result)+' '+str(percent))
                    #dictApoCandidates.setdefault(user_structchain, []).append(cndt_chain+' '+cndt_segment+' '+str(cndt_unp_legnth))
                    #dictApoCandidates_b.setdefault(user_structchain, []).append(cndt_chain+' '+cndt_segment+' '+str(cndt_unp_legnth))
                    uniprot_segments_dict.setdefault(user_structchain, []).append(cndt_chain+' '+cndt_segment+' '+str(cndt_unp_legnth))

                    # Build dict with calculated overlap
                    uniprot_overlap_all.setdefault(cndt_chain, []).append(user_structchain + ' ' + str(percent)) # this keeps all calculated overlaps (from candidate to all query chains)
                    uniprot_overlap_invrs.setdefault(user_structchain,  []).append(cndt_chain + ' ' + str(percent))  # arranged by query structchains
        #sys.exit(1)


    # Merge calculated UniProt overlap percentages of same chains into a single percentage
    uniprot_overlap_merged = merge_fragmented_unp_overlaps(uniprot_overlap_all)
    uniprot_overlap_invrs_merged = merge_fragmented_unp_overlaps(uniprot_overlap_invrs)


    # Remove query chains with no candidates (or candidates with zero overlap [in Uniprot overlap merged])
    # Build candidate dicts here
    # Discard: i)candidates with zero overlap, ii)candidates below overlap threshold (not completely) iii)query chains with no candidates
    #if overlap_threshold != 0 and percent >= overlap_threshold or overlap_threshold == 0 and percent > 0:

    for qr_structchain, cndt_structchains_w_overlaps in uniprot_overlap_invrs_merged.items():
        for cndt_structchain_w_overlap in cndt_structchains_w_overlaps:
            cndt_structchain = cndt_structchain_w_overlap.split()[0]
            cndt_overlap = cndt_structchain_w_overlap.split()[1]
            if float(cndt_overlap) > 0: # keep any positive candidate (also removes query chains without candidates)
                dictApoCandidates_b.setdefault(qr_structchain, []).append(cndt_structchain_w_overlap)

                if float(cndt_overlap) >= overlap_threshold: # only candidates above user threshold [used for apo query chains]
                    dictApoCandidates.setdefault(qr_structchain, []).append(cndt_structchain_w_overlap)

            else:  # Discard
                #discarded_chains.append(cndt_structchain + '\t' + '"0" UNP overlap in observed residues\n')
                discard_msg = '"0" UNP overlap in observed residues'
                discarded_chains.setdefault(cndt_structchain, []).append(discard_msg)

    #print_dict_readable(uniprot_segments_dict, '\nUniprot_segments_dict')
    #print_dict_readable(uniprot_overlap_merged, '\nUniprot overlap merged')
    #print_dict_readable(uniprot_overlap_invrs, '\nUniprot overlap Inverse')
    #print_dict_readable(uniprot_overlap_invrs_merged, '\nUniprot overlap Inverse merged')
    #print_dict_readable(uniprot_overlap_all, '\nUniprot overlap all')
    #print_dict_readable(user_structchains_unp, '\nuser_structchains_unp')

    #print_dict_readable(qr_uniprot_ids, '\nqr_uniprot_ids')
    #print_dict_readable(dictApoCandidates, '\ndictApoCandidates')
    #print_dict_readable(dictApoCandidates_b, 'dictApoCandidates_b (>0% overlap)') # Test print
    #print_dict_readable(own_chains, '\nOwn chains dict')


    # Print info for detected candidate chains
    for qr_structchain, cndt_structchains_w_overlaps in dictApoCandidates_b.items():

        #if own_chains.get(qr_structchain) is not None:
            #print(f'\nTotal chains for {qr_structchain}:{qr_uniprot_ids[qr_structchain]} (including/excluding query structure): [{len(dict_rSIFTS[qr_uniprot_ids[qr_structchain]])}]/[{len(dict_rSIFTS[uniprot_id]) - len(own_chains[qr_structchain])}]')
            #print(f'\nMax candidate chains for {qr_structchain}:{qr_uniprot_ids[qr_structchain]} (excluding/including query structure): [{len(dictApoCandidates_b[qr_structchain])}]/[{len(dictApoCandidates_b[qr_structchain]) + len(own_chains[qr_structchain])}]')
        #else:
            #print(f'\nMax candidate chains for {qr_structchain}:{qr_uniprot_ids[qr_structchain]} (excluding/including query structure): [{len(dictApoCandidates_b[qr_structchain])}]/[{len(dictApoCandidates_b[qr_structchain]) + 0}]')

        try:
            print(f'Candidate chains over [0%] overlap threshold for [{qr_structchain}]:\t{len(dictApoCandidates_b[qr_structchain])}')
        except Exception as ex:
            print(f'Candidate chains over [0%] overlap threshold for [{qr_structchain}]:\t{dictApoCandidates_b.get(qr_structchain)},\t{ex}')
        if overlap_threshold != 0: # otherwise no need to re-print
            try:
                print(f'Candidate chains over [{overlap_threshold}%] (user-specified) overlap threshold for [{qr_structchain}]:\t{len(dictApoCandidates[qr_structchain])}') # - {dictApoCandidates[dict_key]}') #[dict_key][0].split()[0]}]')
            except Exception as ex:
                print(f'Candidate chains over [{overlap_threshold}%] (user-specified) overlap threshold for [{qr_structchain}]:\t{dictApoCandidates.get(qr_structchain)},\t{ex}') # - {dictApoCandidates[dict_key]}') #[dict_key][0].split()[0]}]')


    # Calculate num of total candidate chains for all query chains (use larger dict (dictApoCandidates_b))
    total_chains = sum([len(dictApoCandidates_b[x]) for x in dictApoCandidates_b if isinstance(dictApoCandidates_b[x], list)])
    #total_chains = sum([len(dictApoCandidates[x]) for x in dictApoCandidates if isinstance(dictApoCandidates[x], list)])
    #print(f'\nTotal (non-unique) candidate chains over user-specified overlap threshold [{overlap_threshold}%]:\t{total_chains}\n') # this if for over user threshold
    print(f'\nTotal (non-unique) candidate chains over [0%] overlap:\t{total_chains}')

    # End script if there are 0 candidate chains OR no UniProt is assigned to the query structure
    if total_chains == 0:
        print('\n=== Ending program ===')
        sys.exit(1)

    '''
    # Remove duplicate values from dictApoCandidates & dictApoCandidates_b
    for key, values in dictApoCandidates.items():
        dictApoCandidates[key] = list(dictApoCandidates.fromkeys(values))   # preserves the order of values

    for key, values in dictApoCandidates_b.items():
        dictApoCandidates_b[key] = list(dictApoCandidates_b.fromkeys(values))   # preserves the order of values
    '''


    ## Candidate structure evaluation

    def structures_from_structchains(candidate_dict_in):  # (dict in, set out) function to create a set of PDB codes to download from a dict of structchains
        list_out = list()
        for key, values in candidate_dict_in.items():
            for i in values:    # iterate over the values in each key, i = struct/chain combo & SP_BEG SP_END etc
                candidate_struct = i.split()[0][:4]  # split value strings to get structure only
                list_out.append(candidate_struct)
        return list_out  #list(dict.fromkeys(list_out))  # remove redundancy & preserve order

    # Put all structures for downloading into set [include rsd mapping candidates]
    apo_candidate_structs1 = structures_from_structchains(dictApoCandidates)
    apo_candidate_structs2 = structures_from_structchains(dictApoCandidates_b)  # Add the rest UniProt structures
    #print(apo_candidate_structs1)
    #print(apo_candidate_structs2)
    apo_candidate_structs = apo_candidate_structs1 + apo_candidate_structs2
    apo_candidate_structs = list(dict.fromkeys(apo_candidate_structs))  # remove redundancy & preserve order
    #print(apo_candidate_structs)

    print('Total structures to be parsed: ', len(apo_candidate_structs), '\n')


    # Download candidate structures to specified directory [this should be replaced by load]
    for apo_candidate_structure in apo_candidate_structs[:]:
        try:
            download_mmCIF_gz2(apo_candidate_structure, pathSTRUCTS)
        except Exception:
            # Don't fail, instead remove candidate structure from queue
            #discarded_chains.append(apo_candidate_structure + '\t' + 'PDB structure not found\n')
            discard_msg = 'PDB structure not found'
            discarded_chains.setdefault(apo_candidate_structure, []).append(discard_msg)
            apo_candidate_structs.remove(apo_candidate_structure)
            print(f'*file {apo_candidate_structure} not found, removing from candindates list')


    # Parse (mmCIF) structures to get resolution & method. Apply cut-offs
    resolution_dict = dict()
    exp_method_dict = dict()
    r_free_dict = dict()
    print('Checking resolution and experimental method of candidate structures')
    for apo_candidate_struct in apo_candidate_structs:
        apo_candidate_structPath = download_mmCIF_gz2(apo_candidate_struct, pathSTRUCTS)
        resolution = '?\t'
        # Check if resolution is already measured:
        #resolution = resolution_dict.get(apo_candidate_struct)
        #method = exp_method_dict.get(apo_candidate_struct)
        #if resolution is not None and method is not None:
            #resolution = float(resolution)
        #else:
        with gzip.open (apo_candidate_structPath, 'rt') as mmCIFin:
            method = '?' # reset method/res/r free
            resolution = '-'
            r_free = '-'
            for line in mmCIFin:
                try:
                    if method == '?':
                        if line.split()[0] == '_exptl.method' and line.split("'")[1] is not None:
                            method = line.split("'")[1]  # capture experimental method #method = ' '.join(line.split()[1:])
                            exp_method_dict[apo_candidate_struct] = method
                            if method == 'SOLUTION NMR':  # break fast if method is 'NMR'
                                method = '\t ' + method + '\t\t'
                                break
                            elif method == 'SOLID-STATE NMR':
                                method = '\t ' + method + '\t'
                                break
                            continue
                        if method == '?' and line.split()[0] == '_refine.pdbx_refine_id': # second attempt to get method [+EPR]
                            method = line.split("'")[1] + '*'
                            exp_method_dict[apo_candidate_struct] = method
                        if method == '?' and line.split()[0] == '_refine_hist.pdbx_refine_id': # third attempt to get method [+neutron diffraction]
                            method = line.split("'")[1] + '**'
                            exp_method_dict[apo_candidate_struct] = method

                    elif line.split()[0] == '_refine.ls_d_res_high' and float(line.split()[1]):
                        resolution = round(float(line.split()[1]), 2)  # X-ray highest resolution
                        resolution_dict[apo_candidate_struct] = str(resolution)
                    elif line.split()[0] == '_refine_hist.d_res_high' and float(line.split()[1]): # X-ray res in neutron scattering methods
                        resolution = round(float(line.split()[1]), 2)  # X-ray highest resolution in neutron diffraction
                        resolution_dict[apo_candidate_struct] = str(resolution)
                    elif line.split()[0] == '_refine.ls_R_factor_R_free':# and float(line.split()[1]):
                        if float(line.split()[1]):
                            r_free = round(float(line.split()[1]), 3)
                            r_free_dict[apo_candidate_struct] = str(r_free)
                        break
                    elif line.split()[0] == '_em_3d_reconstruction.resolution' and float(line.split()[1]):
                        resolution = round(float(line.split()[1]), 2)  # EM resolution
                        resolution_dict[apo_candidate_struct] = str(resolution)
                        break
                except Exception:# as ex: # getting weird but harmless exceptions
                    #print('Problem parsing structure: ', apo_candidate_struct)
                    pass

            #resolution = float(resolution_dict[apo_candidate_struct])
            #method = exp_method_dict[apo_candidate_struct]
            try:
                if include_nmr == 1 and method == '\t SOLUTION NMR\t\t' or method == '\t SOLID-STATE NMR\t' and xray_only == 0 or xray_only == 1 and method == 'X-RAY DIFFRACTION' and resolution <= res_threshold or xray_only == 0 and resolution <= res_threshold:
                    print(apo_candidate_struct, ' resolution:\t', resolution, '\t', method, ' \tPASS')  # Xray/EM
                else:
                    #discarded_chains.append(apo_candidate_struct + '\t' + 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']\n')
                    discard_msg = 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']'
                    discarded_chains.setdefault(apo_candidate_struct, []).append(discard_msg)
                    print(apo_candidate_struct, ' resolution:\t', resolution, '\t', method, ' \tFAIL')  # Xray/EM
            except Exception as ex:
                #discarded_chains.append(apo_candidate_struct + '\t' + 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']\n')
                discard_msg = 'Resolution/exp. method\t[' + str(resolution) + ' ' + method + ']'
                discarded_chains.setdefault(apo_candidate_struct, []).append(discard_msg)
                print('*Exception', ex, '\n', apo_candidate_struct, ' resolution:\t', resolution, '\t', method, '\t\tFAIL')  # NMR

    print('Done')


    # Compile set of structures to be discarded (remove comments)
    discard_structs = set()
    for i in discarded_chains:
        if len(i.split()[0]) == 4:
            discard_structs.add(i.split()[0])

    # Discard candidates below threshold(s). Put remainder into new dict, discard UniProt numbering
    print('\nStructures to discard:\t', len(discard_structs))
    if len(discard_structs) > 0:
        print('Discarding structures\t', discard_structs)


    def remove_low_res_structs(candidates_dict, low_res_structs_set, cleave=0, chain_limit=apo_chain_limit): # return dict with good candidates only
        dict_out = dict()
        for key, values in candidates_dict.items():
            for i in values[:chain_limit]:
                if i.split()[0][:4] in low_res_structs_set:
                    print(f'removing candidate chain [{i.split()[0]}] from query chain [{key.split()[0]}]')
                else:
                    if cleave == 0:
                        dict_out.setdefault(key.split()[0], []).append(i)
                    else:
                        dict_out.setdefault(key.split()[0], []).append(i.split()[0]) # same dict without calculated UNP coverage or candidate percentage
        return dict_out

    # Discard candidate structures below threshold(s). Put remainder into new dict
    dictApoCandidates_1 = remove_low_res_structs(dictApoCandidates, discard_structs, cleave=1) # cleave dict (discard UniProt numbering)
    dictApoCandidates_b1 = remove_low_res_structs(dictApoCandidates_b, discard_structs)

    #print_dict_readable(dictApoCandidates_1, '\ndictApoCandidates_1')
    #print_dict_readable(dictApoCandidates_b1, '\ndictApoCandidates_b1')


    eligible_chains1 = sum([len(dictApoCandidates_1[x]) for x in dictApoCandidates_1 if isinstance(dictApoCandidates_1[x], list)])
    eligible_chainsb1 = sum([len(dictApoCandidates_b1[x]) for x in dictApoCandidates_b1 if isinstance(dictApoCandidates_b1[x], list)])

    print(f'Candidate chains (over 0% UNP overlap - for UNP residue mapping) satisfying structure quality requirements (method/resolution) [{res_threshold} Å]:\t{eligible_chainsb1}')
    if overlap_threshold != 0: # otherwise no need to re-print
        print(f'Candidate chains (over user-specified [{overlap_threshold}%] UNP overlap - for apo query chains) satisfying structure quality requirements (method/resolution) [{res_threshold} Å]:\t{eligible_chains1}')
    #print_dict_readable(exp_method_dict, '\nexp_method_dict')
    #sys.exit(1)

    # Make dict with query struct:chains
    dictQueryChains = dict()
    for key in dictApoCandidates_1.keys():
        qr_struct = key[:4]
        qr_chain = key[4:]
        dictQueryChains.setdefault(qr_struct, []).append(qr_chain)



    # Configure autodetect lig search expression according to variable
    if autodetect_lig == 1:
        autodetect_lig_expression = ' or (hetatm and not solvent and not polymer)'
    else:
        autodetect_lig_expression = ''


    #######################################
    # Define query ligand search expression


    if position is not None: # (4 args) assumes that everything is specified (chains(taken care of), 1 ligand, position) ignore_autodetect_lig
        
        if ligand_names[0] == 'HOH': # mark selection & change lig scan radius -handled later
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
            else:
                search_term = 'hetatm and not solvent near_to ' + lig_scan_radius + ' of (resi ' + position + ' and resn ' + ligand_names_bundle + ')'

        else: # find ligands # normal ligand? Error: ligand has different chain
            #print('\nUnaccounted-for query selection case, using default search\n') # TODO examine several examples
            print('\n', ligand_names, query)
            search_term = '(resn ' + ligand_names_bundle + ' and resi ' + position + autodetect_lig_expression + ')'

    else: # position == None (3 or less args)

        if ligand_names is not None: # ligands specified

            if ligand_names[0] in nonstd_rsds:
                search_term = '(resn ' + ligand_names_bundle + autodetect_lig_expression + ')'

            else: # assumes ligand is a real ligand
                search_term = '(hetatm and resn ' + ligand_names_bundle + autodetect_lig_expression + ')'

        elif ligand_names is None:# and autodetect_lig == 1: 
            search_term = 'hetatm and not solvent and not polymer' # water as ligand should be ignored here without a position


    print('\nSearch term = ', search_term)

    if autodetect_lig == 1:
        print('Search mode: Broad (auto-detecting query ligands)')
    else:
        print('Search mode: Focused')



    ##################  Query ligand detection  ##################

    query_lig_positions = dict()
    cndt_lig_positions = dict()
    query_chain_states = dict()
    apo_holo_dict = dict()
    apo_holo_dict_H = dict()
    dict_rsd_map_candidates = dict()
    bad_candidates_rsd_map = dict()
    penultimate_candidates = dict()
    final_candidates = dict()
    candidate_scores = dict()

    #qBS_centerofmass = dict()
    #qBS_coords = dict()

    # TODO: check if next line is correct and if "progress_processed_candidates += 1" is called everywhere is should (moved)
    #progress_total_candidates = sum([len(lst) for lst in dictApoCandidates_1.values()])
    #progress_total_candidates = sum([len(lst) for lst in final_candidates.values()])
    #progress_processed_candidates = 0

    def track_progress(write_results: bool = False):
        if args.track_progress:
            write_file(path_results + '/.progress', f"{progress_processed_candidates}/{progress_total_candidates}")
            if write_results:
                write_results_apo_csv(apo_holo_dict, path_results)
                write_results_holo_csv(apo_holo_dict_H, path_results)
                write_ligands_csv(query_lig_positions, cndt_lig_positions, path_results)
                write_discarded_chains(discarded_chains, path_results)


    ###########################################################################
    # Find interface ligands for query structure (including nucleic acid chains)
    # when the user-specified chain of a ligand is different than the actual PDB chain
    # but the ligand actually binds the specified chain

    # Allow when ligand is not residue or water and ALSO in broad search
    # Maybe make this conditional with a parameter (would like to avoid that)
    # Note: interface ligands in broad search are still missed if assigned non-protein chain

    #print(dictApoCandidates_1)
    #print(dictQueryChains)    

    search_interface = False

    if ligand_names is not None: # and if interface search allowed by parameter?
        for x in ligand_names:
            if x not in nolig_resn and x != 'HOH': # this maybe should be OR == 'HOH'
                search_interface = True
    elif ligand_names is None: # and if interface search allowed by parameter?
        search_interface = True

    interface_ligands = dict()
    interface_ligands_list = list()
    if search_interface:
        print(f'\n=== Searching query structure [{struct}] for interface ligands in chains {dictQueryChains[struct]} ===')

        query_struct = struct
        query_struct_path = download_mmCIF_gz2(query_struct, pathSTRUCTS)

        # Initialize PyMOL but don't reload Holo if present
        if query_struct in cmd.get_object_list('all'):  # object names
            cmd.delete('not ' + query_struct)
        else:
            cmd.reinitialize('everything')
            cmd.load(query_struct_path)

        # Configure search expression for interface ligands 
        if ligand_names is not None:
            if position is None:
                search_interface_expression = query_struct + ' and hetatm and not solvent and resn ' + ligand_names_bundle
            else:  # TODO we don't want all ligands
                search_interface_expression = query_struct + ' and hetatm and not solvent and resn ' + ligand_names_bundle + ' and resi ' + position
        else:
            search_interface_expression = query_struct + ' and hetatm and not solvent and not polymer'

        #print(struct, dictQueryChains[struct])
        all_ligands_selection = cmd.select('structure_ligands', search_interface_expression)

        if all_ligands_selection != 0:
            myspace_intrfc = {'all_ligs': []}
            cmd.iterate('structure_ligands', 'all_ligs.append( (resn+"_"+chain+"_"+resi,ID) )', space=myspace_intrfc)

            # Transfer list of ligand positions to dict
            all_qr_ligands = dict()
            for i in myspace_intrfc['all_ligs']:

                qr_lig_position = i[0]
                qr_lig_atomid = i[1]
                all_qr_ligands.setdefault(qr_lig_position, []).append(str(qr_lig_atomid))
            #print(all_qr_ligands_dict)

            all_qr_ligands_chains = dict()
            for intrfc_position, atom_ids in all_qr_ligands.items():

                # Join ligand atoms
                atom_ids = '+'.join(atom_ids)

                # Select atoms around ligand atoms
                s1 = 'ID ' + atom_ids
                around_all_ligands_sele = cmd.select('around_' + intrfc_position, 'polymer near_to ' + intrfc_lig_radius + ' of ' + s1)
                if around_all_ligands_sele == 0:
                    continue
                else:
                    all_qr_ligands_chains[intrfc_position] = cmd.get_chains('around_' + intrfc_position)
                #print(intrfc_position, cmd.get_chains('around_' + intrfc_position))
            #print('all_qr_ligands_chains DICT', all_qr_ligands_chains)

            # Find interface ligands in query chains
            for intrfc_position, chains in all_qr_ligands_chains.items():
                if len(chains) > 1:
                    for chain in chains:
                        if chain in user_chains: # this can find all interface ligands (that have a common chain with user_chains) at once, faster but has to be moved upstream
                        #if query_chain in chain:

                            # Keep only relevant interface ligands
                            #if lig_chain not in query_chain: # this would be detected anyway, we want interface ligands with non-user chains
                            interface_ligands[intrfc_position] = chains
                            interface_ligands_list.append(intrfc_position.replace("_", " "))
                            #print(f'Interface ligand detected [{intrfc_position.replace("_", " ")}] added to query selection')

            if len(interface_ligands) == 0:
                print('No interface ligands found')
            else:
                print('Interface ligand(s) detected around query chains:', interface_ligands_list)
                print('On chains:', interface_ligands)

    # end interface ligand search loop


    ### Start query chain loop

    for query_structchain, candidates_structchains in dictApoCandidates_1.items():
        print('')
        print(f'=== Processing query chain {query_structchain} ===')
        print('Finding query ligands')

        query_struct = query_structchain[:4]
        query_chain = query_structchain[4:]
        query_struct_path = download_mmCIF_gz2(query_struct, pathSTRUCTS)

        # Initialize PyMOL but don't reload Query if present
        if query_struct in cmd.get_object_list('all'):  # object names
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
            if ligands_selection == 0:# and autodetect_lig == 0: # This should not occur anymore (unless the ligand belongs to non-protein chain)
                print('No ligands found in PDB chain')


        # Search & select interface ligands (if any)
        interface_ligands_selection = 0
        if len(interface_ligands) > 0:
            print('Looking for interface ligands on chain', query_structchain)
            #interface_ligands_list.append(intrfc_position.replace("_", " "))
            #for chains in interface_ligands[query_struct]:
            for intrfc_lig_position, chains in interface_ligands.items():
                for chain in chains:
                    if query_chain in chain: # TODO decide whether to remove ligand after using
                        print('Adding interface ligand to selection: ', intrfc_lig_position)

                        # PyMOL selection (to be done within loop)
                        lig_resn = intrfc_lig_position.split('_')[0]
                        lig_chain = intrfc_lig_position.split('_')[1]
                        lig_resi = intrfc_lig_position.split('_')[2]
                        interface_lig_selection = query_struct + ' and resn ' + lig_resn + ' and chain ' + lig_chain + ' and resi ' + lig_resi
                        interface_ligands_selection = cmd.select('query_ligands', interface_lig_selection, merge=1)


        # Annotate query chain as apo or holo or skip it
        if ligands_selection == 0 and interface_ligands_selection == 0 and autodetect_lig == 1:
                print('No ligands found in broad search - continuing search')
                query_chain_states[query_structchain] = 'apo'
        elif ligands_selection == 0 and interface_ligands_selection == 0 and autodetect_lig == 0:
            print('No ligands found in focused search mode - skipping chain')
            continue
        elif ligands_selection > 0 or interface_ligands_selection != 0:
            query_chain_states[query_structchain] = 'holo'


        # If query is holo, identify the ligands
        ligands_atoms = list()
        if query_chain_states[query_structchain] == 'holo':

            # Identify atom IDs of selected ligand atoms
            ligands_atoms = cmd.identify('query_ligands', mode=0)  # bulk atoms of all query ligands

            # Get positions of specified ligands (iterate through atoms)
            myspace = {'positions': []}
            cmd.iterate('query_ligands', 'positions.append(resi +" "+ chain +" "+ resn)', space = myspace) # iterates through atoms of selection

            # Transfer temporary list with positions to dict
            for i in myspace['positions']:
                query_lig_positions.setdefault(query_structchain, []).append(i)

            # Remove duplicate values from query_lig_positions
            for key, value in query_lig_positions.items():
                query_lig_positions[key] = list(query_lig_positions.fromkeys(value))   # preserves the order of values

            # Name query ligands as seperate selections per "residue"/position. Put real (detected) ligand names into set
            query_lig_names = set()
            '''
            for ligand in query_lig_positions[query_structchain]: # try merge this loop function with the one below
                resi = ligand.split()[0]
                chain = ligand.split()[1]
                resn = ligand.split()[2]
                ligand_ = ligand.replace(' ', '_')
                query_lig_names.add(resn)
                cmd.select('holo_' + ligand_, query_struct + '& resi ' + resi + '& chain ' + chain + '& resn ' + resn) # s1
            '''
            #print_dict_readable(query_lig_positions, '\nQuery lig positions dict')

            # Find binding residues (protein only) for all query chains
            binding_res_dict = dict()

            for structchain_i, ligands in query_lig_positions.items():
                for ligand in ligands:

                    resi = ligand.split()[0]
                    chain = ligand.split()[1]
                    resn = ligand.split()[2]

                    # Put real (detected) ligand names into set
                    if structchain_i == query_structchain:
                        query_lig_names.add(resn)

                    ligand = ligand.replace(' ', '_')

                    # Find ligand
                    cmd.select(ligand + '_' + struct, struct + ' and chain ' + chain + ' and resi ' + resi + ' and resn ' + resn)

                    # Find & select binding residues
                    s1 = struct + ' and polymer.protein'
                    s2 = ligand + '_' + struct
                    cmd.select('arnd_' + ligand, s1 + ' near_to ' + lig_scan_radius + ' of ' + s2)

                    # Iterate and identify binding residues
                    myspace_positions = {'binding_resis': []}
                    cmd.iterate('arnd_' + ligand, 'binding_resis.append(segi+"_"+chain+"_"+resn+"_"+resi)', space=myspace_positions)

                    # Remove duplicate values from query_lig_positions & Transfer binding sites to global dictionary
                    for key, value in myspace_positions.items():
                        binding_res_dict[ligand] = list(myspace_positions.fromkeys(value))

            # Unpack binding residues into a single list
            binding_res_unpacked = []
            for sublist in binding_res_dict.values():
                binding_res_unpacked.extend(sublist)
            binding_res_unpacked = list(dict.fromkeys(binding_res_unpacked)) # Remove duplicates

            # Pass auto-detected ligands as user ligands
            if autodetect_lig == 1 and ligand_names is None:
                ligand_names = query_lig_names.copy() # ligand_names = user-specified ligands. When no ligands specified, they might be undefined


        # Print universal ligand report for query chain
        print('\nQuery ligand information')
        print('Total atoms: ', len(ligands_atoms))
        print('Atom IDs: ', ligands_atoms)
        print('Position/chain/name: ', query_lig_positions.get(query_structchain))
        print('State:', query_chain_states[query_structchain])  #, '/', len(query_lig_positions.get(query_structchain)))
        #print('\nquery_lig_positions', query_lig_positions)


        # Print binding site report for query chain
        if query_chain_states[query_structchain] == 'holo':
            total1 = sum(len(v) for v in binding_res_dict.values())
            total2 = len(binding_res_unpacked)
            #print(binding_res_dict)
            print('\nQuery binding site information')
            #print(f'Parsing query structure [{struct}]')
            #print('Detecting binding residues of query chain (segment, chain, pdb residue number)')
            print(f'Binding residues clustered per ligand:\n{binding_res_dict}')
            #print(f'Bulk unique binding residues [segment, chain, residue, position]:\n{binding_res_unpacked}')
            print(f'Total/unique binding residues: [{total1}]/[{total2}]') # TODO these seem to be cummulative binding residues for all previous chains, check it



        ###### Residue mapping section ######

        if query_chain_states[query_structchain] == 'holo':
            try:
                pdb_xml = download_sifts_xml_gz(query_struct, pathXML)
            except Exception:
                print('SIFTS server unavailable')  # Handle exception (this should not happen, consider exiting w error)
                #sys.exit(0)

            # Map PDB (binding) residues of query to UniProt residue numbers
            print('Mapping query chain binding residues from PDB to UniProt numbering')
            bndgres_pdb_to_unp =  map_pdb_resnum_to_uniprot(binding_res_unpacked, pdb_xml)
            #print(bndgres_pdb_to_unp)

            # Group (UNP num) binding residues by chain (list in, dict out)
            bndgres_pdb_to_unp_chains = group_mapped_res_by_chain(bndgres_pdb_to_unp)
            print(f'Binding residues [UNP] grouped by query structure chain:\n{query_structchain}:{bndgres_pdb_to_unp_chains.get(query_structchain[4:])}')

            # Find whether mapped binding residues are present in each candidate chain
            candidate_hits = examine_cndt_mapped_bs_res(bndgres_pdb_to_unp_chains, query_structchain, uniprot_segments_dict)# dictApoCandidates_b1)# candidates_unp_dict) 

            # Intermediate step to remove negative score when positive is present for the same position
            candidate_metahits = remove_negative_duplicate_cndt_bs_res_pos(candidate_hits)

            # Count how many binding residues (out of total) are present in candidate
            candidate_scores_chain = evaluate_candidate_bs_rsds(candidate_metahits)

            # Append chain scores to main dict with all scores
            #candidate_scores = candidate_scores | candidate_scores_chain # works in python 3.9+
            candidate_scores.update(candidate_scores_chain) # works for python 3.8

            # Put candidates over/under certain threshold to dicts for further processing (applies on precalculated % scores)
            #dict_rsd_map_candidates[query_structchain] = good_candidates_from_residue_mapping(candidate_scores, bndgrsds_threshold)[query_structchain]
            dict_rsd_map_candidates[query_structchain] = good_candidates_from_residue_mapping(candidate_scores, bndgrsds_threshold).get(query_structchain)
            #bad_candidates_rsd_map[query_structchain] = bad_candidates_from_residue_mapping(candidate_scores, bndgrsds_threshold)[query_structchain]
            bad_candidates_rsd_map[query_structchain] = bad_candidates_from_residue_mapping(candidate_scores, bndgrsds_threshold).get(query_structchain)
            #print(bad_candidates_from_residue_mapping(candidate_scores, bndgrsds_threshold))
            
            # Put bad candidates to discarded chains list
            if bad_candidates_rsd_map.get(query_structchain) is not None:
                for bad_cndt in bad_candidates_rsd_map[query_structchain]:
                    #discarded_chains.append(bad_cndt + '\tBelow mapped binding residue threshold ['+ str(bndgrsds_threshold) + '%]\n')
                    discard_msg = 'Below mapped binding residue threshold ['+ str(bndgrsds_threshold) + '%]'
                    discarded_chains.setdefault(bad_cndt, []).append(discard_msg)

            #print_dict_readable(candidate_hits,'candidate_hits')
            #print_dict_readable(candidate_metahits, 'candidate_metahits')
            #print_dict_readable(candidate_scores_chain, 'candidate_scores_chain') # candidate scores of current query chain
            #print_dict_readable(candidate_scores, 'candidate_scores')
            #print_dict_readable(dictApoCandidates_b1, 'dictApoCandidates_b1')
            #print_dict_readable(dictApoCandidates_1, 'dictApoCandidates_1')
            #print_dict_readable(dict_rsd_map_candidates, 'dict_rsd_map_candidates')
            #print_dict_readable(bad_candidates_rsd_map, 'bad_candidates_rsd_map')

    # end query chain loop


    # Assemble final dictionary of candidates and start candidate alignment loop
    # For holo query chain, get rsd_map candidates, otherwise get UNP overlap candidates
    #print_dict_readable(query_chain_states, 'query_chain_states')
    for key, value in query_chain_states.items():
        if value == 'holo':
            penultimate_candidates[key] = dict_rsd_map_candidates[key]
        elif value == 'apo':
            #final_candidates[key] = dictApoCandidates_1[key]
            penultimate_candidates[key] = sorted(list(dict.fromkeys(dictApoCandidates_1[key]))) # remove duplicates, reorder


    # Remove query chains with no candidates here (#TODO move this upstream)
    for qr_chain in penultimate_candidates:
        if penultimate_candidates.get(qr_chain) is not None:
            final_candidates[qr_chain] = penultimate_candidates[qr_chain]
    
    #print_dict_readable(dictApoCandidates_1,'\ndictApoCandidates_1')
    #print_dict_readable(penultimate_candidates,'\nPenultimate candidates')
    #print_dict_readable(dict_rsd_map_candidates,'\ndict_rsd_map_candidates')
    #print_dict_readable(final_candidates,'\nFinal candidates')
    #sys.exit(0)

    # Move total progress calculation here
    progress_total_candidates = sum([len(lst) for lst in final_candidates.values()])
    progress_processed_candidates = 0

    for query_structchain, candidates_structchains in final_candidates.items():
        print('')
        print(f'=== Processing final candidates for query chain {query_structchain} ===')
        #continue
        if final_candidates[query_structchain] is None: # this has to be skipped upstream to save runtime (also avoid possible error in tests)
            print('chain has no candidates, skipping..')
            continue



        def try_candidate_chain(cmd, query_structchain: str, candidate_structchain: str) -> CandidateChainResult:
            # TODO(chris): collect/move all function side effects to returned CandidateChainResult
            # TODO(rdk): make independent of nonlocal/global variables

            # Global variables that are used or updated in this function
            #apo_holo_dict: dict() # updated
            #apo_holo_dict_H: dict() # updated
            #query_lig_positions: dict() # used
            #cndt_lig_positions: dict() # updated

            apo_holo_dict_instance = dict()
            apo_holo_dict_H_instance = dict()
            cndt_lig_positions_instance = dict()

            nonlocal progress_processed_candidates
            progress_processed_candidates += 1
            track_progress(write_results=True)

            query_struct = query_structchain[:4]
            query_chain = query_structchain[4:]
            candidate_struct = candidate_structchain[:4]
            candidate_chain = candidate_structchain[4:]

            candidate_result = CandidateChainResult(query_struct=query_struct, query_chain=query_chain, candidate_struct=candidate_struct,
                                                candidate_chain=candidate_chain)  # ,query_lig_positions=query_lig_positions

            # Load PyMOL objects (if not loaded)
            query_struct_path = download_mmCIF_gz2(query_struct, pathSTRUCTS)
            if query_struct in cmd.get_object_list('all'):
                pass
            else:
                cmd.load(query_struct_path)
            cmd.select(query_structchain, query_struct + '& chain ' + query_chain)

            candidate_struct_path = download_mmCIF_gz2(candidate_struct, pathSTRUCTS)
            if candidate_struct in cmd.get_object_list('all'):
                pass
            else:
                cmd.load(candidate_struct_path)
            cmd.select(candidate_structchain, candidate_struct + '& chain ' + candidate_chain)

            # Get mapped binding residue scores for query/candidate chain combo
            #print_dict_readable(candidate_scores, 'candidate_scores') # All candidate scores
            #sys.exit(1)
            if query_chain_states[query_structchain] == 'holo':
                bndg_rsd_scores = get_scores_from_residue_mapping(candidate_scores, query_structchain, candidate_structchain)
                if bndg_rsd_scores != '-':
                    bndg_rsd_ratio = bndg_rsd_scores.split()[0]
                    bndg_rsd_percent = bndg_rsd_scores.split()[1][:-1]
                else:
                    bndg_rsd_ratio = bndg_rsd_scores
                    bndg_rsd_percent = bndg_rsd_scores
            elif query_chain_states[query_structchain] == 'apo': # for Apo query chains
                bndg_rsd_ratio = '-'
                bndg_rsd_percent = '-'
            #print(f'\nBinding residue scores for {candidate_structchain}.{query_structchain}: {bndg_rsd_scores}')
            #print(f'Binding residue scores for {candidate_structchain}.{query_structchain}: {bndg_rsd_ratio}, {bndg_rsd_percent}')
            # Pass scores to results
            #candidate_result.bndg_rsd_ratio = bndg_rsd_ratio
            #candidate_result.bndg_rsd_percent = bndg_rsd_percent

            cndtlig_scan_radius = lig_scan_radius

            # Align candidate to query chain with TM-align
            #print(f'\n{candidate_structchain} -> {query_structchain}')

            try:
                aln_tm_scores = tmalign2(cmd, candidate_struct + '& chain ' + candidate_chain, query_struct + '& chain ' + query_chain, object='alnobj', transform=1, quiet=1)
                #aln_rms = cmd.align(candidate_struct + '& chain ' + candidate_chain, query_struct + '& chain ' + query_chain, cutoff=2.0, object='alnobj', cycles=0) # sequence based alignment
                #aln_tm = tmalign2(cmd, candidate_struct + '& chain ' + candidate_chain, query_struct + '& chain ' + query_chain, args='a', transform=1, quiet=1) # return average TM-score
                #rms_cur = cmd.rms_cur(candidate_struct + '& chain ' + candidate_chain, query_struct + '& chain ' + query_chain, cutoff=2.0, object='alnobject')
                #print('rms_cur', rms_cur)

                # Parse alignment scores from TM-align
                if len(aln_tm_scores) == 3:
                    aln_tm = float(aln_tm_scores[1])
                    aln_tm_i = float(aln_tm_scores[0])
                    aln_rms = float(aln_tm_scores[2])  # rms_cur from TM-align
                    aln_tm = round(aln_tm, 3)
                    aln_tm_i = round(aln_tm_i, 3)
                    aln_rms = round(aln_rms, 2)
                elif len(aln_tm_scores) == 2:
                    aln_tm = float(aln_tm_scores[1])
                    aln_tm_i = float(aln_tm_scores[0])
                    aln_rms = '-'
                    aln_tm = round(aln_tm, 3)
                    aln_tm_i = round(aln_tm_i, 3)

                # Print alignment scores from TM-align
                print(f'\n{candidate_structchain} -> {query_structchain}')
                print(f'Alignment scores (RMSD/TM-score/inverse TM-score): [{aln_rms} / {aln_tm} / {aln_tm_i}]')
                #print('aln_tm_scores:', aln_tm_scores)

                # TODO(rdk): which alignment is visualized? And what numbers are reported? - TM-align visualized, TM-scores & RMSD reported from the same TM alignment
                candidate_result.rmsd = aln_rms
                candidate_result.tm_score = aln_tm
                candidate_result.tm_score_i = aln_tm_i

            except Exception as ex:
                print(f'\n{candidate_structchain} -> {query_structchain}')
                print('*Exception: ', ex)
                print('Alignment RMSD/TM score: ERROR')
                print('*poor alignment (error), discarding chain ', candidate_structchain)

                discard_msg = 'Alignment error'
                discarded_chains.setdefault(candidate_structchain, []).append(discard_msg)
                candidate_result.discard_reason = "alignment error"
                return candidate_result


            # Discard poor alignments
            if aln_tm != '-':

                if aln_tm < min_tmscore and aln_tm_i < min_tmscore:
                    print('*poor alignment (below threshold), discarding chain ', candidate_structchain)

                    #discarded_chains.append(candidate_structchain + '\t' + 'Poor alignment (below threshold) [RMSD/TM/iTM]: ' + str(aln_rms) +'/'+ str(aln_tm) +'/'+ str(aln_tm_i) + '\n')
                    discard_msg = 'Poor alignment (below threshold) [RMSD/TM/iTM]: ' + str(aln_rms) +'/'+ str(aln_tm) +'/'+ str(aln_tm_i)
                    discarded_chains.setdefault(candidate_structchain, []).append(discard_msg)
                    candidate_result.discard_reason = "poor alignment (below threshold)"
                    return candidate_result


            # Look for ligands in candidate chain

            # Holo query
            if query_chain_states[query_structchain] == 'holo':
                found_ligands = set()
                found_ligands_xtra = set()
                #print('All query ligands:', query_lig_positions[query_structchain])
                #print(f'Assessing detected ligands in {query_lig_names} and {ligand_names}')
                found_cndt_bs = 0 # keep track of how many binding sites are present in cndt chain out of total
                for ligand in query_lig_positions[query_structchain]:
                    #print('scanning ligand:', ligand)
                    resi = ligand.split()[0]
                    chain = ligand.split()[1]
                    resn = ligand.split()[2]
                    ligand_ = ligand.replace(' ', '_') # remove spaces for selection name


                    # Check whether the candidate covers the area around the superimposed query ligand [Superposition]
                    s1cbs = candidate_structchain + ' and polymer'
                    #s2cbs = 'qBS_CoM_' + ligand_ # pseudoatom center of mass for query ligand binding site
                    s2cbs = '(' + query_struct + ' and chain ' + chain + ' and resi ' + resi + ' and resn ' + resn + ')'
                    cBS_sele = cmd.select('cndtBS_arnd' + ligand_, s1cbs + ' near_to ' + lig_scan_radius + ' of ' + s2cbs)
                    if cBS_sele != 0:
                        found_cndt_bs += 1
                        #print('*candidate residues present around query ligand superposition')
                    #else:
                        #print('*no candidate residues found around query ligand superposition, skipping ligand')
                        #continue # skip this ligand

                    # Look for candidate ligands in the superimposed sites of the aligned query ligands # TODO possible to extract binding site residues here
                    if resn == 'HOH':
                        cndtlig_scan_radius = hoh_scan_radius  # Decide whether to switch back to default after scanning once
                    cndt_sele_expression = '(' + candidate_struct + ' and hetatm)' # only limit search to candidate structure (not chain) #cndt_sele_expression = candidate_struct + ' and chain ' + candidate_chain + ' and hetatm'
                    qr_lig_sele_expression = '(' + query_struct + ' and chain ' + chain + ' and resi ' + resi + ' and resn ' + resn + ')'
                    cmd.select(candidate_structchain + '_arnd_' + ligand_, cndt_sele_expression + ' near_to ' + cndtlig_scan_radius + ' of ' + qr_lig_sele_expression)


                    # If cndt ligand belongs to different PDB chain than candidate structchain, it will not be saved, we need to merge the two selections here
                    cmd.select(candidate_structchain, candidate_structchain + '_arnd_' + ligand_, merge=1)

                    # Put selected atoms in a list, check their name identifiers to see if query ligand name is present
                    myspace_cndt = {'cndt_positions': []}
                    cmd.iterate(candidate_structchain + '_arnd_' + ligand_, 'cndt_positions.append(resi +" "+ chain +" "+ resn)', space = myspace_cndt)

                    # Remove duplicate values from myspace_cndt
                    for key, value in myspace_cndt.items():
                        myspace_cndt[key] = list(myspace_cndt.fromkeys(value))   # preserves the order of values - better than set()

                    #print('scanning ligand:', ligand)
                    print(f'-candidate ligands in query ligand binding site [{ligand}]: {myspace_cndt["cndt_positions"]}')

                    # Find which ligands are present in query by matching their names (without set)
                    for cndt_position in myspace_cndt['cndt_positions']:
                        cndt_lig_resn = cndt_position.split()[2]

                        if cndt_lig_resn in query_lig_names or cndt_lig_resn in ligand_names:
                            found_ligands.add(cndt_lig_resn)
                            cndt_lig_positions_instance.setdefault(candidate_structchain, []).append(cndt_position)
                        elif cndt_lig_resn not in nolig_resn:
                            found_ligands_xtra.add(cndt_lig_resn)
                            if lig_free_sites == 1:
                                cndt_lig_positions_instance.setdefault(candidate_structchain, []).append(cndt_position)

                # end holo ligand loop

                # Count number of binding sites occupied by candidate residues in candidate chain
                found_cndt_bs_ratio = 100 * found_cndt_bs / len(query_lig_positions[query_structchain]) # calculate ratio of found ligands


                # Print verdict for chain & save as ".cif.gz"
                print('\t\t\t\t\t\t\t\t*** Chain evaluation ***')
                print(f'[UNP res. mapping]\tQuery chain binding residues mapped to candidate chain\t[{bndg_rsd_scores}]')
                print(f'[UNP seq. overlap]\tOverall UniProt sequence overlap with query chain\t\t\t[{uniprot_overlap_merged[candidate_structchain][0].split()[1]}%]')
                print(f'[Superposition]\t\tQuery binding sites occupied by candidate residues\t\t[{found_cndt_bs}/{len(query_lig_positions[query_structchain])} {found_cndt_bs_ratio}%]')
                print(f'*specified query ligand(s)/position: {ligand_names}/[{position}]\t verified query ligands: {query_lig_names}\t found query ligands: {found_ligands}\t found non-query ligands: {found_ligands_xtra}')

                #if found_cndt_bs != 0: # Condition to keep/discard candidate # TODO decide whether to keep or not

                ligands_str = join_ligands(found_ligands.union(found_ligands_xtra))
                add_res = ' ' + resolution_dict.get(candidate_struct, '-') + ' ' + r_free_dict.get(candidate_struct, '-') + ' '
                append_expression = candidate_structchain + add_res + uniprot_overlap_merged[candidate_structchain][0].split()[1] + ' ' + bndg_rsd_ratio + ' ' + bndg_rsd_percent + ' ' + str(aln_rms) + ' ' + str(aln_tm) + ' ' + str(aln_tm_i) + ' ' + ligands_str

                # Save apo result
                if lig_free_sites == 1 and len(found_ligands_xtra) == 0 and len(found_ligands) == 0 or lig_free_sites == 0 and len(found_ligands) == 0:
                    apo_holo_dict_instance.setdefault(query_structchain, []).append(append_expression)
                    if len(found_ligands_xtra) > 0:
                        print('APO*')
                    else:
                        print('APO')
                    if save_separate == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results_structs + '/query_' + query_struct + '.cif.gz', query_struct) # save query structure
                        cmd.save(path_results_structs + '/apo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain) # save apo chain

                # Save holo result
                else:
                    apo_holo_dict_H_instance.setdefault(query_structchain, []).append(append_expression)
                    if len(found_ligands) > 0:
                        print('HOLO')
                    else:
                        print('HOLO*')
                    if save_separate == 1 and save_oppst == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results_structs + '/query_' + query_struct + '.cif.gz', query_struct) # save query structure
                        #save_sele = (candidate_structchain + '_arnd_' + ligand_) for ligand in query_lig_positions[query_structchain])
                        cmd.save(path_results_structs + '/holo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain) # save holo chain
                #else:
                    #print(f'Discarding candidate chain [{candidate_structchain}] - no overlapping binding sites found')
                    #discarded_chains.append(candidate_structchain + '\t' + 'no overlapping binding sites' + '\n')
                    #candidate_result.discard_reason = "no overlapping binding sites"
                    #return candidate_result



            # Apo query
            elif query_chain_states[query_structchain] == 'apo':

                found_ligands_r = set()
                found_ligands_xtra = set()
                #found_cndt_bs = 0 # in an apo query, the cndt ligands should be checked to be in areas with query coverage

                # Find ligands in candidate with broad search
                cmd.select('cndt_ligands_' + candidate_structchain, candidate_struct + '& chain ' + candidate_chain + '& hetatm & not (polymer or solvent)')

                # If cndt ligand belongs to different PDB chain than candidate structchain, it will not be found
                # Expand selection to include interface HETATMS (according to superposition)! of both detectable ligands and candidate chain selection
                cndt_sele_expression = '(' + candidate_struct + ' and hetatm and not solvent and not polymer)' # only limit search to candidate structure (not chain) #cndt_sele_expression = candidate_struct + ' and chain ' + candidate_chain + ' and hetatm'
                qr_chain_sele_expression = '(' + query_struct + ' and chain ' + query_chain + ')'
                cmd.select(candidate_structchain + '_HET_arnd_' + query_structchain, cndt_sele_expression + ' near_to ' + cndtlig_scan_radius + ' of ' + qr_chain_sele_expression)

                cmd.select('cndt_ligands_' + candidate_structchain, candidate_structchain + '_HET_arnd_' + query_structchain, merge=1)
                cmd.select(candidate_structchain, candidate_structchain + '_HET_arnd_' + query_structchain, merge=1)

                # Put selected atoms into list
                myspace_r = {'r_positions': []}
                cmd.iterate('cndt_ligands_' + candidate_structchain, 'r_positions.append(resi +" "+ chain +" "+ resn)', space = myspace_r)

                # Remove duplicate values from myspace_r
                for key, value in myspace_r.items():
                    myspace_r[key] = list(myspace_r.fromkeys(value))   # preserves the order of values
                print(f'-candidate ligands found: {myspace_r["r_positions"]}') # use set to remove redundant positions

                # check if binding sites of cndt ligands exist in query
                for ligand in myspace_r['r_positions']:
                    #print('scanning candidate ligand:', ligand)
                    resi = ligand.split()[0]
                    chain = ligand.split()[1]
                    resn = ligand.split()[2]
                    ligand_ = ligand.replace(' ', '_') # remove spaces for selection name

                    # Make selection of cndt binding site
                    s1qchain = query_structchain + ' and polymer'
                    s2cndtlig = '(' + candidate_struct + ' and chain ' + chain + ' and resi ' + resi + ' and resn ' + resn + ')'
                    cndt_bs = cmd.select('clig_on_qbs_' + ligand_, s1qchain + ' near_to ' + lig_scan_radius + ' of ' + s2cndtlig)

                    if cndt_bs == 0:
                        print('*ligand validity: [invalid]\nno query chain residues around candidate ligand superposition, ignoring ligand')
                        continue # skip this ligand
                    else:
                        #print('*ligand validity: [valid]')  # query residues present around candidate ligand superposition

                        # Register valid ligands
                        r_atom_lig_name = ligand.split()[2]

                        if r_atom_lig_name not in nolig_resn:  # exclude non ligands
                            found_ligands_r.add(r_atom_lig_name)
                            cndt_lig_positions_instance.setdefault(candidate_structchain, []).append(ligand)


                # Print verdict for chain & save it as ".cif.gz"
                print('\t\t\t\t\t\t\t\t*** Chain evaluation ***')
                print(f'[UNP seq. overlap] Percentage of overall UniProt sequence overlap with query chain: [{uniprot_overlap_merged[candidate_structchain][0].split()[1]}]')
                print(f'*found ligands: {found_ligands_r}')

                add_res = ' ' + resolution_dict.get(candidate_struct, '-') + ' ' + r_free_dict.get(candidate_struct, '-') + ' '
                # Save holo result
                if len(found_ligands_r) > 0:
                    ligands_str = join_ligands(found_ligands_r)
                    append_expression = candidate_structchain + add_res + uniprot_overlap_merged[candidate_structchain][0].split()[1] + ' ' + bndg_rsd_ratio + ' ' + bndg_rsd_percent + ' ' + str(aln_rms) + ' ' + str(aln_tm) + ' ' + str(aln_tm_i) + ' ' + ligands_str
                    apo_holo_dict_H_instance.setdefault(query_structchain, []).append(append_expression)
                    print('HOLO')
                    if save_separate == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results_structs + '/query_' + query_struct + '.cif.gz', query_struct) # save query structure
                        cmd.save(path_results_structs + '/holo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain) # save apo chain

                # Save apo result
                else:
                    ligands_str = join_ligands(found_ligands_r.union(found_ligands_xtra))
                    append_expression = candidate_structchain + add_res + uniprot_overlap_merged[candidate_structchain][0].split()[1] + ' ' + bndg_rsd_ratio + ' ' + bndg_rsd_percent + ' ' + str(aln_rms) + ' ' + str(aln_tm) + ' ' + str(aln_tm_i) + ' ' + ligands_str
                    apo_holo_dict_instance.setdefault(query_structchain, []).append(append_expression)
                    if len(found_ligands_xtra) > 0:
                        print('APO*')
                    else:
                        print('APO')
                    if save_separate == 1 and save_oppst == 1:
                        if not os.path.isfile(path_results + '/query_' + query_struct + '.cif.gz'):
                            cmd.save(path_results_structs + '/query_' + query_struct + '.cif.gz', query_struct)  # save query structure
                        cmd.save(path_results_structs + '/apo_' + candidate_structchain + '_aligned_to_' + query_structchain + '.cif.gz', candidate_structchain)  # save holo chain

            # end apo/holo candidate verdict


            candidate_result.apo_holo_dict_instance = apo_holo_dict_instance
            candidate_result.apo_holo_dict_H_instance = apo_holo_dict_H_instance
            candidate_result.cndt_lig_positions_instance = cndt_lig_positions_instance

            candidate_result.passed = True
            return candidate_result


        def try_candidate_chains(cmd, query_structchain, candidates_structchains: list) -> list: # list[CandidateChainResult]:
            query_parallelism = args.query_parallelism

            results = []
            if query_parallelism == 1:
                # Start candidate chain loop
                # Align candidate chain to query chain and mark atom selections around superimposed query ligand binding sites
                for cand in candidates_structchains:
                    # TODO(rdk): make cmd independent for each run
                    results.append(try_candidate_chain(cmd, query_structchain, cand))
            else:
                def _try_chain(cand_strchain):
                    # TODO(rdk): make cmd independent for each run

                    local_pm = pymol2.PyMOL()
                    local_pm.start()
                    local_cmd = local_pm.cmd
                    local_cmd.load(query_struct_path)

                    return try_candidate_chain(local_cmd, query_structchain, cand_strchain)

                with ThreadPoolExecutor(max_workers=query_parallelism) as pool:
                # with ProcessPoolExecutor(max_workers=query_parallelism) as pool:
                    results = list(pool.map(_try_chain, candidates_structchains))

            return results


        def update_dict_of_lists(modified_dict: dict, to_add: dict):
            for key, value in to_add.items():
                modified_dict.setdefault(key, []).extend(value)


        candidate_results = try_candidate_chains(cmd, query_structchain, candidates_structchains)
        #TODO integrate results to global state here, write to disk
        passed_results = [cr for cr in candidate_results if cr.passed]

        # Merge/update dictionaries
        for res in passed_results:
            update_dict_of_lists(apo_holo_dict, res.apo_holo_dict_instance)            # key = query_structchain
            update_dict_of_lists(apo_holo_dict_H, res.apo_holo_dict_H_instance)        # key = query_structchain
            update_dict_of_lists(cndt_lig_positions, res.cndt_lig_positions_instance)  # key = candidate_structchain

        # Remove duplicate values from dictionaries (cndt_lig_positions dict may have duplicates)
        for key, value in cndt_lig_positions.items():
            cndt_lig_positions[key] = list(cndt_lig_positions.fromkeys(value))  


    track_progress()
    # end for

    write_pymol_script(path_results)

    print('')

    # Print chains that were discarded
    if len(discarded_chains) > 0:
        print('Discarded chains: ', len(discarded_chains))
        #print(f"{' '.join(map(str, discarded_chains))}\n")


    # Universal results
    
    # Remove duplicate discard msgs from dict "discarded_chains" 
    for key, value in discarded_chains.items():
        discarded_chains[key] = list(discarded_chains.fromkeys(value))
    write_ligands_csv(query_lig_positions, cndt_lig_positions, path_results)
    write_discarded_chains(discarded_chains, path_results)

    # Apo results
    write_results_apo_csv(apo_holo_dict, path_results)
    num_apo_chains = sum([len(apo_holo_dict[x]) for x in apo_holo_dict if isinstance(apo_holo_dict[x], list)])  # number of found APO chains
    print('\nApo chains: ', num_apo_chains)
    for key in apo_holo_dict:
        print(key, apo_holo_dict.get(key))

    # Holo results
    write_results_holo_csv(apo_holo_dict_H, path_results)
    num_holo_chains = sum([len(apo_holo_dict_H[x]) for x in apo_holo_dict_H if isinstance(apo_holo_dict_H[x], list)])  # number of found HOLO chains
    print('\nHolo chains: ', num_holo_chains)
    for key in apo_holo_dict_H:
        print(key, apo_holo_dict_H.get(key))


    if len(apo_holo_dict) == 0 and len(apo_holo_dict_H) == 0:
        print('\nConsider reversing the search or revising the input query')
        # Note: we don't want to delete empty job folder but keep it for potential further processing by the webserver
        # Delete empty job folder
        # print('\nDeleting empty job folder')    
        # try:
        #     os.rmdir(pathRSLTS)
        # except OSError as error:
        #     print('Job folder not empty, job ID: ', job_id, error)

    # Print states of query chains (apo or holo)
    print(f'\nQuery chain states:\n{query_chain_states}')

    # Test print query & candidate ligand positions
    print(f'Query ligands\n{query_lig_positions}\n')
    #print(f'\nCandidate ligands\n{cndt_lig_positions}')

    # Append the name of the query and the job_id in the queries.txt
    if job_id:
        print(f'Saving query search: {[user_query_parameters]} for query {[query]}') #query_full)  # Don't show full query string
        with open(pathQRS + '/' + 'queries.txt', 'a') as out_q:
            out_q.write(query_full + '\n')

    # Print argparse arguments
    #print('\n', args)
    #print(vars(args))
    # Print no lig list
    #print('\nNo lig list:', nolig_resn)

    pm.stop()

    print('Results saved to directory:\t', path_results)
    print('Done processing query: ', query)

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

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
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
    #parser.add_argument('--query', type=str,   default='1a73')
    #parser.add_argument('--query', type=str,   default='1a73 A zn',    help='main input query') # OK apo 0, holo 16
    #parser.add_argument('--query', type=str,   default='1a73 A,B zn',  help='main input query') # OK apo 0, holo 32
    #parser.add_argument('--query', type=str,   default='1a73 * zn',    help='main input query') # OK apo 0, holo 32
    #parser.add_argument('--query', type=str,   default='1a73 E mg 205',help='main input query') # fail, ligand assigned non-polymer chain
    #parser.add_argument('--query', type=str,   default='1a73 A',       help='main input query') # apo 0, holo 16
    #parser.add_argument('--query', type=str,   default='1a73 * *',     help='main input query') # OK apo 0, holo 32
    #parser.add_argument('--query', type=str,   default='1a73 * mg',    help='main input query') # one MG is on non-protein chain
    #parser.add_argument('--query', type=str,   default='1a73 ! mg',    help='main input query') # apo 8, holo 24. one MG is on non-protein chain
    #parser.add_argument('--query', type=str,   default='1a73 ! mg,zn',    help='main input query') # apo
    #parser.add_argument('--query', type=str,   default='3vro ! ptr',    help='main input query') # apo 7, holo 0. 1 PTR between 2 chains (assigned chain is tiny fragment)
    #parser.add_argument('--query', type=str,   default='5aep ! ptr',    help='main input query') # 2 PTRs in same chain non-interface
    #parser.add_argument('--query', type=str,   default='5aep ! hem',    help='main input query') 
    parser.add_argument('--query', type=str,   default='3n7y ! ptr',    help='main input query') # apo 21, holo 270. good test for "!"
    #parser.add_argument('--query', type=str,   default='5j72 A na 703',help='main input query') # apo 0, holo 0 (no UniProt chains)
    #parser.add_argument('--query', type=str,   default='1a73 b mg 206',help='main input query') # OK, apo 4, holo 12
    #parser.add_argument('--query', type=str,   default='1a73 b mg 206',help='main input query') # water_as_ligand=1 OK, apo 4, holo 12
    #parser.add_argument('--query', type=str,   default='7s4z A *',     help='main input query') # apo 103, holo 104 *many irrelevant ligands
    #parser.add_argument('--query', type=str,   default='3fav all zn',  help='main input query')
    #parser.add_argument('--query', type=str,   default='3fav all',     help='main input query')
    #parser.add_argument('--query', type=str,   default='1y57 A mpz',   help='main input query') # apo 5, holo 29
    #parser.add_argument('--query', type=str,   default='6h3c B,G zn',  help='main input query') # OK apo 0, holo 4
    #parser.add_argument('--query', type=str,   default='2v0v',         help='main input query') # Fully apo, apo 8, holo 24
    #parser.add_argument('--query', type=str,   default='2v0v A,B',     help='main input query') # apo 4, holo 12
    #parser.add_argument('--query', type=str,   default='2v0v A',       help='main input query') # apo 2, holo 6
    #parser.add_argument('--query', type=str,   default='2hka all c3s', help='main input query') # OK apo 2, holo 0
    #parser.add_argument('--query', type=str,   default='2v57 A,C prl', help='main input query') # OK apo 4, holo 0
    #parser.add_argument('--query', type=str,   default='3CQV all hem', help='main input query') # OK apo 6, holo 5
    #parser.add_argument('--query', type=str,   default='3CQV A hem', help='main input query') # OK apo 6, holo 5
    #parser.add_argument('--query', type=str,   default='3CQV ! hem', help='main input query') #
    #parser.add_argument('--query', type=str,   default='2npq A bog',   help='main input query') # long, apo 149, holo 114, p38 MAP kinase cryptic sites
    #parser.add_argument('--query', type=str,   default='1ksw A NBS',   help='main input query') # apo 4, holo 28 Human c-Src Tyrosine Kinase (Thr338Gly Mutant) in Complex with N6-benzyl ADP
    #parser.add_argument('--query', type=str,   default='1ai5',         help='main input query') # negative uniprot overlap (fixed)
    #parser.add_argument('--query', type=str,   default='2hka all c3s', help='main input query') # first chain is apo, it was ignored before now works
    #parser.add_argument('--query', type=str,   default='2hka ALL',     help='main input query') #
    #parser.add_argument('--query', type=str,   default='6j19 all atp', help='main input query') # problematic case, 6j19B has wrong UNP mapping in SIFTS 
    #parser.add_argument('--query', type=str,   default='1aro P HG 904',   help='main input query') # fragmented UniProt candidates, to use for testing UNP overlap calculation
    #parser.add_argument('--query', type=str,   default='4V51 BA MG 3327', help='main input query') # ribosome protein binding to nucleic acid only
    #parser.add_argument('--query', type=str,   default='4V51 BA MG 3328', help='main input query') # ribosome protein binding also protein (ok)
    #parser.add_argument('--query', type=str,   default='4V51 ! MG', help='main input query') # too slow for "!"
    #parser.add_argument('--query', type=str,   default='1GB1',         help='main input query') # apo 20, holo 16, apo struct, has solid state nmr candidate "2K0P"
    #parser.add_argument('--query', type=str,   default='6hwv A BOG 402',  help='main input query') # apo 134, holo 128
    #parser.add_argument('--query', type=str,   default='6hwv A BOG',   help='main input query') # apo 76, holo 186. Bug with residue mapping section (maps only first ligand, then transfers the binding residues to rest ligands)

    # Issue: Ligands bound to query protein chain (interaface) but annotated to different chain (either of the protein or the polymer/nucleic acid)
    #parser.add_argument('--query', type=str,   default='6XBY A adp,mg',  help='main input query') # apo 4, holo 2
    #parser.add_argument('--query', type=str,   default='6XBY * adp,mg',  help='main input query') # ATPase, big query
    #parser.add_argument('--query', type=str,   default='6XBY s nag')
    #parser.add_argument('--query', type=str,   default='6XBY b pov')
    #parser.add_argument('--query', type=str,   default='6XBY A thr 257', help='main input query')
    #parser.add_argument('--query', type=str,   default='1a73 e mg 205')
    #parser.add_argument('--query', type=str,   default='1a73 A mg,zn')
    #parser.add_argument('--query', type=str,   default='1a73 * mg,zn')
    #parser.add_argument('--query', type=str,   default='1a73 * mg')

    # Issue part B: ligand specified to a correct but non-protein chain
    #parser.add_argument('--query', type=str,   default='1a73 E mg 205', help='main input query') # apo 4, holo 12 (fixed) fail, ligand assigned non-polymer chain

    # Residue
    #parser.add_argument('--query', type=str,   default='1a73 A ser',    help='main input query') # expected parsing fail
    #parser.add_argument('--query', type=str,   default='1a73 A ser 97', help='main input query') # OK apo 4, holo 12

    # SARS CoV 2
    #parser.add_argument('--query', type=str,   default='7aeh A R8H',     help='main input query') # apo 116, holo 431, job 314 H PC, CoV2 Mpro caspase-1 inhibitor SDZ 224015
    #parser.add_argument('--query', type=str,   default='2amq A his 164', help='main input query') # apo 41 (43*), holo 66, job 318 H, Crystal Structure Of SARS_CoV Mpro in Complex with an Inhibitor N3
    #parser.add_argument('--query', type=str,   default='6lu7 A R8H',     help='main input query')
    #parser.add_argument('--query', type=str,   default='7krn A adp',     help='main input query')
    #parser.add_argument('--query', type=str,   default='7krn ! adp',     help='main input query')

    # Water
    #parser.add_argument('--query', type=str,   default='1a73 B hoh 509',  help='main input query') # OK apo 6, holo 10 (previous apo 9, holo 7)
    #parser.add_argument('--query', type=str,   default='1a73 * hoh 509',  help='main input query') # OK apo 9, holo 7
    #parser.add_argument('--query', type=str,   default='1a73 * hoh',      help='main input query') # expected parsing fail
    #parser.add_argument('--query', type=str,   default='3i34 X hoh 311',  help='main input query') # apo 113, holo 94 *many irrelevant ligands show up
    #parser.add_argument('--query', type=str,   default='1pkz A tyr 9',    help='main input query') # apo 7, holo 93, water as lig, marian, allosteric effect of hoh
    #parser.add_argument('--query', type=str,   default='1fmk A HOH 1011', help='main input query') # Issue related (query longer than candidate seq, poor one-way TM score, hit 4hxj is discarded)
    #parser.add_argument('--query', type=str,   default='4hxj A,B',        help='main input query')
    #parser.add_argument('--query', type=str,   default='2v7d B hoh 2026', help='main input query')

    # Non standard residues
    #parser.add_argument('--query', type=str,   default='6sut A tpo',     help='main input query') # OK apo 0, holo 3
    #parser.add_argument('--query', type=str,   default='6sut ! tpo',     help='main input query') # 
    #parser.add_argument('--query', type=str,   default='6sut A',         help='main input query') # OK apo 0, holo 3
    #parser.add_argument('--query', type=str,   default='6sut A tpo 285', help='main input query') # OK apo 0, holo 3
    #parser.add_argument('--query', type=str,   default='1a73 A zn 201',  help='main input query') # OK apo 0, holo 16
    #parser.add_argument('--query', type=str,   default='1a37 P sep 259', help='main input query') # OK apo 6, holo 4 [with nonstd_rsds_as_lig=1]
    #parser.add_argument('--query', type=str,   default='1a37 ! sep', help='main input query') #
    #parser.add_argument('--query', type=str,   default='1a37 + sep', help='main input query') # Wrong input format
    #parser.add_argument('--query', type=str,   default='1a37 P sep',     help='main input query') # OK apo 6, holo 4 [with nonstd_rsds_as_lig=1], without now also works
    #parser.add_argument('--query', type=str,   default='1apm E sep',     help='main input query') # OK apo 5, holo 76
    #parser.add_argument('--query', type=str,   default='1apm E sep,tpo', help='main input query') # OK apo 2, holo 79
    #parser.add_argument('--query', type=str,   default='1apm E sep,ser', help='main input query') # wrong query, should get error

    # Long queries
    #parser.add_argument('--query', type=str,   default='1cim A', help='main input query')
    #parser.add_argument('--query', type=str,   default='1h1p A', help='main input query')
    #parser.add_argument('--query', type=str,   default='1k1j A', help='main input query')

    # D amino acids
    #parser.add_argument('--query', type=str,   default='148l S DAL 170', help='main input query') # not working, not registered as a ligand, chain S is non-UniProt, no other UniProt chains around
    #parser.add_argument('--query', type=str,   default='148l E ARG 137', help='main input query') # (750 structures) does not seem to pick up DAL as ligand even with setting turned on
    #parser.add_argument('--query', type=str,   default='1cfa', help='main input query') # works and gives results
    #parser.add_argument('--query', type=str,   default='1cfa B dar', help='main input query') # works after handling exception @ line 922 (many negative unp overlaps?)

    # Non-UniProt query structure (7MJB)
    #parser.add_argument('--query', type=str,   default='7mjb',       help='main input query')
    #parser.add_argument('--query', type=str,   default='5IBO * DKA', help='main input query')  # This is in UniProt, same sequence to previous one (we don't find it)
    #parser.add_argument('--query', type=str,   default='5IBO ! DKA', help='main input query')  # 

    # Test new UNP overlap computation
    #parser.add_argument('--query', type=str,   default='7khr B')  # apo7, holo 2 (previous error? apo 30, holo 4)
    #parser.add_argument('--query', type=str,   default='3fav all zn',  help='main input query')

    # Basic
    parser.add_argument('--res_threshold',     type=float, default=3.8,   help='Lowest allowed resolution for result structures (applies to highest resolution value for scattering methods, expressed in angstroms), condition is <=')
    parser.add_argument('--include_nmr',       type=int,   default=1,     help='0/1: Discard/include NMR structures')
    parser.add_argument('--xray_only',         type=int,   default=0,     help='0/1: Only consider X-ray structures')
    parser.add_argument('--lig_free_sites',    type=int,   default=1,     help='0/1: Ligand-free binding sites. When on, resulting apo sites will be free of any other known ligands in addition to specified ligands')
    #parser.add_argument('--autodetect_lig',    type=int,   default=0,     help='0/1: This will find and consider any non-protein and non-solvent heteroatoms as ligands and mark their binding sites, in addition to any specified ligands (useful when the user does not know the ligand)')
    #parser.add_argument('--reverse_search',    type=int,   default=0,     help='0/1: Start the search with an apo structure that does not bind any ligands')
    parser.add_argument('--water_as_ligand',   type=int,   default=0,     help='0/1: When examining the superimposed binding sites of candidate structures, consider HOH molecules as ligands and show them in the results')

    # Advanced
    parser.add_argument('--overlap_threshold', type=float, default=0.0,   help='Minimum % of sequence overlap between query and result chains (using the SIFTS residue-level mapping with UniProt), condition is ">="')
    parser.add_argument('--bndgrsds_threshold',type=float, default=1.0,   help='Percentage of binding residues of the query that have to be present in the candidate according to UNP residue mapping, for the candidate to be considered, condition is ">="')
    parser.add_argument('--lig_scan_radius',   type=float, default=4.5,   help='Angstrom radius to look around the query ligand(s) superposition (needs to be converted to str)')

    parser.add_argument('--min_tmscore',       type=float, default=0.5,   help='Minimum acceptable TM score for apo-holo alignments (condition is ">")')
    parser.add_argument('--nonstd_rsds_as_lig',type=int,   default=0,     help='0/1: Ignore/consider non-standard residues as ligands')
    parser.add_argument('--d_aa_as_lig',       type=int,   default=0,     help='0/1: Ignore/consider D-amino acids as ligands')

    # Experimental
    #parser.add_argument('--beyond_hetatm',     type=int,   default=0,     help='0/1: when enabled, does not limit holo ligand detection to HETATM records for specified ligand/residue')  # [might need to apply this to apo search too
    parser.add_argument('--look_in_archive',   type=int,   default=0,     help='0/1: Search if the same query has been processed in the past (can give very fast results)')

    # Internal
    parser.add_argument('--apo_chain_limit',   type=int,   default=9999,  help='limit number of apo chains to consider when aligning (for fast test runs)')
    parser.add_argument('--work_dir',          type=str,   default=None,  help='global root working directory for pre-computed and intermediary data')
    parser.add_argument('--out_dir',           type=str,   default=None,  help='explicitly specified output directory')
    parser.add_argument('--threads',           type=int,   default=4,     help='number of concurrent threads for processing multiple queries')
    parser.add_argument('--query_parallelism', type=int,   default=2,     help='number of concurrent threads for processing single query')
    parser.add_argument('--track_progress',    type=bool,  default=False, help='track the progress of long queries in .progress file, update result csv files continually (not just at the end)')
    parser.add_argument('--intrfc_lig_radius', type=float, default=4.5,   help='Angstrom radius to look around atoms of ligand for interactions with protein atoms')
    parser.add_argument('--hoh_scan_radius',   type=float, default=2.5,   help='Angstrom radius to look around the query ligand(s) superposition (needs to be converted to str, applies to water query ligands only)')

    # Saving
    parser.add_argument('--save_oppst',        type=int,   default=1,     help='0/1: also save chains same with query (holo chains when looking for apo, and apo chains when looking for holo)')
    parser.add_argument('--save_separate',     type=int,   default=1,     help='0/1: save each chain object in a separate file (default save)')
    #parser.add_argument('--save_session',      type=int,   default=0,     help='0/1: save each result as a PyMOL ".pse" session (zipped, includes annotations -less recommended)')
    #parser.add_argument('--multisave',         type=int,   default=0,     help='0/1: save each result in a .pdb file (unzipped, no annotations -least recommended)')


    return parser.parse_args(argv)


def main(argv):
    print(f'APO-HOLO JUXTAPOSITION (v{VERSION})')
    time_local = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    print(time_local)  # print(f'Local time [{time_local}]')
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
