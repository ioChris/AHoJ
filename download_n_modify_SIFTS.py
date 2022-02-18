# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 18:27:43 2022

@author: ChrisX
"""
'''
Download new SIFTS file with PDB chain mapping and unzip it in specified folder
Create 2 dictionaries: 
    i) a shorter version of the SIFTS file with joined structure and chain
    ii) a reverse version with structchains as values of UniProt ID(s)
    
    -backups and deletes old SIFTS and dict files

    
'''
import os
import gzip
import wget
import shutil
import time
import argparse

# Input arguments

parser = argparse.ArgumentParser()

parser.add_argument('--work_directory',  type=str,   default=None,   help='root working directory for pre-computed and intermediary data')
parser.add_argument('--pdb_uniprot_url', type=str,   default='http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz')

args = parser.parse_args()


# Saving options (default = 1)
del_old_sifts = 0 # 0 means processing existing SIFTS file and not downloading the new one
save_dicts = 0 # if 0, overrides following settings

save_r_spnum = 1    # save dict_R_SPnum dictionary to file (includes readable version w headers)
save_plaindict = 1  # save dict_SIFTS dictionary to file (includes readable version w headers)

if save_dicts == 0:
    save_r_spnum = 0
    save_plaindict = 0

use_old_sifts = False


##########################################################################################################
# Define functions
##########################################################################################################

def get_2level_cwd():
    npath = os.path.normpath(os.getcwd())   # Normalize the path string for the OS
    path0 = os.path.join(npath.split(os.sep)[0], '/', npath.split(os.sep)[1], npath.split(os.sep)[2])
    if os.path.exists(path0):
        memo = "Detected current working >> " + path0
    else:
        memo = 'Error detecting current working dir:' + npath
    print(memo)
    return path0


def work_directory(args):
    if (args.work_directory):
        return args.work_directory
    else:
        return get_2level_cwd() + r'\Documents\Bioinfo_local\Ions\datasets_local\APO_candidates\webserver'  # default work directory

##########################################################################################################


path_root = work_directory(args)
pathSIFTS = path_root + r'/SIFTS'           # Pre compiled files with UniProt PDB mapping
pathOLD = path_root + r'/oldSIFTS'          # old SIFTS files, keep as backup
# pathSIFTS = get_2level_cwd() + r'\ownCloud\Bioinfo_ownCloud\Projects\Ions\Uniprot_PDBchain\autodownload'
# pathSIFTS = r'C:\Users\TopOffice\Documents\GitHub\workDir\files'
pdb_uniprot_url = args.pdb_uniprot_url
filename = pdb_uniprot_url.split('/')[-1]

print('Downloading SIFTS files from EBI ftp to local folder')
print('Download URL:', pdb_uniprot_url, '\n') # Print the URL that is being checked/downloaded
print('Target filename: "%s"' % filename) # Print name of file to download
print('Saving directory:', pathSIFTS)

# Check if SIFTS directory exists, if not, create it
if not os.path.exists(pathSIFTS):
    print('Working directory does not exist, creating new directory\n')
    try:
        os.makedirs(pathSIFTS)
    except:
        print('Could not create directory: ', pathSIFTS)

# Check if (old) SIFTS file exists and/or download new file    
if os.path.exists(pathSIFTS + "/" + filename):
    if del_old_sifts == 1:

        # Archive old SIFTS files before deleting/replacing
        print('Backing up current SIFTS and dict files')
        if not os.path.exists(pathOLD):            os.makedirs(pathOLD)
        time_str = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        output_archive_pathname = pathOLD + '/oldSIFTS_' + time_str
        shutil.make_archive(output_archive_pathname, 'zip', pathSIFTS)

        print('Deleting existing SIFTS file')
        os.remove(pathSIFTS + "/" + filename)

        print('Downloading new SIFTS file into ', pathSIFTS)
        wget.download(pdb_uniprot_url, pathSIFTS)  # Download new file
    else:
        print('Using existing SIFTS file')
        use_old_sifts = True
else:
    print('Downloading new SIFTS file into ', pathSIFTS)
    wget.download(pdb_uniprot_url, pathSIFTS)  # Download new file



# Unzip new .gz file (and overwrite old)
if not use_old_sifts:
    with gzip.open(pathSIFTS + "/" + filename ,"rb") as infile, open(pathSIFTS + "/" + filename[:-3] ,"wb") as outfile:
        print('Unzipping new file')
        for line in infile:
            outfile.write(line)
    print('Done\n')


## Create reverse version with SPnum and also a simple version of the dict with structchain:uniprot_id
fileSIFTS = pathSIFTS + "/" + filename[:-3]
filenameSIFTS = filename[:-3]
counter_zeroes = 0
counter_neg = 0
dict_R_SPnum = dict()   # The dictionary with the reversed SIFTS chains/structures and SPnum
dict_SIFTS = dict()     # Normal forward dict, structure -> structchain:uniprot_id


with open(fileSIFTS, 'r') as infile1:

    # Get/skip first two header lines
    header1 = infile1.readline()
    header2 = infile1.readline()
    print("\nBuilding dictionary")

    for line in infile1:

        structure = line.split("\t")[0]
        chain = line.split("\t")[1]
        uniprot_id = line.split("\t")[2]    #uniprot_id = line.split()[2]
        SIFTSstructchain = structure + chain

        SP_BEG = int(line.split("\t")[7])
        SP_END = int(line.split("\t")[8][:-1])    # includes \n cause of EOL
        sp_length = SP_END - SP_BEG     # Calculate length of chain

        # Build basic dict with structchain as keys, uniprot IDs as values
        dict_SIFTS[SIFTSstructchain] = uniprot_id

        # Append values to the R_SPnum dict (SP_BEG, SP_END, sp_length)
        dict_R_SPnum.setdefault(uniprot_id, []).append(structure+chain+' '+str(SP_BEG)+' '+str(SP_END))
        #dict_R_SPnum.setdefault(uniprot_id, []).append(structure + chain)  # basic version with UniProt ID as the key and corresponding PDB structures+chains as values

        # Count the cases where the length/coverage is zero
        if sp_length == 0:
            counter_zeroes += 1
        elif sp_length < 0:
            counter_neg += 1
            #print(uniprot_id, structure+chain+' '+str(SP_BEG)+' '+str(SP_END)+' '+str(sp_length))

    print('Done')


print('Dict keys, values:')
print(len(dict_R_SPnum.keys()), sum(len(v) for v in dict_R_SPnum.values()))
print('Entries with 0 length:\t', counter_zeroes)
print('Entries with negative length:\t', counter_neg)

'''
# Discard duplicate values in SIFTS_dict (values)
print('Removing duplicate entries')
for key,value in SIFTS_dict.items():
    SIFTS_dict[key] = list(SIFTS_dict.fromkeys(value))   # preserves the order of values
print('Done\n')
'''

## Save dicts to file (include readable versions)

# Reverse dict (dict_R_SPnum)
if save_r_spnum == 1:
    outfile1 = filenameSIFTS[:-4] + "_REVERSE_SPnum.txt"
    outfile2 = outfile1[:-4] + "_readable.txt"
    print('Saving output of dict_R_SPnum as:\n', outfile1, '\t', outfile2)

    with open (pathSIFTS + '/' + outfile1, 'wt') as out1:
        out1.write(str(dict_R_SPnum))

    # Write same dict in a more readable format (new line per key)
    with open (pathSIFTS + '/' + outfile2, 'wt') as out2:
        for key, value, in dict_R_SPnum.items():
            out2.write('{0}: {1}\n'.format(key, value))

# Forward simple dict (dict_SIFTS)
if save_plaindict == 1:
    outfile_SIFTS = filenameSIFTS[:-4] + '_dict.txt'
    outfile_SIFTS_readable = outfile_SIFTS[:-4] + "_readable.txt"
    print('Saving output of dict_SIFTS as:\n', outfile_SIFTS, '\t', outfile_SIFTS_readable)

    with open (pathSIFTS + '/' + outfile_SIFTS, 'wt') as out3:
        out3.write(str(dict_SIFTS))

    # Write  dict2 in a more readable format (new line per key)
    with open (pathSIFTS + '/' + outfile_SIFTS_readable, 'wt') as out4:
        #header = '#HEADER: HOLO_chain SP_BEG SP_END : APO_chain SP_BEG SP_END SP_LEN overlap_len %overlap\n'
        out4.write(header1)
        out4.write(header2)
        for key, value, in dict_SIFTS.items():
            out4.write('{0}: {1}\n'.format(key, value))

print('All done')
