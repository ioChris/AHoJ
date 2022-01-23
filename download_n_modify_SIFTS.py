# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 18:27:43 2022

@author: TopOffice
"""
'''
Download new SIFTS file with PDB chain mapping and unzip it in specified folder
Creates 2 dictionaries: 
    i) one short version of the SIFTS file with joint structure and chain
    ii) one reverse version with structchains as values of UniProt ID(s)
    
    -deletes old SIFTS file
    -deletes old dict files
    
'''

# Saving options
save_r_spnum = 0    # save dict_R_SPnum dictionary to file (includes readable version w headers)
save_plaindict = 0  # save dict_SIFTS dictionary to file (includes readable version w headers)

import os
import gzip
import wget
#from get_root_path import root_path

url1 = "http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
   
#path1 = root_path()  + '\\' + r'ownCloud\Bioinfo_ownCloud\Projects\Ions\Uniprot_PDBchain\autodownload'
path1 = r'C:\Users\TopOffice\Documents\GitHub\workDir\files'
filename = url1.split('/')[-1]

print('Downloading SIFTS files from EBI ftp to local folder') # Print the description of the script as a message
print("Download URL:", url1, '\n') # Print the URL that is being checked/downloaded
print('Checking file: "%s"' % filename) # Print name of file to download
print('Saving directory:', path1)

# Clean old .gz file
if os.path.exists(path1 + "\\" + filename):
    print('Deleting existing file')
    os.remove(path1 + "\\" + filename)
        
# Download new file
print('Downloading SIFTS file into ', path1)
wget.download(url1, path1)

# Unzip new .gz file (and overwrite old)
with gzip.open(path1 + "\\" + filename ,"rb") as infile, open(path1 + "\\" + filename[:-3] ,"wb") as outfile:
    print('Unzipping new file')
    for line in infile:
        outfile.write(line)
print('Done\n')


## Create reverse version with SPnum and also a simple version of the dict with structchain:uniprot_id
fileSIFTS = path1 + "\\" + filename[:-3]
filenameSIFTS = filename[:-3]
counter_zeroes = 0  #
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
        uniprot_id = line.split("\t")[2]
        #uniprot_id = line.split()[2]
        SIFTSstructchain = line.split('\t')[0] + line.split('\t')[1]
        
        SP_BEG = int(line.split("\t")[7])
        SP_END = int(line.split("\t")[8][:-1])    # includes \n cause of EOL     
        sp_length = SP_END - SP_BEG     # Calculate length of chain
        
        # Build basic dictionary with UniProt ID as the key and corresponding PDB structures+chains as values
        #dict_R_SPnum.setdefault(uniprot_id, []).append(structure + chain)
        
        # Build basic dict with structchain as keys, uniprot IDs as values
        dict_SIFTS[SIFTSstructchain] = uniprot_id
        
        # Append values to the R_SPnum dict (SP_BEG, SP_END, sp_length)
        dict_R_SPnum.setdefault(uniprot_id, []).append(structure+chain+' '+str(SP_BEG)+' '+str(SP_END))
        
        # Count the cases where the length/coverage is zero
        if sp_length == 0:
            counter_zeroes += 1
            #print(uniprot_id, structure+chain+' '+str(SP_BEG)+' '+str(SP_END)+' '+str(sp_length))
    
    print('Done')
    

print('Dict keys, values:')
print(len(dict_R_SPnum.keys()), sum(len(v) for v in dict_R_SPnum.values()))
print('Entries with 0 length:\t', counter_zeroes)

## Save dicts to file (include readable versions)

# R_SPnum dict
if save_r_spnum == 1:
    outfile1 = filenameSIFTS[:-4] + "_REVERSE_SPnum.txt"
    outfile2 = outfile1[:-4] + "_readable.txt"
    print('Saving output of dict_R_SPnum as ', outfile1, outfile2)
    
    with open (path1 + '\\' + outfile1, 'wt') as out1:
        out1.write(str(dict_R_SPnum))
    
    # Write same dict in a more readable format (new line per key)
    with open (path1 + '\\' + outfile2, 'wt') as out2:   
        for key, value, in dict_R_SPnum.items():
            out2.write('{0}: {1}\n'.format(key, value))
        
# Forward simple dict
if save_plaindict == 1:
    outfile_SIFTS = filenameSIFTS[:-4] + '_dict.txt'
    outfile_SIFTS_readable = outfile_SIFTS[:-4] + "_readable.txt"
    print('Saving output of dict_SIFTS as ', outfile_SIFTS, outfile_SIFTS_readable)
    
    with open (path1 + '\\' + outfile_SIFTS, 'wt') as out3:
        out3.write(str(dict_SIFTS))   
    
    # Write  dict2 in a more readable format (new line per key)
    with open (path1 + '\\' + outfile_SIFTS_readable, 'wt') as out4:   
        #header = '#HEADER: HOLO_chain SP_BEG SP_END : APO_chain SP_BEG SP_END SP_LEN overlap_len %overlap\n'
        out4.write(header1)
        out4.write(header2)
        for key, value, in dict_SIFTS.items():
            out4.write('{0}: {1}\n'.format(key, value))

print('Done')