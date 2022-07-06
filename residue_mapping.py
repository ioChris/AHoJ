# -*- coding: utf-8 -*-
"""
Created on Fri May 13 18:33:16 2022

@author: ChrisH
"""

''' Functions for residue mapping '''


from lxml import etree


def print_dict_readable(input_dict, header_msg):
    #print('')
    print(header_msg)
    for i, j in input_dict.items():
        print(i, j)


def map_uniprot_resnum_to_pdb2(uniprot_resnum_list, sifts_xml_file):

    # Load the xml with lxml
    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(sifts_xml_file, parser)
    root = tree.getroot()
    my_pdb_resnum = None
    # TODO: "Engineered_Mutation is also a possible annotation, need to figure out what to do with that
    my_pdb_annotation = False

    # Parse list of residues and their chains
    list_out = list()            
    for resichain in uniprot_resnum_list:
        resichain.replace('_', '.')
        chain_id = resichain.split('.')[0]
        uniprot_resnum = resichain.split('.')[1]
        print(chain_id, uniprot_resnum)

        # Find the right chain (entities in the xml doc)
        ent = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity'
        for chain in root.findall(ent):
            # IMPORTANT - entityId is not the chain ID
            if chain.attrib['entityId'] == chain_id:
                # Find the "crossRefDb" tag that has the attributes dbSource="UniProt" and  dbResNum="your_resnum_here"
                # Then match it to the crossRefDb dbResNum that has the attribute dbSource="PDBresnum"

                # Check if uniprot + resnum even exists in the sifts file (it won't if the pdb doesn't contain the residue)
                ures = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="UniProt"][@dbResNum="%s"]' % uniprot_resnum
                my_uniprot_residue = chain.findall(ures)

                if len(my_uniprot_residue) == 1: # this is prolly a list and we need it to have 1 element
                    # Get crossRefDb dbSource="PDB"
                    parent = my_uniprot_residue[0].getparent()
                    pres = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="PDB"]'
                    my_pdb_residue = parent.findall(pres)
                    my_pdb_resnum = int(my_pdb_residue[0].attrib['dbResNum'])
                    my_pdb_chain = str(my_pdb_residue[0].attrib['dbChainId'])
                    resichain_pdb = my_pdb_chain + '.' + str(my_pdb_resnum)

                    # Get <residueDetail dbSource="PDBe" property="Annotation">
                    # Will be Not_Observed if it is not seen in the PDB
                    anno = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail[@dbSource="PDBe"][@property="Annotation"]'
                    my_pdb_annotation = parent.findall(anno)
                    if len(my_pdb_annotation) == 1:
                        my_pdb_annotation = my_pdb_annotation[0].text
                        if my_pdb_annotation == 'Not_Observed':
                            my_pdb_annotation = False
                            list_out.append(resichain + '>*' + resichain_pdb)
                    else:
                        my_pdb_annotation = True
                        list_out.append(resichain + '>' + resichain_pdb)
                else:
                    list_out.append(resichain + '>None')#  + None)
    return list_out


def map_pdb_resnum_to_uniprot(pdb_resnum_list, sifts_xml_file):

    # Load the xml with lxml
    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(sifts_xml_file, parser)
    root = tree.getroot()
    my_unp_resnum = None
    # TODO: "Engineered_Mutation is also a possible annotation, need to figure out what to do with that
    my_pdb_annotation = False

    # Parse list of residues and their chains
    list_out = list()            
    for resichain in pdb_resnum_list:
        segi_id = resichain.split('_')[0]
        chain_id = resichain.split('_')[1]
        pdb_resnum = resichain.split('_')[3]
        resichain_frmt = '.'.join([segi_id, chain_id, pdb_resnum]) # Reformat the string to more compact form
        #print(segi_id, chain_id, pdb_resnum)

        # Find the right segment (entities in the xml doc)
        ent = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity'
        for segment in root.findall(ent):
            # IMPORTANT - entityId is not the chain ID
            if segment.attrib['entityId'] == segi_id:
                # Find the "crossRefDb" tag that has the attributes dbSource="UniProt" and  dbResNum="your_resnum_here"
                # Then match it to the crossRefDb dbResNum that has the attribute dbSource="PDBresnum"

                # Check if uniprot + resnum even exists in the sifts file (it won't if the pdb doesn't contain the residue)
                ures = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="PDB"][@dbResNum="%s"]' % pdb_resnum
                my_pdb_residue = segment.findall(ures)

                if len(my_pdb_residue) == 1: # this is prolly a list and we need it to have 1 element
                    # Get crossRefDb dbSource="PDB"
                    parent = my_pdb_residue[0].getparent()
                    pres = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="UniProt"]'
                    my_uniprot_residue = parent.findall(pres)
                    try:
                        my_unp_resnum = int(my_uniprot_residue[0].attrib['dbResNum'])
                    except Exception: # Ignore atom residues annotated in non-protein chains (i.e. 2amq)
                        continue
                    #my_pdb_chain = str(my_uniprot_residue[0].attrib['dbChainId']) # UniProt field has no chain in xml
                    my_pdb_chain = chain_id
                    resichain_pdb = my_pdb_chain + '.' + str(my_unp_resnum)

                    # Get <residueDetail dbSource="PDBe" property="Annotation">
                    # Will be Not_Observed if it is not seen in the PDB
                    anno = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail[@dbSource="PDBe"][@property="Annotation"]'
                    my_pdb_annotation = parent.findall(anno)
                    if len(my_pdb_annotation) == 1:
                        my_pdb_annotation = my_pdb_annotation[0].text
                        if my_pdb_annotation == 'Not_Observed':
                            my_pdb_annotation = False
                            list_out.append(resichain_frmt + '>*' + resichain_pdb) # Not observed in structure
                    else:
                        my_pdb_annotation = True
                        list_out.append(resichain_frmt + '>' + resichain_pdb)
                else:
                    list_out.append(resichain_frmt + '>None')#  + None)
    return list_out


# Group (UNP num) binding residues by chain (list in, dict out)
def group_mapped_res_by_chain(mapped_res_list):
    bndgres_pdb_to_unp_chains = dict()
    for unp_binding_res in mapped_res_list:
        unp_part = unp_binding_res.split('>')[1]
        chain = unp_part.split('.')[0]
        unpnum = unp_part.split('.')[1]
        bndgres_pdb_to_unp_chains.setdefault(chain, []).append(unpnum)
    return bndgres_pdb_to_unp_chains


# Find whether mapped binding residues are present in each candidate chain
def examine_cndt_mapped_bs_res(dict_of_bndgres_pdb_to_unp, query_structchain, candidates_unp_dict):
    candidate_hits = dict()
    query_struct = query_structchain[:4]
    query_chain = query_structchain[4:]

    for chain, positions in dict_of_bndgres_pdb_to_unp.items():
        if chain == query_chain: # Only examine residues of query structchain (ignore interface chains/resis)
            candidates = candidates_unp_dict[query_structchain] # Get candidates from UniProt ID

            for candidate_entry in candidates:
                candidate_structchain = candidate_entry.split()[0]
                candidate_struct = candidate_structchain[:4]

                if candidate_struct != query_struct: # eliminate candidates of the same structure

                    # Make dict with candidate structchain and UniProt range(s) as values, then check overlap
                    cndt_SP_BEG = int(candidate_entry.split()[1])
                    cndt_SP_END = int(candidate_entry.split()[2])

                    # Loop through positions, look if position is within UniProt range
                    for position in positions:
                        if int(position) >= cndt_SP_BEG and int(position) <= cndt_SP_END:
                            candidate_hits.setdefault(candidate_structchain+'.'+query_structchain, []).append(chain + '.' + position + ' ' + str(1))
                        else:
                            candidate_hits.setdefault(candidate_structchain+'.'+query_structchain, []).append(chain + '.' + position + ' ' + str(0))
    # Remove duplicates
    for key, value in candidate_hits.items():
        candidate_hits[key] = list(candidate_hits.fromkeys(value))

    return candidate_hits


# Intermediate step to remove negative score when positive is present for the same position
def remove_negative_duplicate_cndt_bs_res_pos(dict_of_candidate_bs_rsds_assessment):
    candidate_metahits = dict()
    for candidate, scores in dict_of_candidate_bs_rsds_assessment.items():
        score_list = list()
        for score in scores:
            position = score.split()[0]
            result = score.split()[1]
            bad_result = position + ' ' + '0'
            good_result = position + ' ' + '1'
            if bad_result in score_list and result == '1':
                score_list.remove(bad_result)
                score_list.append(position + ' ' + result)
            elif good_result in score_list and result == '0':
                pass
            else:
                score_list.append(position + ' ' + result)
        candidate_metahits[candidate] = score_list
    return candidate_metahits


# Count number of mapped binding residues (out of total) present in candidate
def evaluate_candidate_bs_rsds(dict_of_candidate_bs_rsds_scores):
    # dict e.g. "1lbaA.1aroL ['L.19 1', 'L.20 1', 'L.131 1']"    
    candidate_scores = dict()
    for candidate, scores in dict_of_candidate_bs_rsds_scores.items():
        current_score = 0  # reset
        for score in scores:
            #position = score.split()[0]
            result = score.split()[1]
            if result == '1':
                current_score += 1
        ratio = 100 * current_score / len(scores)
        ratio = str(round(ratio)) + '%'
        candidate_scores.setdefault(candidate, []).append(str(current_score) + '/' + str(len(scores)) + ' ' + ratio)
    return candidate_scores


# Put candidates over certain threshold to dict for further processing (applies on precalculated % scores)
def good_candidates_from_residue_mapping(candidate_score_dict, binding_residue_threshold):
    dict_rsd_map_candidates = dict()
    for key, value in candidate_score_dict.items():
        cndt_structchain_part = key.split('.')[0]
        qr_structchain_part = key.split('.')[1]
        brsds_percent = int(value[0].split()[1][:-1])
        if brsds_percent >= binding_residue_threshold:
            #print(cndt_structchain_part, brsds_percent)
            dict_rsd_map_candidates.setdefault(qr_structchain_part, []).append(cndt_structchain_part)
    return dict_rsd_map_candidates


# Put candidates under certain threshold to dict for further processing (applies on precalculated % scores)
def bad_candidates_from_residue_mapping(candidate_score_dict, binding_residue_threshold):
    bad_candidates = dict()
    for key, value in candidate_score_dict.items():
        cndt_structchain_part = key.split('.')[0]
        qr_structchain_part = key.split('.')[1]
        brsds_percent = int(value[0].split()[1][:-1])
        if brsds_percent < binding_residue_threshold:
            #print(cndt_structchain_part, brsds_percent)
            bad_candidates.setdefault(qr_structchain_part, []).append(cndt_structchain_part)
    return bad_candidates


def get_scores_from_residue_mapping(candidate_score_dict, query_struchain, candidate_structchain):
    key_expression = candidate_structchain + '.' + query_struchain
    try:
        scores = candidate_score_dict[key_expression]
        scores = scores[0]
    except Exception:
        scores = '-'
    return scores
