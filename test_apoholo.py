import unittest

import prepare
import apoholo
import os
import sys
import shlex

from apoholo import load_precompiled_data, Query, QueryResult, parse_query
from common import get_default_workdir, read_file

# TODO always run tests with clean temporary workdir e.g. 'tmp_test/workdir'
workdir = get_default_workdir()


def exists_not_empty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def non_blank_lines(path):
    """
    Number of non-blank lines in a text file, returns 0 if file doesn't exist
    """
    count = 0
    if os.path.exists(path):
        with open(path) as file:
            for line in file:
                if line.strip():
                    count += 1
    return count


def count_files(dir, suffix):
    return len([f for f in os.listdir(dir) if f.endswith(suffix) and exists_not_empty(os.path.join(dir, f))])


# T01 prefix to ensure it is run first
# TODO find better way to enforce test order
class T01_Prepare(unittest.TestCase):

    def test_prepare_script(self):
        exit_code = prepare.main([])
        self.assertEqual(0, exit_code)
        #self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_dict.txt'))
        #self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_REVERSE_SPnum.txt'))
        #self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_dict.bin'))
        #self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_REVERSE_SPnum.bin'))
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/uniprot_segments_observed_dict.txt'))
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/uniprot_segments_observed_REVERSE_SPnum.txt'))
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/uniprot_segments_observed_dict.bin'))
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/uniprot_segments_observed_REVERSE_SPnum.bin'))


class T02_Apoholo(unittest.TestCase):

    def setUp(self):
        self.precompiled_data = load_precompiled_data(workdir)

    # Test successful query
    def tst_query(self, args_str, expect_apo=0, expect_holo=0) -> QueryResult:
        argv = shlex.split(args_str)  # split but preserve '...' as substrings
        print("Testing with args:", argv)

        args = apoholo.parse_args(argv)
        res = apoholo.process_query(args.query, workdir, args, self.precompiled_data)

        print("Query result:", res)

        if expect_apo > 0:
            assert res.num_apo_chains >= expect_apo, f"Found less APO chains ({res.num_apo_chains} vs {expect_apo})"  # may find more due to future data
            assert non_blank_lines(res.result_dir + '/results_apo.csv') == 1 + expect_apo, "Bad results_apo.csv"

        if expect_holo > 0:
            assert res.num_holo_chains >= expect_holo, f"Found less HOLO chains ({res.num_holo_chains} vs {expect_holo})"  # may find more due to future data
            assert non_blank_lines(res.result_dir + '/results_holo.csv') == 1 + expect_holo, "Bad results_holo.csv"

        # Test produced structure files
        if expect_apo + expect_holo > 0:
            assert count_files(res.result_dir + '/structure_files', '.cif.gz') >= 1 + expect_apo + expect_holo, "Failed to produce right number of .cif.gz files"

        return res

    # TODO instead of running main, test by running subprocess and capture stdout/stderr, test for presence of expected error messages
    def tst_main_fail(self, args_str):
        argv = shlex.split(args_str)  # split but preserve '...' as substrings
        print("Testing with args:", argv)
        with self.assertRaises(BaseException):
            exit_code = apoholo.main(argv)

    def test_parse_query(self):
        assert parse_query('1a73',           autodetect_lig=True,  water_as_ligand_auto=False, nonstd_rsds_as_lig_auto=False, d_aa_as_lig_auto=False) == Query(struct='1a73', chains='ALL', ligands=None,   position=None, autodetect_lig=1, water_as_ligand_auto=0, nonstd_rsds_as_lig_auto=0, d_aa_as_lig_auto=0)
        assert parse_query('1a73 A,B ZN',    autodetect_lig=False, water_as_ligand_auto=False, nonstd_rsds_as_lig_auto=False, d_aa_as_lig_auto=False) == Query(struct='1a73', chains='A,B', ligands='ZN',   position=None, autodetect_lig=0, water_as_ligand_auto=0, nonstd_rsds_as_lig_auto=0, d_aa_as_lig_auto=0)
        assert parse_query('1a73 ALL ZN,MG', autodetect_lig=False, water_as_ligand_auto=False, nonstd_rsds_as_lig_auto=False, d_aa_as_lig_auto=False) == Query(struct='1a73', chains='ALL', ligands='ZN,MG',position=None, autodetect_lig=0, water_as_ligand_auto=0, nonstd_rsds_as_lig_auto=0, d_aa_as_lig_auto=0)
        assert parse_query('1a73 A',         autodetect_lig=True,  water_as_ligand_auto=False, nonstd_rsds_as_lig_auto=False, d_aa_as_lig_auto=False) == Query(struct='1a73', chains='A',   ligands=None,   position=None, autodetect_lig=1, water_as_ligand_auto=0, nonstd_rsds_as_lig_auto=0, d_aa_as_lig_auto=0)

    def test_simple_queries(self):
        self.tst_query("--query '1a73 A,B ZN' ",  expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73 ALL ZN' ",  expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73 * ZN' ",    expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73' ",         expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73 A' ",       expect_apo=0, expect_holo=16)
        self.tst_query("--query '1a73 A ZN,MG' ", expect_apo=0, expect_holo=16)

    def test_advanced_queries(self):
        self.tst_query("--query '3CQV A HEM' ",      expect_apo=6,  expect_holo=5)
        self.tst_query("--query '3fav all zn' ",     expect_apo=2,  expect_holo=0)
        self.tst_query("--query '2hka all c3s' ",    expect_apo=2,  expect_holo=0)  # bovine NPC2 complex with cholesterol sulfate
        self.tst_query("--query '2v57 A,C prl' ",    expect_apo=4,  expect_holo=0)  # SS changes in transcriptional regulator LfrR in complex with proflavine
        self.tst_query("--query '1fmk A HOH 1011' ", expect_apo=18, expect_holo=3)  # (hard target) Water molecule in the interface of two domains that undergo extensive conformational changes upon ligand binding (v0.4.6 apo 25)
        self.tst_query("--query '1aro P HG 904' ",   expect_apo=45, expect_holo=0)  # (previously apo 19-21) fragmented UniProt candidates, to use for testing UNP overlap calculation
        self.tst_query("--query '1pkz A tyr 9' --water_as_ligand 1",   expect_apo=9, expect_holo=88) # marian example, allosteric effects of HOH

    def test_advanced_queries_B(self):
        self.tst_query("--query '1a37 P sep 259' ",     expect_apo=6,  expect_holo=4)  # Phosphoserine on non-uniprot chain (hard target)
        self.tst_query("--query '1a37 P sep' ",         expect_apo=6,  expect_holo=4)  # same without specifying position
        self.tst_query("--query '1gb1' ",               expect_apo=22, expect_holo=16) # v0.4.6 apo 20, holo 16. Apo NMR structure with solid state NMR candidates (threw error in previous versions)

    def test_broad_search(self):
        self.tst_query("--query '2v0v' ",     expect_apo=8, expect_holo=24)    # test for reverse search (this is a fully apo structure)
        self.tst_query("--query '2v0v A,B' ", expect_apo=4, expect_holo=12)    #

    def test_interface_ligands_search(self):
        self.tst_query("--query '1a73 E mg 205' ", expect_apo=4, expect_holo=12) # non-protein query chain
        self.tst_query("--query '1a73 A mg' ",     expect_apo=4, expect_holo=12) # interface ligand assigned nucleic acid chain in the PDB file
        self.tst_query("--query '6XBY A adp,mg' ", expect_apo=7, expect_holo=2)  # ADP is an interface ligand on non-query chain
        self.tst_query("--query '3k0v B NAG 1' ",  expect_apo=2, expect_holo=86) # ligand is on a non-uniprot chain that is remapped, it was not detected in previous versions (< 0.4.9)

    def test_lig_bndng_chains_parsing(self):
        self.tst_query("--query '1a73 ! MG' ",    expect_apo=8, expect_holo=24)     # "!" means that only query chains that bind the ligand(s) will be considered for the search
        self.tst_query("--query '3vro ! PTR' ",   expect_apo=7, expect_holo=0)
        self.tst_query("--query '3n7y ! PTR' ",   expect_apo=21, expect_holo=270)   # local-only error (python 3.9.9): total number of cif files is less than total apo + holo results [case-insensitivity in Windows overwrites upper/lower case chains resulting in less cif files (e.g, 5cdwC, 5cdwc)]
        self.tst_query("--query '1a73 ! MG,ZN' ", expect_apo=0, expect_holo=32)

    def test_expected_failures(self):
        self.tst_main_fail("--invalid_param ")
        self.tst_main_fail("--query 'INVALID_QUERY X X X' ")
        self.tst_main_fail("--query 'XXXX' ")                # XXXX is not in the PDB
        self.tst_main_fail("--query '1a73 \n XXXX' ")        # should fail if at least one is not valid
        self.tst_main_fail("--query '1apm E sep,ser' ")      # standard residue without positional argument, should fail
        self.tst_main_fail("--query '5aep ! hem' ")          # "!" operator with ligand that does not bind any query chain
        # TODO add more tests cases


    def assert_progress(self, res, expected_content=None):
        progress_file = res.result_dir + '/.progress'
        assert os.path.exists(progress_file), ".progress file was not created"
        assert non_blank_lines(progress_file) == 1, ".progress doesn't have exactly one line"
        content = read_file(progress_file)
        assert content == expected_content, f".progress file has invalid value. actual: '{content}' expected: '{expected_content}'"


    def test_track_progress(self):
        res = self.tst_query("--query '3CQV A HEM' --track_progress 1", expect_apo=6, expect_holo=5)
        self.assert_progress(res, expected_content="11/11")
        res = self.tst_query("--query '3fav all zn' --track_progress 1",  expect_apo=2, expect_holo=0)
        self.assert_progress(res, expected_content="2/2")  # change expected content from "4/4")


if __name__ == '__main__':
    sys.exit(unittest.main())
