import unittest

import prepare
import apoholo_J
import os
import sys
import shlex
from common import get_default_workdir


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
class T01_PrepareScriptTests(unittest.TestCase):

    def test_prepare_script(self):
        exit_code = prepare.main([])
        self.assertEqual(0, exit_code)
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_dict.txt'))
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_REVERSE_SPnum.txt'))


class T02_ApoholoTests(unittest.TestCase):

    def tst_query(self, args_str, expect_apo=0, expect_holo=0):
        argv = shlex.split(args_str)  # split but preserve '...' as substrings
        print("Testing with args:", argv)

        args = apoholo_J.parse_args(argv)
        res = apoholo_J.process_query(args.query, workdir, args)

        print("Query result:", res)

        if expect_apo > 0:
            assert res.num_apo_chains >= expect_apo, "Found less APO chains"  # may be > due to future data
            assert non_blank_lines(res.result_dir + '/results_apo.csv') == 1 + expect_apo, "Bad results_apo.csv"

        if expect_holo > 0:
            assert res.num_holo_chains >= expect_holo, "Found less HOLO chains"  # may be > due to future data
            assert non_blank_lines(res.result_dir + '/results_holo.csv') == 1 + expect_holo, "Bad results_holo.csv"

        # Test produced structure files
        if expect_apo + expect_holo > 0:
            assert count_files(res.result_dir, '.cif.gz') >= 1 + expect_apo + expect_holo, "Failed to produce right number of .cif.gz files"

    # TODO instead of running main, test by running subprocess and capture stdout/stderr, test for presence of expected error messages
    def tst_main_fail(self, args_str):
        argv = shlex.split(args_str)  # split but preserve '...' as substrings
        print("Testing with args:", argv)
        with self.assertRaises(SystemExit):
            exit_code = apoholo_J.main(argv)
            # self.assertNotEqual(0, exit_code)

    def test_successful_runs(self):
        self.tst_query("--query '1a73 A,B ZN' ",  expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73 ALL ZN' ",  expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73 ZN' ",      expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73' ",         expect_apo=0, expect_holo=32)
        self.tst_query("--query '1a73 A' ",       expect_apo=0, expect_holo=0)   # TODO really 0?
        self.tst_query("--query '1a73 A ZN,MG' ", expect_apo=0, expect_holo=16)

        self.tst_query("--query '3CQV A HEM' ",   expect_apo=6, expect_holo=5)
        self.tst_query("--query '3fav all zn' ",  expect_apo=2, expect_holo=0)
        self.tst_query("--query '2hka all c3s' ", expect_apo=2, expect_holo=0)  # bovine NPC2 complex with cholesterol sulfate
        self.tst_query("--query '2v57 a,c prl' ", expect_apo=4, expect_holo=0)  # SS changes in transcriptional regulator LfrR in complex with proflavine

        # TODO add expected numbers
        self.tst_query("--reverse_search 1 --query '2v0v' ")   # test for reverse search (this is a fully apo structure)
        self.tst_query("--reverse_search 1 --query '2v0v a,b' ")

    def test_expected_failures(self):
        self.tst_main_fail("--invalid_param ")
        self.tst_main_fail("--query 'INVALID_QUERY X X X' ")
        # self.do_test("--query 'XXXX' ",                expect_failure=True)  # XXXX is not in PDB    TODO fix failing test
        # TODO add more tests cases


    #def tst_main(self, args_str, expect_failure=False):
    #    argv = shlex.split(args_str)
    #    print("Testing with args:", argv)
    #
    #    if expect_failure:
    #        with self.assertRaises(SystemExit):
    #            exit_code = apoholo_J.main(argv)
    #            # self.assertNotEqual(0, exit_code)
    #    else:
    #        exit_code = apoholo_J.main(argv)
    #        self.assertEqual(0, exit_code)
    #        # TODO test correctness of actual output

    # test if github action fails correctly on test error
    # def test_fail_1(self):
    #     self.assertEqual(1, 2)


if __name__ == '__main__':
    sys.exit(unittest.main())
