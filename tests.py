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


# T01 prefix to ensure it is run first
# TODO find better way to enforce test order
class T01_PrepareScriptTests(unittest.TestCase):

    def test_prepare_script(self):
        exit_code = prepare.main([])
        self.assertEqual(0, exit_code)
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_dict.txt'))
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_REVERSE_SPnum.txt'))


class T02_ApoholoTests(unittest.TestCase):

    def do_test(self, args_str, expect_failure=False):
        argv = shlex.split(args_str)
        print("Testing with args:", argv)

        if expect_failure:
            with self.assertRaises(SystemExit):
                exit_code = apoholo_J.main(argv)
                # self.assertNotEqual(0, exit_code)
        else:
            exit_code = apoholo_J.main(argv)
            self.assertEqual(0, exit_code)
            # TODO test correctness of actual output

    def test_successful_runs(self):
        self.do_test("--query '1a73 A,B ZN' ")
        self.do_test("--query '1a73 ALL ZN' ")
        self.do_test("--query '1a73 ZN' ")
        self.do_test("--query '1a73' ")
        self.do_test("--query '1a73 A' ")
        self.do_test("--query '1a73 A ZN,MG' ")

        self.do_test("--query '3CQV A HEM' ")
        self.do_test("--query '3fav all zn' ")
        self.do_test("--query '2hka all c3s' ")  # bovine NPC2 complex with cholesterol sulfate
        self.do_test("--query '2v57 a,c prl' ")  # SS changes in transcriptional regulator LfrR in complex with proflavine

        self.do_test("--reverse_search 1 --query '2v0v' ")   # test for reverse search (this is a fully apo structure)
        self.do_test("--reverse_search 1 --query '2v0v a,b' ")

    def test_expected_failures(self):
        self.do_test("--invalid_param ",               expect_failure=True)
        self.do_test("--query 'INVALID_QUERY X X X' ", expect_failure=True)
        # self.do_test("--query 'XXXX' ",                expect_failure=True)  # XXXX is not in PDB    TODO fix failing test
        # TODO add more tests cases
        # TODO test presence of error messages

    # test if github action fails correctly on test error
    # def test_fail_1(self):
    #     self.assertEqual(1, 2)


if __name__ == '__main__':
    sys.exit(unittest.main())
