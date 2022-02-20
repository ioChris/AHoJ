import unittest

import prepare
import apoholo_J
import os
import shlex
from common import get_default_workdir


# TODO always test in clean temporary workdir e.g. 'tmp_test/workdir'
workdir = get_default_workdir()


def exists_not_empty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


class PrepareTests(unittest.TestCase):

    def test_prepare(self):
        exit_code = prepare.main([])
        self.assertEqual(0, exit_code)
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_dict.txt'))
        self.assertTrue(exists_not_empty(workdir + '/SIFTS/pdb_chain_uniprot_REVERSE_SPnum.txt'))


class ApoholoTests(unittest.TestCase):

    def do_test(self, args_str, expect_failure=False):
        argv = shlex.split(args_str)
        print("Testing with args:", argv)
        exit_code = apoholo_J.main(argv)
        if expect_failure:
            self.assertNotEqual(0, exit_code)
        else:
            self.assertEqual(0, exit_code)
            # TODO test correctness of actual output

    def test_successful_runs(self):
        self.do_test("--save_oppst 1 --save_separate 1 --query '1a73 A,B ZN' ")
        self.do_test("--save_oppst 1 --save_separate 1 --query '1a73 ALL ZN' ")
        self.do_test("--save_oppst 1 --save_separate 1 --query '1a73 ZN' ")
        self.do_test("--save_oppst 1 --save_separate 1 --query '1a73' ")
        self.do_test("--save_oppst 1 --save_separate 1 --query '1a73 A' ")
        self.do_test("--save_oppst 1 --save_separate 1 --query '1a73 A ZN,MG' ")
        # TODO add more tests cases

    def test_expected_failures(self):
        self.do_test("--invalid_param ",               expect_failure=True)
        self.do_test("--query 'INVALID_QUERY X X X' ", expect_failure=True)
        # TODO add more tests cases


if __name__ == '__main__':
    unittest.main()
