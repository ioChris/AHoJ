import os
import sys
import pickle


def get_2level_cwd():
    npath = os.path.normpath(os.getcwd())   # Normalize the path string for the OS
    path0 = os.path.join(npath.split(os.sep)[0], '/', npath.split(os.sep)[1], npath.split(os.sep)[2])
    if os.path.exists(path0):
        memo = "Root path found >> " + path0
    else:
        memo = 'Error finding root path in working dir:' + npath
    print(memo)
    return path0


def get_script_directory():
    """ Get the directory the current script (or interpreter) is running in  """
    # path = os.path.realpath(sys.argv[0]) # doesn't work in tests
    path = os.path.realpath(__file__)      # doesn't work in the interpreter
    if os.path.isdir(path):
        return path
    else:
        return os.path.dirname(path)


def get_default_workdir():
    if os.name == 'nt':
        return get_2level_cwd() + '/Documents/Bioinfo_local/Ions/datasets_local/APO_candidates/webserver'
    else:
        return os.path.normpath(os.path.join(get_script_directory(), '..', 'ahoj_workdir'))


def get_workdir(args):
    """ Get path to global work directory  """
    if args.work_directory is not None:
        return args.work_directory
    else:
        return get_default_workdir()
    # TODO add an override from local config or ENV variable


def load_dict_binary(path):
    return pickle.load(open(path, "rb"))


def save_dict_binary(dict, path):
    return pickle.dump(dict, open(path, "wb"))
