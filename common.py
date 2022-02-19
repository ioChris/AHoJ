import os
import sys


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
    """ Returns the directory the current script (or interpreter) is running in  """
    path = os.path.realpath(sys.argv[0])
    if os.path.isdir(path):
        return path
    else:
        return os.path.dirname(path)


def get_default_workdir():
    if os.name == 'nt':
        return get_2level_cwd() + '/Documents/Bioinfo_local/Ions/datasets_local/APO_candidates/webserver'
    else:
        return os.path.join(get_script_directory(), '..', 'ahoj_workdir')


def get_workdir(args):
    # TODO add override from local config
    if args.work_directory is not None:
        return args.work_directory
    else:
        return get_default_workdir()

