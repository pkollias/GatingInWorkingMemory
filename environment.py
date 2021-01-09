from os import path

class Environment:

    fsroot = '/Volumes'

    preproc_src_env = path.join(fsroot, 'buschman/Users/Pavlos/Code')
    preproc_dat_env = path.join(fsroot, 'buschman/Projects/GatingInWorkingMemory/Data')
    preproc_dest_env = path.join(fsroot, 'buschman/Projects/GatingInWorkingMemory/db/IO')
    proc_dest_env = path.join(fsroot, 'buschman/Projects/GatingInWorkingMemory/db/IO')
