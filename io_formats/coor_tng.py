from frame import *
from numpy import *
from ctypes import *
from numpy.ctypeslib import ndpointer

import os
cwd=os.path.dirname(os.path.abspath(__file__))

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame

try: 
    tng=cdll.LoadLibrary(cwd+'/libtngfile.so')
except:
    raise IOError("libtngfile.so can't be loaded")

def open_traj_read(file_name):

    io_vars=[0 for ii in range(30)]
    io_pos=0
    io_err=0

    funit=tng.TrajngOpenRead(file_name)
    if not funit:
        io_err=1
    else:
        natoms=tng.TrajngNatoms(funit)
        io_vars[0]=natoms

    return funit,io_vars,io_pos,io_err

def read_next(file_unit,io_vars=None,io_pos=None):

    io_err=0
    io_end=0

    

    while 1

def nose1(file_unit):
    return tng.TrajngGetProgramInfo(file_unit)

def nose2(file_unit):
    return tng.TrajngGetAtomLabels(file_unit)

def close_traj(file_unit):

    io_err=0
    tng.TrajngClose(file_unit)
