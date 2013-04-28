from 
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
    tng=cdll.LoadLibrary(cwd+'/libtrajng.so')
except:
    raise IOError("libtrajng.so can't be loaded")


