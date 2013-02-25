##############################
######## Internal variables:


###############################
####### External modules:

from os import system
from os import path
import copy
from numpy import *
import pickle as pic
import sys

workdir=path.dirname(path.abspath(__file__))
sys.path.append(workdir+'/top_par/')
sys.path.append(workdir+'/io_formats/')

##############################
######## Internal modules:

# Base:
from pyn_cl_set import *
from pyn_cl_coors import *
import pyn_fort_general as f
from pyn_cl_anal_trajs import *
from pyn_cl_kin_anal import *

# Water:
#from pyn_cl_water import *

# Networks:
from pyn_cl_net import *
#import pyn_fort_net as f_net

# GNM and ANM:
#from pyn_cl_enm import *
#import pyn_fort_enm as f_enm

# General Hbonds:
#from pyn_cl_hbonds import *
