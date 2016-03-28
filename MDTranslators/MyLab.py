import MDTranslators as MDTranslators
import MyLab_geometry as geometry
import MyLab_math as math
#from MyLab_clases import *

def typename(atype):

    #if not isinstance(atype, type):
    #    raise Exception('Argument is not a type')

    modulename = atype.__module__
    try:
        typename   = atype.__name__
    except:
        typename   = atype.__class__.__name__

    if modulename != '__builtin__':
        typename = modulename + '.' + typename

    return typename

def info(obj):
    
    if typename(obj)=='MDAnalysis.core.AtomGroup.Universe':

        print '# TYPE: Universe MDAnalysis'
        print '#','System created from the file',None
        print '#',len(obj.atoms),' atoms'
        print '#',len(obj.residues),' residues'
        print '#',None,' chains'
        print '#',len(obj.segments),' segments'
        print '#',None,' waters'
        print '#',None,' ions'
        print '#',None,' lipids'
        print '#',None,' proteins'
        print '#',None,' ligands'
        print '#',None,' other molecules'

    else:

        print '# Object not recognised'


