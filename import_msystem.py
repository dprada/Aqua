from cl_set import *

import parmed as pmd
import MDAnalysis as mda

def MDAnalysis(universe):

    temp=msystem()
    
    for atom in universe.atoms:

        temp_atom=cl_unit()
        temp_atom.name=atom.name
        temp_atom.index=atom.id
        temp_atom.pdb_index=atom.index
        temp_atom.resid.name=atom.resname
        temp_atom.resid.pdb_index=atom.resid
        
        temp.atom.append(temp_atom)

    return temp
