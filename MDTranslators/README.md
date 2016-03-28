# lab
ESTO ES UN PROYECTO DE LAB PY

##
Hay que estandarizar los nombres de las cosas... para que cuando uno recurre a viejos scripts
pueda cortar y pegar con menos riesgo. Por ejemplo:

    topfile  = U_top_trj_out[0]
    coorfile = U_top_trj_out[1]
    outfile  = U_top_trj_out[2]
    U_mda    = mda.Universe(topfile)
    U        = aa.MDTranslators.MDAnalysis2Aqua(U_mda)
    num_atoms= len(U_mda.atoms)

Tampoco esta demas intentar hacer los parrafos modulares, de tal manera que en el corta pega sea facil enchufar las cosas.

----------------
# Standard:

import MyLab as mylab
import MDAnalysis as mda
import parmed as pmd
import aqua as aqa
