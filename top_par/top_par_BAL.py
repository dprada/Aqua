############################################################
### Template to create a new topology with new parameters. #
### ------------------------------------------------------ #
###   This is an example with the molecule Fluoroethanol   #
###              thanks to Tristan Bereau                  #
############################################################

### Name of molecule in pdb and type (residue_types.py)

residue_name='BAL'        # Same residue or molecule name found in your pdb or gro file. 
residue_type='Ion'   # Molecule, Ligand, Protein, Water, Ion, DNA, ...

### Atoms:

atoms={
#'atom_name' : 'atom_type' ,
'BAL'      : 'BAL'         
}

charge={
#'atom_name' : 'charge'    ,
'BAL'         :   0.0   
}

donors=[

]

acceptors=[
]

### Bonds:

covalent_bonds=[
]

