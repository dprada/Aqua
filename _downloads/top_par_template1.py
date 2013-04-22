############################################################
### Template to create a new topology with new parameters. #
### ------------------------------------------------------ #
###   This is an example with the molecule Fluoroethanol   #
###              thanks to Tristan Bereau                  #
############################################################

### Name of molecule in pdb and type (residue_types.py)

residue_name='LIG'        # Same residue or molecule name found in your pdb or gro file. 
residue_type='Molecule'   # Molecule, Ligand, Protein, Water, Ion, DNA, ...

### Atoms:

atoms={
#'atom_name' : 'atom_type' ,
'C1'      : 'C'         ,
'C2'      : 'C'         ,
'H1'      : 'H'         ,
'H2'      : 'H'         ,
'H3'      : 'H'         ,
'H4'      : 'H'         ,
'H5'      : 'H'         ,
'O1'      : 'O'         ,
'F1'      : 'F'         
}

charge={
#'atom_name' : 'charge'    ,
'C1'         :   0.306055  , 
'C2'         :   0.115513  , 
'F1'         :  -0.258437  , 
'O1'         :  -0.640003  , 
'H1'         :  -0.019003  , 
'H2'         :  -0.019003  , 
'H3'         :   0.055159  , 
'H4'         :   0.055159  , 
'H5'         :   0.404558 
}

donors=[
'O1'
]

acceptors=[
'O1',
'F1'
]

### Bonds:

covalent_bonds=[
['C1'   ,'C2'    ],
['C1'   ,'O1'    ],
['C1'   ,'H1'    ],
['C1'   ,'H2'    ],
['C2'   ,'F1'    ],
['C2'   ,'H3'    ],
['C2'   ,'H4'    ],
['O1'   ,'H5'    ]
]

