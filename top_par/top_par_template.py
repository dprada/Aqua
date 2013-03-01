###########################################################
### Template to create a new topology with new parameters.
###########################################################

## Atom parameters:

atoms={
#'atom_name' : 'atom_type' ,
'HA1'       : 'H'         ,
'HA2'       : 'H'         ,
'OD1'       : 'C'         
}

charge={
#'atom_name' : 'charge'    ,
'HA1'       :  0.2        ,
'OD1'       : -1.2
}

donors=[
#'atom_name',
'HA1'      ,
'HA2'
]

acceptors=[
#'atom_name',
'OD1'
]

## Molecule or residue:

residue={
#'name'      : 'type'
'LIG'       : 'Ligand'
}

atom_list=[
'C1' 
'C2'    
'F1'    
'O1'    
'H1'    
'H2'    
'H3'    
'H4'    
'H5' 
]

covalent_bonds=[
['atN'   ,'atH'    ],
['atN'   ,'atCA'   ],
['atCA'  ,'atHA'   ],
['atCA'  ,'atCB'   ],
['atCA'  ,'atC'    ],
['atCB'  ,'atHB1'  ],
['atCB'  ,'atHB2'  ],
['atCB'  ,'atHB3'  ],
['atC'   ,'atO'    ] 
]

