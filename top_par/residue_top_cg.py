### Peptidic Bond
peptidic_bond=['atBB','atBB']

### Alias for residues

residue={
## residues
'ALA'     : 'ALA'    ,
'ARG'     : 'ARG'    ,
'ASN'     : 'ASN'    ,
'ASP'     : 'ASP'    ,
'ASPH'    : 'ASP'    ,
'CYS'     : 'CYS'    ,
'GLU'     : 'GLU'    ,
'GLN'     : 'GLN'    ,
'GLY'     : 'GLY'    ,
'HIS'     : 'HIS'    ,
'HSD'     : 'HIS'    ,
'HSE'     : 'HIS'    ,
'HSP'     : 'HIS'    ,
'HYP'     : 'HYP'    ,
'ILE'     : 'ILE'    ,
'LEU'     : 'LEU'    ,
'LYS'     : 'LYS'    ,
'LYSH'    : 'LYS'    ,
'MET'     : 'MET'    ,
'PHE'     : 'PHE'    ,
'PRO'     : 'PRO'    ,
'SER'     : 'SER'    ,
'THR'     : 'THR'    ,
'TRP'     : 'TRP'    ,
'TYR'     : 'TYR'    ,
'VAL'     : 'VAL'    ,
# terminals
'ACE'     : 'ACE'    ,
'NME'     : 'NME'    ,
'NAC'     : 'NME'    ,
'NHE'     : 'NHE'    ,
'NH2'     : 'NHE'    ,
# lipids
'AOT'     : 'AOT'    ,
# water
'SOL'     : 'SOL3'    ,
'HOH'     : 'SOL3'    ,
'TIP'     : 'SOL3'    ,
'TIP3'    : 'SOL3'    ,
'HO4'     : 'SOL4'    ,
'TIP4'    : 'SOL4'    ,
'HO5'     : 'SOL5'    ,
'TIP5'    : 'SOL5'    ,
'SWM'     : 'SOL5'    ,
# ions
'NA'      : 'NA'      ,
'SOD'     : 'NA'      ,
'K'       : 'K'       ,
'LI'      : 'LI'      ,
'CL'      : 'CL'      ,
'CLA'     : 'CL'      
}

#residue_atoms={'amino':[A,B,C],...}
#covalent_bonds={'amino':[[A,B],[A,C]],...}

residue_atoms={}
covalent_bonds={}
terminal_bonds={}

###############################
###############################

############ PEPTIDES::

######## Aminoacids

### ALA:

residue_atoms['ALA']=[
'atBB'
]
 
covalent_bonds['ALA']=[
]

### ARG:

residue_atoms['ARG']=[
'atBB',
'atSC1',
'atSC2',
]

covalent_bonds['ARG']=[
['atBB'   ,'atSC1'    ], 
['atSC1'  ,'atSC2'    ]
]

### ASN:

residue_atoms['ASN']=[
'atBB',
'atSC1'
]

covalent_bonds['ASN']=[
['atBB'   ,'atSC1'     ]
]

### ASP:

residue_atoms['ASP']=[
'atBB',
'atSC1'
]

covalent_bonds['ASP']=[
['atBB'   ,'atSC1'    ]
]

### CYS and CYSH [HG1]:

residue_atoms['CYS']=[
'atBB',
'atSC1'
]
 
covalent_bonds['CYS']=[
['atBB'   ,'atSC1'    ]
]

### GLU:

residue_atoms['GLU']=[
'atBB',
'atSC1'
]

covalent_bonds['GLU']=[
['atBB'   , 'atSC1'   ]
]

### GLN:

residue_atoms['GLN']=[
'atBB',
'atSC1'
]

covalent_bonds['GLN']=[
['atBB'   , 'atSC1'   ]
]

### GLY:

residue_atoms['GLY']=[
'atBB'
]

covalent_bonds['GLY']=[
]

### HISE (ND1 no H, NE2 with H),
### HISD (ND1 with H, NE2 no H),
### HISH (ND1 with H, NE2 with H),
### All included in HIS:

residue_atoms['HIS']=[
'atBB',
'atSC1',
'atSC2',
'atSC3'
]

covalent_bonds['HIS']=[
['atBB'   ,'atSC1'    ], 
['atSC1'   ,'atSC2'   ], 
['atSC1'  ,'atSC3'   ], 
['atSC2'  ,'atSC3'   ]
]

### HYP:

residue_atoms['HYP']=[
'atN',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atHG1',
'atOD1',
'atHD1',
'atCD2', 
'atHD21',
'atHD22',
'atC',
'atO'
]

covalent_bonds['HYP']=[
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'  ], 
['atCG',  'atOD1'  ],
['atOD1',  'atHD1'  ],
['atCG',  'atCD2'   ], 
['atCD2',  'atHD21' ], 
['atCD2',  'atHD22' ], 
['atCD2',  'atN'   ], 
['atC',   'atO'    ] 
]



### ILE:

residue_atoms['ILE']=[
'atBB', 
'atSC1'
]

covalent_bonds['ILE']=[
['atBB',   'atSC1'    ]
]

### LEU:

residue_atoms['LEU']=[
'atBB',
'atSC1'
]

covalent_bonds['LEU']=[
['atBB',   'atSC1'    ]
]

### LYS with LYSH (HZ3):

residue_atoms['LYS']=[
'atBB',
'atSC1',
'atSC2'
]

covalent_bonds['LYS']=[
['atBB',   'atSC1'   ], 
['atSC1',  'atSC2'  ]
]

### MET:

residue_atoms['MET']=[
'atBB',
'atSC1' 
]

covalent_bonds['MET']=[
['atBB',   'atSC1'    ]
]

### PHE:

residue_atoms['PHE']=[
'atBB',
'atSC1',
'atSC2',
'atSC3'
]

covalent_bonds['PHE']=[
['atBB',   'atSC1'    ], 
['atSC1',   'atSC2'   ], 
['atSC1',  'atSC3'   ], 
['atSC2',  'atSC3'   ]
]

### PRO:

residue_atoms['PRO']=[
'atBB',
'atSC1'
]

covalent_bonds['PRO']=[
['atBB',   'atSC1'   ]
]

### SER:

residue_atoms['SER']=[
'atBB',
'atSC1'
]

covalent_bonds['SER']=[
['atBB',  'atSC1'    ]
]

### THR:

residue_atoms['THR']=[
'atBB',
'atSC1'
]

covalent_bonds['THR']=[
['atBB',   'atSC1'     ]
]

### TRP:

residue_atoms['TRP']=[
'atBB',
'atSC1',
'atSC2',
'atSC3',
'atSC4'
]

covalent_bonds['TRP']=[
['atBB',   'atSC1'   ], 
['atSC1',   'atSC2'  ], 
['atSC2',  'atSC3'  ], 
['atSC1',  'atSC3'  ], 
['atSC2',  'atSC4'   ], 
['atSC3',  'atSC4' ]
]

### TYR:

residue_atoms['TYR']=[
'atBB',
'atSC1',
'atSC2', 
'atSC3' 
]

covalent_bonds['TYR']=[
['atBB',   'atSC1'    ], 
['atSC1',   'atSC2'   ], 
['atSC1',  'atSC3'   ], 
['atSC2',  'atSC3'   ]
]

### VAL:

residue_atoms['VAL']=[
'atBB',
'atSC1'
]

covalent_bonds['VAL']=[
['atBB',   'atSC1'   ]
]


## Terminals

terminal_bonds=[]
#terminal_bonds[]={}



##### LIPIDS:

### AOT:


##### WATER:

### SOL3:


### SOL4:


### SOL5:



##### IONS:

## NA:

## K:

## LI:

## CL:




