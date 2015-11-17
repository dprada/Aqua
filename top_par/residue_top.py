### Peptidic Bond
peptidic_bond=['atC','atN']

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
'POPC'    : 'POPC'   ,
'CHL1'    : 'CHL1'   ,
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
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atHB3',
'atC',
'atO'
]
 
covalent_bonds['ALA']=[
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

### ARG:

residue_atoms['ARG']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atHD1', 
'atHD2', 
'atNE',
'atHE',
'atCZ',
'atNH1',  
'atHH11', 
'atHH12', 
'atNH2',  
'atHH21', 
'atHH22', 
'atC',
'atO'
]

covalent_bonds['ARG']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atHG1'  ], 
['atCG'  ,'atHG2'  ], 
['atCG'  ,'atCD'   ], 
['atCD'  ,'atHD1'  ], 
['atCD'  ,'atHD2'  ], 
['atCD'  ,'atNE'   ], 
['atNE'  ,'atHE'   ], 
['atNE'  ,'atCZ'   ], 
['atCZ'  ,'atNH1'  ], 
['atCZ'  ,'atNH2'  ], 
['atNH1' ,'atHH11' ], 
['atNH1' ,'atHH12' ], 
['atNH2' ,'atHH21' ], 
['atNH2' ,'atHH22' ], 
['atC'   ,'atO'    ] 
]

### ASN:

residue_atoms['ASN']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atOD1',  
'atND2',  
'atHD21', 
'atHD22', 
'atC',
'atO'
]

covalent_bonds['ASN']=[
['atN'   ,'atH'     ], 
['atN'   ,'atCA'    ], 
['atCA'  ,'atHA'    ], 
['atCA'  ,'atCB'    ], 
['atCA'  ,'atC'     ], 
['atCB'  ,'atHB1'   ], 
['atCB'  ,'atHB2'   ], 
['atCB'  ,'atCG'    ], 
['atCG'  ,'atOD1'   ], 
['atCG'  ,'atND2'   ], 
['atND2' ,'atHD21'  ], 
['atND2' ,'atHD22'  ], 
['atC'   ,'atO'     ] 
]

### ASP and ASPH [HD2]:

residue_atoms['ASP']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atOD1', 
'atOD2',
'atHD2',
'atC',
'atO',
]

covalent_bonds['ASP']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atOD1'  ], 
['atCG'  ,'atOD2'  ],
['atOD2' ,'atHD2'  ],
['atC'   ,'atO'    ] 
]

### CYS and CYSH [HG1]:

residue_atoms['CYS']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atSG',
'atHG1',
'atC',
'atO'
]
 
covalent_bonds['CYS']=[
['atN'   ,'atH'    ],
['atN'   ,'atCA'   ],
['atCA'  ,'atHA'   ],
['atCA'  ,'atCB'   ],
['atCA'  ,'atC'    ],
['atCB'  ,'atHB1'  ],
['atCB'  ,'atHB2'  ],
['atCB'  ,'atSG'   ],
['atSG'  ,'atHG1'  ],
['atC'   ,'atO'    ] 
]

### GLU:

residue_atoms['GLU']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atOE1', 
'atOE2', 
'atC',
'atO'
]

covalent_bonds['GLU']=[
['atN'   , 'atH'   ], 
['atN'   , 'atCA'  ], 
['atCA'  , 'atHA'  ], 
['atCA'  , 'atCB'  ], 
['atCA'  , 'atC'   ], 
['atCB'  , 'atHB1' ], 
['atCB'  , 'atHB2' ], 
['atCB'  , 'atCG'  ], 
['atCG'  , 'atHG1' ], 
['atCG'  , 'atHG2' ], 
['atCG'  , 'atCD'  ], 
['atCD'  , 'atOE1' ], 
['atCD'  , 'atOE2' ], 
['atC'   , 'atO'   ] 
]

### GLN:

residue_atoms['GLN']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atOE1', 
'atNE2', 
'atHE21',
'atHE22',
'atC',
'atO'
]

covalent_bonds['GLN']=[
['atN'   , 'atH'   ], 
['atN'   , 'atCA'  ], 
['atCA'  , 'atHA'  ], 
['atCA'  , 'atCB'  ], 
['atCA'  , 'atC'   ], 
['atCB'  , 'atHB1' ], 
['atCB'  , 'atHB2' ], 
['atCB'  , 'atCG'  ], 
['atCG'  , 'atHG1' ], 
['atCG'  , 'atHG2' ], 
['atCG'  , 'atCD'  ], 
['atCD'  , 'atOE1' ], 
['atCD'  , 'atNE2' ], 
['atNE2' , 'atHE21'], 
['atNE2' , 'atHE22'], 
['atC'   , 'atO'   ] 
]

### GLY:

residue_atoms['GLY']=[
'atN',
'atH',
'atCA',
'atHA1',
'atHA2',
'atC',
'atO'
]

covalent_bonds['GLY']=[
['atN'   , 'atH'   ],
['atN'   , 'atCA'  ],
['atCA'  , 'atHA1' ],
['atCA'  , 'atHA2' ],
['atCA'  , 'atC'   ],
['atC'   , 'atO'   ] 
]

### HISE (ND1 no H, NE2 with H),
### HISD (ND1 with H, NE2 no H),
### HISH (ND1 with H, NE2 with H),
### All included in HIS:

residue_atoms['HIS']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2',
'atCG',
'atND1', 
'atHD1',
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atNE2', 
'atHE2', 
'atC',
'atO'
]

covalent_bonds['HIS']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atND1'  ], 
['atCG'  ,'atCD2'  ], 
['atND1' ,'atHD1'  ], 
['atND1' ,'atCE1'  ], 
['atCD2' ,'atHD2'  ], 
['atCD2' ,'atNE2'  ], 
['atCE1' ,'atHE1'  ], 
['atCE1' ,'atNE2'  ], 
['atNE2' ,'atHE2'  ], 
['atC'   ,'atO'    ] 
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
'atN', 
'atH', 
'atCA', 
'atHA', 
'atCB', 
'atHB', 
'atCG1', 
'atHG11', 
'atHG12', 
'atCG2', 
'atHG21', 
'atHG22', 
'atHG23', 
'atCD', 
'atHD1', 
'atHD2', 
'atHD3', 
'atC', 
'atO' 
]

covalent_bonds['ILE']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB'   ], 
['atCB',  'atCG1'  ], 
['atCB',  'atCG2'  ], 
['atCG1', 'atHG11' ], 
['atCG1', 'atHG12' ], 
['atCG1', 'atCD'   ], 
['atCG2', 'atHG21' ], 
['atCG2', 'atHG22' ], 
['atCG2', 'atHG23' ], 
['atCD',  'atHD1'  ], 
['atCD',  'atHD2'  ], 
['atCD',  'atHD3'  ], 
['atC',   'atO'    ] 
]

### LEU:

residue_atoms['LEU']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atCG',
'atHG1',
'atCD',
'atHD11', 
'atHD12', 
'atHD13', 
'atCD2',
'atHD21', 
'atHD22', 
'atHD23', 
'atC',
'atO'
]

covalent_bonds['LEU']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'   ], 
['atCG',  'atCD'  ], 
['atCG',  'atCD2'  ], 
['atCD', 'atHD11' ], 
['atCD', 'atHD12' ], 
['atCD', 'atHD13' ], 
['atCD2', 'atHD21' ], 
['atCD2', 'atHD22' ], 
['atCD2', 'atHD23' ], 
['atC',   'atO'    ] 
]

### LYS with LYSH (HZ3):

residue_atoms['LYS']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atHG1',
'atHG2',
'atCD', 
'atHD1',
'atHD2',
'atCE', 
'atHE1',
'atHE2',
'atNZ', 
'atHZ1',
'atHZ2',
'atHZ3',
'atC',
'atO'
]

covalent_bonds['LYS']=[
['atN',   'atH'   ], 
['atN',   'atCA'  ], 
['atCA',  'atHA'  ], 
['atCA',  'atCB'  ], 
['atCA',  'atC'   ], 
['atCB',  'atHB1' ], 
['atCB',  'atHB2' ], 
['atCB',  'atCG'  ], 
['atCG',  'atHG1' ], 
['atCG',  'atHG2' ], 
['atCG',  'atCD'  ], 
['atCD',  'atHD1' ], 
['atCD',  'atHD2' ], 
['atCD',  'atCE'  ], 
['atCE',  'atHE1' ], 
['atCE',  'atHE2' ], 
['atCE',  'atNZ'  ], 
['atNZ',  'atHZ1' ], 
['atNZ',  'atHZ2' ], 
['atNZ',  'atHZ3' ], 
['atC',   'atO'   ] 
]

### MET:

residue_atoms['MET']=[
'atN',
'atH', 
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atHG1',
'atHG2',
'atSD', 
'atCE', 
'atHE1',
'atHE2',
'atHE3',
'atC',
'atO'
]

covalent_bonds['MET']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'  ], 
['atCG',  'atHG2'  ], 
['atCG',  'atSD'   ], 
['atSD',  'atCE'   ], 
['atCE',  'atHE1'  ], 
['atCE',  'atHE2'  ], 
['atCE',  'atHE3'  ], 
['atC',   'atO'    ] 
]

### PHE:

residue_atoms['PHE']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atCD', 
'atHD1', 
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atCE2', 
'atHE2', 
'atCZ',
'atHZ',
'atC',
'atO'
]

covalent_bonds['PHE']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atCD'  ], 
['atCG',  'atCD2'  ], 
['atCD', 'atHD1'  ], 
['atCD', 'atCE1'  ], 
['atCD2', 'atHD2'  ], 
['atCD2', 'atCE2'  ], 
['atCE1', 'atHE1'  ], 
['atCE1', 'atCZ'   ], 
['atCE2', 'atHE2'  ], 
['atCE2', 'atCZ'   ], 
['atCZ',  'atHZ'   ], 
['atC',   'atO'    ] 
]

### PRO:

residue_atoms['PRO']=[
'atN',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atHG1',
'atHG2',
'atCD', 
'atHD1',
'atHD2',
'atC', 
'atO'
]

covalent_bonds['PRO']=[
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'  ], 
['atCG',  'atHG2'  ], 
['atCG',  'atCD'   ], 
['atCD',  'atHD1'  ], 
['atCD',  'atHD2'  ], 
['atCD',  'atN'    ], 
['atC',   'atO'    ] 
]

### SER:

residue_atoms['SER']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atOG', 
'atHG1', 
'atC',
'atO'
]

covalent_bonds['SER']=[
['atN',  'atH'    ], 
['atN',  'atCA'   ], 
['atCA', 'atHA'   ], 
['atCA', 'atCB'   ], 
['atCA', 'atC'    ], 
['atCB', 'atHB1'  ], 
['atCB', 'atHB2'  ], 
['atCB', 'atOG'   ], 
['atOG', 'atHG1'  ], 
['atC',  'atO'    ] 
]

### THR:

residue_atoms['THR']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB',
'atOG1',
'atHG1',
'atCG2',
'atHG21',
'atHG22', 
'atHG23', 
'atC',
'atO'
]

covalent_bonds['THR']=[
['atN',   'atH'     ], 
['atN',   'atCA'    ], 
['atCA',  'atHA'    ], 
['atCA',  'atCB'    ], 
['atCA',  'atC'     ], 
['atCB',  'atHB'    ], 
['atCB',  'atOG1'   ], 
['atCB',  'atCG2'   ], 
['atOG1', 'atHG1'   ], 
['atCG2', 'atHG21'  ], 
['atCG2', 'atHG22'  ], 
['atCG2', 'atHG23'  ], 
['atC',   'atO'     ] 
]

### TRP:

residue_atoms['TRP']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atCD', 
'atHD1', 
'atCD2', 
'atNE1', 
'atHE1', 
'atCE2', 
'atCE3', 
'atHE3', 
'atCZ2', 
'atHZ2', 
'atCZ3', 
'atHZ3', 
'atCH2', 
'atHH2', 
'atC',
'atO'
]

covalent_bonds['TRP']=[
['atN',   'atH'   ], 
['atN',   'atCA'  ], 
['atCA',  'atHA'  ], 
['atCA',  'atCB'  ], 
['atCA',  'atC'   ], 
['atCB',  'atHB1' ], 
['atCB',  'atHB2' ], 
['atCB',  'atCG'  ], 
['atCG',  'atCD' ], 
['atCG',  'atCD2' ], 
['atCD', 'atHD1' ], 
['atCD', 'atNE1' ], 
['atCD2', 'atCE2' ], 
['atCD2', 'atCE3' ], 
['atNE1', 'atHE1' ], 
['atNE1', 'atCE2' ], 
['atCE2', 'atCZ2' ], 
['atCE3', 'atHE3' ], 
['atCE3', 'atCZ3' ], 
['atCZ2', 'atHZ2' ], 
['atCZ2', 'atCH2' ], 
['atCZ3', 'atHZ3' ], 
['atCZ3', 'atCH2' ], 
['atCH2', 'atHH2' ], 
['atC',   'atO'   ] 
]

### TYR:

residue_atoms['TYR']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atCD',
'atHD1', 
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atCE2', 
'atHE2', 
'atCZ', 
'atOH', 
'atHH', 
'atC',
'atO'
]

covalent_bonds['TYR']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atCD'  ], 
['atCG',  'atCD2'  ], 
['atCD', 'atHD1'  ], 
['atCD', 'atCE1'  ], 
['atCD2', 'atHD2'  ], 
['atCD2', 'atCE2'  ], 
['atCE1', 'atHE1'  ], 
['atCE1', 'atCZ'   ], 
['atCE2', 'atHE2'  ], 
['atCE2', 'atCZ'   ], 
['atCZ',  'atOH'   ], 
['atOH',  'atHH'   ], 
['atC',   'atO'    ] 
]

### VAL:

residue_atoms['VAL']=[
'atN',
'atH', 
'atCA',
'atHA',
'atCB',
'atHB',
'atCG1', 
'atHG11',
'atHG12',
'atHG13',
'atCG2', 
'atHG21',
'atHG22',
'atHG23',
'atC',
'atO'
]

covalent_bonds['VAL']=[
['atN',   'atH'     ], 
['atN',   'atCA'    ], 
['atCA',  'atHA'    ], 
['atCA',  'atCB'    ], 
['atCA',  'atC'     ], 
['atCB',  'atHB'    ], 
['atCB',  'atCG1'   ], 
['atCB',  'atCG2'   ], 
['atCG1', 'atHG11'  ], 
['atCG1', 'atHG12'  ], 
['atCG1', 'atHG13'  ], 
['atCG2', 'atHG21'  ], 
['atCG2', 'atHG22'  ], 
['atCG2', 'atHG23'  ], 
['atC',   'atO'     ] 
]

######## Terminals

### ACE:

residue_atoms['ACE']=[
'atCH3',
'atHH31',
'atHH32',
'atHH33',
'atC',
'atO'
]

covalent_bonds['ACE']=[
['atCH3' , 'atHH31' ],
['atCH3' , 'atHH32' ],
['atCH3' , 'atHH33' ],
['atCH3' , 'atC'    ],
['atC'   , 'atO'    ] 
]

### NME:

residue_atoms['NME']=[
'atN',
'atH',
'atCH3',
'atHH31',
'atHH32',
'atHH33'
]

covalent_bonds['NME']=[
['atN'   ,'atH'    ],
['atN'   ,'atCH3'  ],
['atCH3' ,'atHH31' ],
['atCH3' ,'atHH32' ],
['atCH3' ,'atHH33' ] 
]

### NHE:

residue_atoms['NHE']=[
'atN',
'atH1',
'atH2'
]

covalent_bonds['NHE']=[
['atN'   ,'atH1'   ],
['atN'   ,'atH2'   ]
]

## Terminals

terminal_bonds['atO']={}
terminal_bonds['atH']={}
terminal_bonds['none']={}

terminal_bonds['atO']['atOT1']  = 'atC'
terminal_bonds['atO']['atOT2']  = 'atC'

terminal_bonds['none']['atOXT']  = 'atC'

terminal_bonds['none']['atNT']  = 'atC'
terminal_bonds['none']['atHT1']  = 'atNT'
terminal_bonds['none']['atHT2']  = 'atNT'

terminal_bonds['atH']['atH1'] = 'atN'
terminal_bonds['atH']['atH2'] = 'atN'
terminal_bonds['atH']['atH3'] = 'atN'


#### NHE for ACEMD
terminal_bonds['atH']['atHT1'] = 'atN'
terminal_bonds['atH']['atHT2'] = 'atN'
terminal_bonds['atH']['atHT3'] = 'atN'

#### ACE for ACEMD
terminal_bonds['none']['atNT']   = 'atC'      
terminal_bonds['none']['atHY1']  = 'atCAY'
terminal_bonds['none']['atHY2']  = 'atCAY'
terminal_bonds['none']['atHY3']  = 'atCAY'
terminal_bonds['none']['atCAY']  = 'atCY'
terminal_bonds['none']['atOY']   = 'atCY'
terminal_bonds['none']['atCY']   = 'atN'

#### NME for ACEMD
terminal_bonds['none']['atNT']   = 'atC'      # 
terminal_bonds['none']['atHT1']  = 'atCAT'    # 
terminal_bonds['none']['atHT2']  = 'atCAT'    # 
terminal_bonds['none']['atHT3']  = 'atCAT'    # 
terminal_bonds['none']['atCAT']  = 'atNT'     # 
terminal_bonds['none']['atHT']   = 'atNT'      # 



##### LIPIDS:

### AOT:

residue_atoms['AOT']=[
'atS',
'atOS1',
'atOS2',
'atOS3',
'atC1',
'atH1',
'atC2',
'atH2',
'atH3',
'atC3',
'atOT1',
'atOT2',
'atC12',
'atO3',
'atO4',
'atC13',
'atH21',
'atH22',
'atC14',
'atH23',
'atC19',
'atH33',
'atH34',
'atC20',
'atH35',
'atH36',
'atH37',
'atC15',
'atH24',
'atH25',
'atC16',
'atH26',
'atH27',
'atC17',
'atH28',
'atH29',
'atC18',
'atH30',
'atH31',
'atH32',
'atC4',
'atH4',
'atH5',
'atC4',
'atC5',
'atC6',
'atC7',
'atC8',
'atC9',
'atC10',
'atC11',
'atH6',
'atH7',
'atH8',
'atH9',
'atH10',
'atH11',
'atH12',
'atH13',
'atH14',
'atH15',
'atH16',
'atH17',
'atH18',
'atH19',
'atH20'
]

covalent_bonds['AOT']=[
['atS'   ,'atOS1'   ],
['atS'   ,'atOS2'   ],
['atS'   ,'atOS3'   ],
['atS'   ,'atC1'    ],
['atC1'  ,'atH1'    ],
['atC1'  ,'atC2'    ],
['atC2'  ,'atH2'    ],
['atC2'  ,'atH3'    ],
['atC2'  ,'atC12'   ],
['atC12' ,'atO3'    ],
['atC12' ,'atO4'    ],
['atO4'  ,'atC13'   ],
['atC13' ,'atH21'   ],
['atC13' ,'atH22'   ],
['atC13' ,'atC14'   ],
['atC14' ,'atH23'   ],
['atC14' ,'atC19'   ],
['atC19' ,'atH33'   ],
['atC19' ,'atH34'   ],
['atC19' ,'atC20'   ],
['atC20' ,'atH35'   ],
['atC20' ,'atH36'   ],
['atC20' ,'atH37'   ],
['atC14' ,'atC15'   ],
['atC15' ,'atH24'   ],
['atC15' ,'atH25'   ],
['atC15' ,'atC16'   ],
['atC16' ,'atH26'   ],
['atC16' ,'atH27'   ],
['atC16' ,'atC17'   ],
['atC17' ,'atH28'   ],
['atC17' ,'atH29'   ],
['atC17' ,'atC18'   ],
['atC18' ,'atH30'   ],
['atC18' ,'atH31'   ],
['atC18' ,'atH32'   ],
['atC1'  ,'atC3'    ],
['atC3'  ,'atOT1'   ],
['atC3'  ,'atOT2'   ],
['atOT2' ,'atC4'    ],
['atC4'  ,'atH4'    ],
['atC4'  ,'atH5'    ],
['atC4'  ,'atC5'    ],
['atC5'  ,'atH6'    ],
['atC5'  ,'atC10'   ],
['atC10' ,'atH16'   ],
['atC10' ,'atH17'   ],
['atC10' ,'atC11'   ],
['atC11' ,'atH18'   ],
['atC11' ,'atH19'   ],
['atC11' ,'atH20'   ],
['atC5'  ,'atC6'    ],
['atC6'  ,'atH7'    ],
['atC6'  ,'atH8'    ],
['atC6'  ,'atC7'    ],
['atC7'  ,'atH9'    ],
['atC7'  ,'atH10'   ],
['atC7'  ,'atC8'    ],
['atC8'  ,'atH11'   ],
['atC8'  ,'atH12'   ],
['atC8'  ,'atC9'    ],
['atC9'  ,'atH13'   ],
['atC9'  ,'atH14'   ],
['atC9'  ,'atH15'   ]
]

### POPC:

residue_atoms['POPC']=[
'atN',
'atC12',  
'atH12A',
'atH12B', 
'atC13',
'atH13A',
'atH13B',
'atH13C',
'atC14',
'atH14A', 
'atH14B', 
'atH14C', 
'atC15' , 
'atH15A', 
'atH15B', 
'atH15C', 
'atC11',  
'atH11A', 
'atH11B', 
'atP',    
'atO13',  
'atO14',  
'atO12',  
'atO11',  
'atC1',   
'atHA',   
'atHB',   
'atC2',   
'atHS',   
'atO21',  
'atC21',  
'atO22',  
'atC22',  
'atH2R',  
'atH2S',  
'atC3',  
'atHX',  
'atHY',  
'atO31',  
'atC31',  
'atO32',  
'atC32',  
'atH2X',  
'atH2Y',  
'atC23',  
'atH3R',  
'atH3S',  
'atC24',  
'atH4R',  
'atH4S',  
'atC25',  
'atH5R',  
'atH5S',  
'atC26',
'atH6R', 
'atH6S', 
'atC27', 
'atH7R', 
'atH7S', 
'atC28', 
'atH8R', 
'atH8S', 
'atC29', 
'atH91', 
'atC210', 
'atH101', 
'atC211', 
'atH11R', 
'atH11S', 
'atC212', 
'atH12R', 
'atH12S', 
'atC213', 
'atH13R', 
'atH13S', 
'atC214', 
'atH14R', 
'atH14S', 
'atC215', 
'atH15R', 
'atH15S', 
'atC216', 
'atH16R', 
'atH16S', 
'atC217', 
'atH17R', 
'atH17S', 
'atC218', 
'atH18R', 
'atH18S', 
'atH18T', 
'atC33', 
'atH3X', 
'atH3Y', 
'atC34', 
'atH4X', 
'atH4Y', 
'atC35', 
'atH5X', 
'atH5Y', 
'atC36', 
'atH6X', 
'atH6Y', 
'atC37', 
'atH7X', 
'atH7Y', 
'atC38', 
'atH8X', 
'atH8Y',
'atC39',
'atH9X',   
'atH9Y',   
'atC310',  
'atH10X',  
'atH10Y',  
'atC311',  
'atH11X',  
'atH11Y',  
'atC312',  
'atH12X',  
'atH12Y',  
'atC313',  
'atH13X',  
'atH13Y',  
'atC314',  
'atH14X',  
'atH14Y',  
'atC315',  
'atH15X',  
'atH15Y',  
'atC316',  
'atH16X',  
'atH16Y',  
'atH16Z'  
]


covalent_bonds['POPC']=[
['atN',    'atC13' ], 
['atN',    'atC14' ], 
['atN',    'atC15' ], 
['atC13',  'atH13A'], 
['atC13',  'atH13B'], 
['atC13',  'atH13C'], 
['atC14',  'atH14A'], 
['atC14',  'atH14B'], 
['atC14',  'atH14C'], 
['atC15',  'atH15A'], 
['atC15',  'atH15B'], 
['atC15',  'atH15C'], 
['atN',    'atC12' ], 
['atC12',  'atH12A'], 
['atC12',  'atH12B'], 
['atC12',  'atC11' ], 
['atC11',  'atH11A'], 
['atC11',  'atH11B'], 
['atC11',  'atO12' ], 
['atO11',  'atC1'  ], 
['atO12',  'atP'   ], 
['atP',    'atO11' ], 
['atP',    'atO13' ], 
['atP',    'atO14' ], 
['atC1',   'atHA'  ], 
['atC1',   'atHB'  ], 
['atC1',   'atC2'  ], 
['atC2',   'atHS'  ], 
['atC2',   'atC3'  ], 
['atC2',   'atO21' ], 
['atC3',   'atHX'  ], 
['atC3',   'atHY'  ], 
['atC3',   'atO31' ], 
['atO21',  'atC21' ], 
['atC21',  'atC22' ], 
['atC22',  'atH2R' ], 
['atC22',  'atH2S' ], 
['atC22',  'atC23' ], 
['atC23',  'atH3R' ], 
['atC23',  'atH3S' ], 
['atC23',  'atC24' ], 
['atC24',  'atH4R' ], 
['atC24',  'atH4S' ], 
['atC24',  'atC25' ], 
['atC25',  'atH5R' ], 
['atC25',  'atH5S' ], 
['atC25',  'atC26' ], 
['atC26',  'atH6R' ], 
['atC26',  'atH6S' ], 
['atC26',  'atC27' ], 
['atC27',  'atH7R' ], 
['atC27',  'atH7S' ], 
['atC27',  'atC28' ], 
['atC28',  'atH8R' ], 
['atC28',  'atH8S' ], 
['atC28',  'atC29' ], 
['atC29',  'atH91' ], 
['atC210', 'atH101'], 
['atC210', 'atC211'], 
['atC211', 'atH11R'], 
['atC211', 'atH11S'], 
['atC211', 'atC212'], 
['atC212', 'atH12R'], 
['atC212', 'atH12S'], 
['atC212', 'atC213'], 
['atC213', 'atH13R'], 
['atC213', 'atH13S'], 
['atC213', 'atC214'], 
['atC214', 'atH14R'], 
['atC214', 'atH14S'], 
['atC214', 'atC215'], 
['atC215', 'atH15R'], 
['atC215', 'atH15S'], 
['atC215', 'atC216'], 
['atC216', 'atH16R'], 
['atC216', 'atH16S'], 
['atC216', 'atC217'], 
['atC217', 'atH17R'], 
['atC217', 'atH17S'], 
['atC217', 'atC218'], 
['atC218', 'atH18R'], 
['atC218', 'atH18S'], 
['atC218', 'atH18T'], 
['atO31',  'atC31' ], 
['atC31',  'atC32' ], 
['atC32',  'atH2X' ], 
['atC32',  'atH2Y' ], 
['atC32',  'atC33' ], 
['atC33',  'atH3X' ], 
['atC33',  'atH3Y' ], 
['atC33',  'atC34' ], 
['atC34',  'atH4X' ], 
['atC34',  'atH4Y' ], 
['atC34',  'atC35' ], 
['atC35',  'atH5X' ], 
['atC35',  'atH5Y' ], 
['atC35',  'atC36' ], 
['atC36',  'atH6X' ], 
['atC36',  'atH6Y' ], 
['atC36',  'atC37' ], 
['atC37',  'atH7X' ], 
['atC37',  'atH7Y' ], 
['atC37',  'atC38' ], 
['atC38',  'atH8X' ], 
['atC38',  'atH8Y' ], 
['atC38',  'atC39' ], 
['atC39',  'atH9X' ], 
['atC39',  'atH9Y' ], 
['atC39',  'atC310'], 
['atC310', 'atH10X'], 
['atC310', 'atH10Y'], 
['atC310', 'atC311'], 
['atC311', 'atH11X'], 
['atC311', 'atH11Y'], 
['atC311', 'atC312'], 
['atC312', 'atH12X'], 
['atC312', 'atH12Y'], 
['atC312', 'atC313'], 
['atC313', 'atH13X'], 
['atC313', 'atH13Y'], 
['atC313', 'atC314'], 
['atC314', 'atH14X'], 
['atC314', 'atH14Y'], 
['atC314', 'atC315'], 
['atC315', 'atH15X'], 
['atC315', 'atH15Y'], 
['atC315', 'atC316'], 
['atC316', 'atH16X'], 
['atC316', 'atH16Y'], 
['atC316', 'atH16Z'], 
['atC21',  'atO22' ], 
['atC29',  'atC210'], 
['atC31',  'atO32' ] 
]

### CHL1:

residue_atoms['CHL1']=[
'atC3'  , 
'atO3'  , 
"atH3'" , 
'atH3'  , 
'atC4'  , 
'atH4A' , 
'atH4B' , 
'atC5'  , 
'atC6'  , 
'atH6'  , 
'atC7'  , 
'atH7A' , 
'atH7B' , 
'atC8'  , 
'atH8'  , 
'atC14' , 
'atH14' , 
'atC15' , 
'atH15A', 
'atH15B', 
'atC16' , 
'atH16A', 
'atH16B', 
'atC17' , 
'atH17' , 
'atC13' , 
'atC18' , 
'atH18A', 
'atH18B', 
'atH18C', 
'atC12' , 
'atH12A', 
'atH12B', 
'atC11' , 
'atH11A', 
'atH11B', 
'atC9'  , 
'atH9'  , 
'atC10' , 
'atC19' , 
'atH19A', 
'atH19B', 
'atH19C', 
'atC1'  , 
'atH1A' , 
'atH1B' , 
'atC2'  , 
'atH2A' , 
'atH2B' , 
'atC20' , 
'atH20' , 
'atC21' , 
'atH21A', 
'atH21B', 
'atH21C', 
'atC22' , 
'atH22A', 
'atH22B', 
'atC23' , 
'atH23A', 
'atH23B', 
'atC24' , 
'atH24A', 
'atH24B', 
'atC25' , 
'atH25' , 
'atC26' , 
'atH26A', 
'atH26B', 
'atH26C', 
'atC27' , 
'atH27A', 
'atH27B', 
'atH27C' 
]


covalent_bonds['CHL1']=[
['atC3'  ,  'atO3'   ], 
['atC3'  ,  'atH3'   ], 
['atO3'  ,  "atH3'"  ], 
['atC3'  ,  'atC2'   ], 
['atC2'  ,  'atH2A'  ], 
['atC2'  ,  'atH2B'  ], 
['atC2'  ,  'atC1'   ], 
['atC1'  ,  'atH1A'  ], 
['atC1'  ,  'atH1B'  ], 
['atC3'  ,  'atC4'   ], 
['atC4'  ,  'atH4A'  ], 
['atC4'  ,  'atH4B'  ], 
['atC4'  ,  'atC5'   ], 
['atC5'  ,  'atC10'  ], 
['atC10' ,  'atC1'   ], 
['atC10' ,  'atC19'  ], 
['atC19' ,  'atH19A' ], 
['atC19' ,  'atH19B' ], 
['atC19' ,  'atH19C' ], 
['atC6'  ,  'atH6'   ], 
['atC6'  ,  'atC7'   ], 
['atC7'  ,  'atH7A'  ], 
['atC7'  ,  'atH7B'  ], 
['atC7'  ,  'atC8'   ], 
['atC8'  ,  'atH8'   ], 
['atC8'  ,  'atC9'   ], 
['atC9'  ,  'atH9'   ], 
['atC9'  ,  'atC10'  ], 
['atC8'  ,  'atC14'  ], 
['atC14' ,  'atH14'  ], 
['atC14' ,  'atC13'  ], 
['atC13' ,  'atC12'  ], 
['atC12' ,  'atH12A' ], 
['atC12' ,  'atH12B' ], 
['atC12' ,  'atC11'  ], 
['atC11' ,  'atH11A' ], 
['atC11' ,  'atH11B' ], 
['atC11' ,  'atC9'   ], 
['atC13' ,  'atC18'  ], 
['atC18' ,  'atH18A' ], 
['atC18' ,  'atH18B' ], 
['atC18' ,  'atH18C' ], 
['atC14' ,  'atC15'  ], 
['atC15' ,  'atH15A' ], 
['atC15' ,  'atH15B' ], 
['atC15' ,  'atC16'  ], 
['atC16' ,  'atH16A' ], 
['atC16' ,  'atH16B' ], 
['atC16' ,  'atC17'  ], 
['atC17' ,  'atH17'  ], 
['atC17' ,  'atC13'  ], 
['atC17' ,  'atC20'  ], 
['atC20' ,  'atH20'  ], 
['atC20' ,  'atC21'  ], 
['atC21' ,  'atH21A' ], 
['atC21' ,  'atH21B' ], 
['atC21' ,  'atH21C' ], 
['atC20' ,  'atC22'  ], 
['atC22' ,  'atH22A' ], 
['atC22' ,  'atH22B' ], 
['atC22' ,  'atC23'  ], 
['atC23' ,  'atH23A' ], 
['atC23' ,  'atH23B' ], 
['atC23' ,  'atC24'  ], 
['atC24' ,  'atH24A' ], 
['atC24' ,  'atH24B' ], 
['atC24' ,  'atC25'  ], 
['atC25' ,  'atH25'  ], 
['atC25' ,  'atC26'  ], 
['atC26' ,  'atH26A' ], 
['atC26' ,  'atH26B' ], 
['atC26' ,  'atH26C' ], 
['atC25' ,  'atC27'  ], 
['atC27' ,  'atH27A' ], 
['atC27' ,  'atH27B' ], 
['atC27' ,  'atH27C' ], 
['atC5'  ,  'atC6'   ] 
]


##### WATER:

### SOL3:

residue_atoms['SOL3']=[
'atOW',
'atHW1',
'atHW2'
]

covalent_bonds['SOL3']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

### SOL4:

residue_atoms['SOL4']=[
'atOW',
'atHW1',
'atHW2',
'atvir'
]

covalent_bonds['SOL4']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

### SOL5:

residue_atoms['SOL5']=[
'atOW',
'atHW1',
'atHW2',
'atvir',
'atvir'
]

covalent_bonds['SOL5']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]


##### IONS:

## NA:
residue_atoms['NA']=[
'itNA'
]

covalent_bonds['NA']=[
]

## K:
residue_atoms['K']=[
'itK'
]

covalent_bonds['K']=[
]

## LI:
residue_atoms['LI']=[
'itLI'
]

covalent_bonds['LI']=[
]

## CL:
residue_atoms['CL']=[
'itCL'
]

covalent_bonds['CL']=[
]



