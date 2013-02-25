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

### ASP:

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
'atO'
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
terminal_bonds['none']['atHT']  = 'atNT'      # 



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


