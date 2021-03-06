####################################
# GENERAL COMMENTS

NAME_VERSION="Aqua 0.1"

#
# This module requires:
#

### External libraries:
from os import system
from os import path
from os import sys
import copy as ccopy
import numpy as numpy

try:
    import cPickle as pic
except:
    import pickle as pic
import datetime as datetime

### pyno libraries:
from cl_coors import *
import top_par as tp_aa
import top_par_cg as tp_cg
from libgeneral import glob as faux
#from libmss import glob as mss_funcs
import libmath as libmath

### topologies all atom by default
### and top added by user
tp=tp_aa
for ii in tp.user_topol:
    tp.add_topol(tp,ii)

#
# Structure of the file:
#
# -Class Set
#            - Instantiation
#            - Functions
# -Class Unit
# -Class Residue
# -Class Water
# -External functions
#

#
# END GENERAL COMMENTS
####################################


#######################################################
#######################################################
#### CLASSES
    
####
#### Labels or common attributes to be inherited
####

class labels_unit():                           # Every unit (atom) has at least these attributes
    def __init__(self):
        self.name=None
        self.index=None
        self.pdb_index=None
        self.covalent_bonds=[]                 # List of atoms covalently bonded.
        self.hbonds=[]                         # Atom h-bonded: [[atom_index][strength]]
        self.type=None                         # for residues: protein,ion,solv. for atom: atom_nature (H,O,...)
        self.__int_name__=None                   # internal name of atom

class labels_set(labels_unit):                 # Every set of units (chain, residue, molecule) has at least these attributes
    def __init__(self):
        self.num_atoms=0
        self.list_atoms=[]

class labels_parent(labels_unit):               # A selection always keep memory of the previous location
    def __init__(self,parent,argument=None,unit=False):
        if unit:
            self.name=parent.name
            self.index=parent.index
            self.pdb_index=parent.pdb_index
            self.__int_name__=parent.__int_name__
        else:
            self.name=parent.name
            self.condition=argument

####
#### Class Unit (atom)
####


class cl_unit(labels_unit):                     # Attributes of an atom

    def __init__(self):
        '''Initialize an atom object'''

        # From labels_unit: .name, .index, .pdb_index, .covalent_bonds    

        # > Topological properties

        self.resid=labels_unit()        # Residue which this atom belongs to. 
        self.chain=labels_unit()        # Chain which this atom belongs to.
        self.parent=None                # Selection which this atom comes from.

        self.alt_loc=0                  # Alternate location (PDB)
        self.code_ins_res=0             # Code of insertion of residues (PDB)
        self.seg_ident=''               # Index segment (PDB)
        self.elem_symb=''               # Element symbol (PDB)
        self.type_pdb=''                # Type of atom for the PDB (ATOM,HETATM)

        self.covalent_bonds=[]          # esto deberia estar heredado de labels_unit @!!!!!!@

        # > Physical and Chemical properties

        self.coors=[]
        self.mass=0.0                   # Mass
        self.charge=0.0                 # Charge
        self.vdw=0.0                    # VdW radius
        self.occup=0.0                  # Occupation (PDB)
        self.bfactor=0.0                # B-Factor
        self.acceptor=False             # True or false 
        self.donor=False                # True or false
        self.polarizability=False       # True of false

    def info(self,resid=True,pdb=True):
        if resid:
            if pdb:
                return self.name+'-'+str(self.pdb_index)+'/'+self.resid.name+'-'+str(self.resid.pdb_index)
            else:
                return self.name+'-'+str(self.index)+'/'+self.resid.name+'-'+str(self.resid.index)
        else:
            if pdb:
                return self.name+'-'+str(self.pdb_index)
            else:
                return self.name+'-'+str(self.index)

####
#### Class residue (set of atoms)
####


class cl_residue(labels_set):           # Attributes of a residue (inherits from label_set)

    def __init__( self ):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms, .__int_type

        self.chain=labels_unit()         # Chain which this residue belongs to.
        self.__int_dict_atoms__={}         # Dictionary with internal atom names and indexes

        pass

    def info(self,pdb=True):
        if pdb:
            return self.name+'-'+str(self.pdb_index)
        else:
            return self.name+'-'+str(self.index)


####
#### Class water (set of atoms)
####

class cl_water(labels_set):             # Attributes of a water molecule
    '''Water specific attributes and methods'''

    def __init__( self ,water_model=None):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        self.model=water_model
                
        self.O=labels_unit()
        self.H1=labels_unit()
        self.H2=labels_unit()
        self.uvect_norm=[]
        self.bisector=[]
        self.microstate=''
        
        pass



####
#### Class set (set of atoms: Molecule)
####


class msystem(labels_set):               # The supra-estructure: System (waters+cofactors+proteins...)

    
    def __init__(self,input_file=None,download=None,coors=False,with_bonds=True,missing_atoms=True,wrap=True,cg=False,verbose=False):


        # Load topology coarse-grained if this is the case
        if cg:
            tp=tp_cg
            self.ff='cg'
        else:
            tp=tp_aa
            self.ff='aa'

        # If download:

        if download:
            if not download.endswith('.pdb'):
                    download=download+'.pdb'
            input_file=download
            if not path.exists(input_file):
                temp='wget -nv http://www.rcsb.org/pdb/files/'+input_file+' 1>/dev/null 2>&1'
                system(temp)
                if path.exists(input_file):
                    print '# File saved as '+input_file
                else:
                    print '# Error downloading http://www.rcsb.org/pdb/files/'+input_file

            else:
                print '# The file '+input_file+' exists in the local folder. Loading it...'

            self.file_topol=input_file


        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        # > Instantation options:
        self.name=None
        self.file_topol=input_file       # input file name
        self.file_topol_type='' 
        if self.file_topol:
            self.file_topol_type=input_file.split('.')[-1] # pdb,gro,psf
        self.file_hbonds=''             # still not useful -do not remove-
        self.file_mss=''                # still not useful -do not remove-
        self.file_shell=''              # still not useful -do not remove-
        self.file_ff=''                 # all-atom (aa) or coarse-grained (cg)
        #self.selection=False           # is the set a selection (a mirror -pointer- of a parent set)
        
        # > Topological properties
        self.atom=[]                    # list of atoms    (objects: cl_unit)
        self.resid=[]                   # list of residues (objects: molecule)
        self.chain=[]                   # list of chains   (objects: molecule)
        self.chains=[]                  # list of chain names (strings)
        self.protein=[]
        self.lipid=[]
        self.ligand=[]
        self.molecule=[]
        self.ion=[]                     # list of ions (objects: molecule)
        self.water=[]                   # list of waters   (objects: cl_water)
        self.water_model=None           # water model
        self.parent=None                # Parent set if it comes from selection (object: labels_parent)

        # > Physical and Chemical properties
        self.acceptors=[]               # list of acceptors
        self.donors=[]             # list of [donors,donor_hydrogen]
        self.dimensionality=0           # dimensionality (num_atoms*3)

        # > Coordinates and Trajectory
        self.which_traj=None
        self.which_frame=None
        self.box=[]
        self.cell=[]
        self.traj=[]
        self.dist_matrix=[]             # distance matrix (This should go to cl_coors?)

        # > Info PDB
        self.pdb_header=[]              # PDB Header (HEAD + TITLE)
        self.pdb_ss=[]                  # PDB Secondary structure

        # > Internal info for internal procedures


        ##################################

        # A SYSTEM CAN BE BUILT FROM A FILE OR FROM A SELECTION




        # BUILDING THE SET FROM A FILE

        if input_file:

            # The file does not exist?
            if not download and not path.exists(input_file):
                print "# The file "+input_file+" does not exist."
                return

            # Reading the file and
            # attaching the atoms: self.atom[] (cl_unit)

            self.load_topol(self.file_topol)

            # Finnal set up of the attributes of the cl_units in self.atom[]:
            
            # Atoms and residues:
            without_hs=True
            before_resid=None
            before_chain=None
            jj=-1    # index residue for pynoramix
            ii=-1    # index atom for pynoramix
            kk=-1    # index chain for pynoramix
            for atom in self.atom:
                if atom.chain.name!=before_chain :   #### Chain
                    #if atom.type_pdb in ['ATOM']:
                    before_chain=atom.chain.name
                    kk+=1
                    self.chains.append(atom.chain.name)
                    temp_chain=labels_set()
                    temp_chain.name=atom.chain.name
                    temp_chain.index=kk
                    self.chain.append(temp_chain)
                if atom.resid.pdb_index!=before_resid :    #### Residue
                    before_resid=atom.resid.pdb_index
                    jj+=1
                    temp_residue=cl_residue()
                    temp_residue.index=jj
                    temp_residue.list_atoms=[]
                    temp_residue.pdb_index=atom.resid.pdb_index
                    temp_residue.name=atom.resid.name
                    temp_residue.chain.name=atom.chain.name
                    temp_residue.chain.index=kk
                    temp_residue.__int_dict_atoms__={}
                    try:
                        temp_residue.type=tp.residue_type[temp_residue.name]
                    except:
                        print temp_residue.name, 'without residue.type'
                    self.resid.append(temp_residue)
                ii+=1                                      #### Atom
                atom.index=ii                   
                atom.resid.index=jj
                try:
                    atom.resid.type=self.resid[jj].type
                except:
                    pass
                atom.chain.index=kk
                atom.hbonds=[]
                self.resid[jj].list_atoms.append(ii)
                self.chain[kk].list_atoms.append(ii)

            if self.file_topol_type not in ['psf']:
                for residue in self.resid:
                    residue.__int_name__=tp.residue[residue.name]
                for atom in self.atom:
                    jj=atom.resid.index
                    atom.resid.__int_name__=self.resid[jj].__int_name__
                    atom.__int_name__=tp.atom[atom.name]
                    atom.type=tp.atom_type[atom.__int_name__]
                    self.resid[jj].__int_dict_atoms__[atom.__int_name__]=atom.index
                    if atom.type=='H':
                        without_hs=False
            else:
                for residue in self.resid:
                    if residue.type=='Water':
                        residue.__int_name__=tp.residue[residue.name]
                        for ii in residue.list_atoms:
                            self.atom[ii].__int_name__=tp.atom[self.atom[ii].name]
                            residue.__int_dict_atoms__[self.atom[ii].__int_name__]=ii

            ### Setting up the subsets.
            num_wat=0
            for residue in self.resid[:]:
                
                if residue.type=='Water':       ### Waters
                    
                    if without_hs:
                        residue.__int_name__='SOL3'
                    else:
                        residue.__int_name__='SOL'+str(len(residue.list_atoms))
             
                    temp_water=cl_water()
                    temp_water.name=residue.name
                    temp_water.index=residue.index
                    temp_water.pdb_index=residue.pdb_index
                    temp_water.list_atoms=residue.list_atoms
                    temp_water.model=self.water_model
                    temp_water.__int_name__=residue.__int_name__
                    for aa in residue.list_atoms:
                        self.atom[aa].resid.__int_name__=residue.__int_name__
                        self.atom[aa].resid.water=num_wat
             
                    if 'atO' in residue.__int_dict_atoms__.keys():
                        xxx=residue.__int_dict_atoms__.pop('atO')
                        residue.__int_dict_atoms__['atOW']=xxx
                        self.atom[xxx].__int_name__='atOW'
                    if 'atH1' in residue.__int_dict_atoms__.keys():
                        xxx=residue.__int_dict_atoms__.pop('atH1')
                        residue.__int_dict_atoms__['atHW1']=xxx
                        self.atom[xxx].__int_name__='atHW1'
                    if 'atH2' in residue.__int_dict_atoms__.keys():
                        xxx=residue.__int_dict_atoms__.pop('atH2')
                        residue.__int_dict_atoms__['atHW2']=xxx
                        self.atom[xxx].__int_name__='atHW2'

                    for aa in [[False,'atOW',temp_water.O],[without_hs,'atHW1',temp_water.H1],[without_hs,'atHW2',temp_water.H2]]:
                        if not aa[0]:
                            xxx=residue.__int_dict_atoms__[aa[1]]
                            aa[2].index=xxx
                            aa[2].pdb_index=self.atom[xxx].pdb_index
                            aa[2].name=self.atom[xxx].name
             
                    self.water.append(temp_water)
                    residue.water=num_wat
                    num_wat+=1
             
                if residue.type=='Ion':        ### Ions
                    temp_residue=cl_residue()
                    temp_residue.list_atoms=residue.list_atoms
                    temp_residue.index=residue.index
                    temp_residue.pdb_index=residue.pdb_index
                    temp_residue.name=residue.name
                    self.ion.append(temp_residue)

            del(num_wat)
            ### Setting up the local attributes
            if self.file_topol_type not in ['psf']:

                # Topology and Covalent bonds

                 if with_bonds:
                  
                     for residue in self.resid[:]:
                  
                         # Found, missing or unknown atoms in the residue and its topology
                         jj=residue.index
                         tp_residue_name=residue.__int_name__
                         tp_residue_atoms=tp.residue_atoms[tp_residue_name]
                         found=residue.__int_dict_atoms__.keys()
                         if without_hs:
                             aux=ccopy.deepcopy(tp_residue_atoms)
                             for ii in tp_residue_atoms:
                                 if tp.atom_type[ii]=='H':
                                     aux.remove(ii)
                             tp_residue_atoms=aux
                         missing=list(set(tp_residue_atoms).difference(found))
                         unknown=list(set(found).difference(tp_residue_atoms))
                  
                         # Covalent bonds: Topology residue
                         for ii in tp.covalent_bonds[tp_residue_name]:
                             try:
                                 aa=residue.__int_dict_atoms__[ii[0]]
                                 bb=residue.__int_dict_atoms__[ii[1]]
                                 self.atom[aa].covalent_bonds.append(bb)
                                 self.atom[bb].covalent_bonds.append(aa)
                             except:
                                 pass
                  
                         # all atom Peptide bond [C(i)->N(i+1)] or Coarse Grained peptide bond [BB(i)->BB(i+1)]
                         try:
                             next_residue=self.resid[residue.index+1]
                             if residue.type=='Protein' and next_residue.type=='Protein':
                                 if residue.chain.name==next_residue.chain.name :
                                     aa=residue.__int_dict_atoms__[tp.peptidic_bond[0]]
                                     bb=next_residue.__int_dict_atoms__[tp.peptidic_bond[1]]
                                     self.atom[aa].covalent_bonds.append(bb)
                                     self.atom[bb].covalent_bonds.append(aa)
                         except:
                             pass
                  
                         # Covalent bonds: Terminals
                         unk2rm={}; miss2rm={}
                         if unknown:
                             for jj in (missing+['none']):
                                 for ii in unknown:
                                     try: 
                                         kk=tp.terminal_bonds[jj][ii]
                                         aa=residue.__int_dict_atoms__[ii]
                                         bb=residue.__int_dict_atoms__[kk]
                                         self.atom[aa].covalent_bonds.append(bb)
                                         self.atom[bb].covalent_bonds.append(aa)
                                         if jj!='none': 
                                             miss2rm[jj]=''
                                         unk2rm[ii]=''
                                     except:
                                         pass
                             for aa in miss2rm.keys():
                                 missing.remove(aa)
                             for aa in unk2rm.keys():
                                 unknown.remove(aa)
                  
                         # Listing missing or unknown atoms
                         if missing_atoms:
                             for ii in missing:
                                 print '# No atom type',ii,'in', residue.name, residue.pdb_index
                             for ii in unknown:
                                 print '# Unknown atom type',ii,'in', residue.name, residue.pdb_index

                 ## Aqui pongo las cadenas reales
                 cov_chains=self.selection_covalent_chains(chain='ALL',select='ALL')
                 for cov_chain in cov_chains:
                     if self.atom[cov_chain[0]].resid.type=='Lipid':
                         self.lipid.append(cov_chain)
                     elif self.atom[cov_chain[0]].resid.type=='Protein':
                         self.protein.append(cov_chain)
                     elif self.atom[cov_chain[0]].resid.type=='Ligand':
                         self.ligand.append(cov_chain)
                     elif self.atom[cov_chain[0]].resid.type=='Molecule':
                         self.molecule.append(cov_chain)
                         
                 del(cov_chains)

                 # Charge
                  
                 for atom in self.atom[:]:
                     if tp.atom[atom.name] in tp.charge:
                         atom.charge=tp.charge[tp.atom[atom.name]]
                  
                 # Acceptors-Donors
                  
                 for atom in self.atom[:]:
                     # Donors default
                     if atom.__int_name__ in tp.donors: 
                         atom.donor=True
                     # Donors exceptions
                     if atom.__int_name__ in tp.donors_exception:
                         try:
                             exception=tp.donors_exception[atom.__int_name__][atom.resid.__int_name__]
                             if exception[0]=='Always':
                                 atom.donor=exception[1]
                             elif exception[0]=='Without H':
                                 if 'H' not in [self.atom[ii].type for ii in atom.covalent_bonds]:
                                     atom.donor=exception[1]
                             else:
                                 print '# Error with donors exceptions:', atom.name, atom.resid.name
                         except:
                             pass
                     # Acceptors default
                     if atom.__int_name__ in tp.acceptors: 
                         atom.acceptor=True
                     # Acceptors exceptions
                     if atom.__int_name__ in tp.acceptors_exception:
                         try:
                             exception=tp.acceptors_exception[atom.__int_name__][atom.resid.__int_name__]
                             if exception[0]=='Always':
                                 atom.acceptor=exception[1]
                             elif exception[0]=='With H':
                                 if 'H' in [self.atom[ii].type for ii in atom.covalent_bonds]:
                                     atom.acceptor=exception[1]
                             else:
                                 print '# Error with acceptors exceptions:', atom.name, atom.resid.name
                         except:
                             pass
                     
                     if atom.acceptor:
                         self.acceptors.append(atom.index)
                  
                     if atom.donor:
                         for ii in atom.covalent_bonds:
                             if self.atom[ii].type=='H':
                                 self.donors.append([atom.index,ii])
                  
                 self.acceptors=numpy.array(self.acceptors,dtype=int,order='Fortran')
                 self.donors=numpy.array(self.donors,dtype=int,order='Fortran')

            ### Setting up the global attributes

            self.name=self.file_topol[:-self.file_topol[::-1].find('.')-1]       # file=FOO.N.pdb -> name=FOO.N
            self.num_atoms=len(self.atom)
            self.dimensionality=self.num_atoms*3
            self.num_residues=len(self.resid)
            self.num_waters=len(self.water)
            self.num_chains=len(self.chains)
            self.num_ions=len(self.ion)
            self.num_proteins=len(self.protein)
            self.num_lipids=len(self.lipid)
            self.num_molecules=len(self.molecule)
            self.num_ligands=len(self.ligand)
            self.list_atoms=[ii for ii in range(self.num_atoms)]


            ### Loading coordinates
            if self.file_topol_type in ['psf']: coors=False
            if coors:
                self.load_traj(self.file_topol,frame='ALL',wrap=wrap,verbose=False)


            if verbose:
                self.info()
                if coors:
                    self.traj[0].info(index=0)

        ## END of IF input_file


    # END OF INSTANTATION

    ###
    ### INTERNAL FUNCTIONS OF A SET
    ###

    # Info function

    def info(self):

        self.num_atoms=len(self.atom)
        self.num_residues=len(self.resid)
        self.num_chains=len(self.chain)
        self.num_waters=len(self.water)
        self.num_trajs=len(self.traj)
        self.num_ions=len(self.ion)
        self.num_proteins=len(self.protein)
        self.num_lipids=len(self.lipid)
        self.num_ligands=len(self.ligand)
        self.num_molecules=len(self.molecule)
        print '#','System created from the file',self.file_topol,':'
        print '#',self.num_atoms,' atoms'
        print '#',self.num_residues,' residues'
        print '#',self.num_chains,' chains'
        if self.num_waters:    print '#',self.num_waters,' waters'
        if self.num_ions:      print '#',self.num_ions,' ions'
        if self.num_lipids:    print '#',self.num_lipids,' lipids'
        if self.num_proteins:  print '#',self.num_proteins,' proteins'
        if self.num_ligands:   print '#',self.num_ligands,' ligands'
        if self.num_molecules: print '#',self.num_molecules,' other molecules'

    # To handle files

    def load_topol (self,name_file):

        if name_file.endswith('pdb'):
            fff=open(name_file,'r')
            for line in fff:
                ss=line.split()[0]
                if ss in ['HEADER','TITLE','CRYST1']: self.pdb_header.append(line)
                if ss.startswith('END'): break  # To read only the 1st model
                if ss in ['HELIX','SHEET','TURN']: self.pdb_ss.append(line)
                if ss in ['ATOM','HETATM']:
                    temp_atom=cl_unit()
                    temp_atom.type_pdb=line[0:6].replace(' ', '')
                    temp_atom.pdb_index=int(line[6:11])
                    temp_atom.name=(line[12:16].split())[0]
                    temp_atom.alt_loc=line[16]
                    temp_atom.resid.name=(line[17:21]).replace(' ', '') # real format: 17:20
                    temp_atom.chain.name=line[21]
                    temp_atom.resid.pdb_index=int(line[22:26])
                    temp_atom.code_ins_res=line[26]
                    temp_atom.occup=float(line[54:60])
                    temp_atom.bfactor=float(line[60:66])
                    temp_atom.seg_ident=line[72:76].replace(' ', '')
                    temp_atom.elem_symb=line[76:78].replace(' ', '')
                    temp_atom.charge=line[78:80].replace(' ', '')
                    temp_atom.index=len(self.atom)
                    self.atom.append(temp_atom)

            fff.close()

        if name_file.endswith('gro'):

            fff=open(name_file,'r')
            line=fff.readline()                                          # Header of the gro file
            line=fff.readline()                                        
            self.num_atoms=int(line)

            for i in range(self.num_atoms):           
            ## Fixed format taken from http://manual.gromacs.org/online/gro.html
                temp_atom=cl_unit()
                line=fff.readline()
                temp_atom.pdb_index=int(line[15:20])
                temp_atom.name=line[10:15].replace(" ", "")
                temp_atom.resid.name=line[5:10].replace(" ", "")
                temp_atom.resid.pdb_index=int(line[0:5]) 
                temp_atom.index=i           
                temp_atom.chain.name='A'
                self.atom.append(temp_atom)

            fff.close()

        if name_file.endswith('psf'):

            fff=open(name_file,'r')

            line=fff.readline()
            if not line.startswith('PSF'):
                print '# Error: uknown PSF format'
                return

            line=fff.readline()
            line=fff.readline()
            num_head_lines=int(line.split()[0])
            for ii in range(num_head_lines):
                line=fff.readline()

            line=fff.readline()

            line=fff.readline()
            self.num_atoms=int(line.split()[0])
            if line.split()[1]!='!NATOM':
                print '# Error: uknown PSF format'
                return

            for ii in range(self.num_atoms):
                temp_atom=cl_unit()
                line=fff.readline()
                xx=line.split()
                temp_atom.pdb_index=int(xx[0])
                temp_atom.chain.name=xx[1]
                temp_atom.resid.pdb_index=int(xx[2]) 
                temp_atom.resid.name=xx[3]
                temp_atom.name=xx[4]
                temp_atom.type_pdb=xx[5]
                temp_atom.type=xx[5]
                temp_atom.charge=float(xx[6])
                temp_atom.mass=float(xx[7])
                self.atom.append(temp_atom)

            line=fff.readline()
            line=fff.readline()
            numbonds=int(line.split()[0])
            if line.split()[1]!='!NBOND:':
                print '# Error: uknown PSF format'
                return
            aa=0
            while aa<numbonds:
                line=fff.readline()
                bb=0
                xx=line.split()
                cc=len(xx)
                while bb<cc:
                    b1=int(xx[bb])-1
                    b2=int(xx[bb+1])-1
                    bb+=2
                    self.atom[b1].covalent_bonds.append(b2)
                    self.atom[b2].covalent_bonds.append(b1)
                aa+=bb

    def write (self,filename=None,traj=0,frame=0,sel='ALL'):
        
        if len(self.traj)==0:
            print '# Error: The system has no coordinates to export. '
            return

        select,nlist,nsys=__read_set_opt__(self,sel)
        
        select.sort()

        if filename==None:

            print 'Enter filename: '
            print '      foo.write("foo.pdb")'

        elif filename.endswith('.pdb'):

            if path.exists(filename): 
                print '# The file '+filename+' already exists.'
                return
     
            fff=open(filename,'w')
     
            a='HEADER    '+'> FILE CREATED BY AQUALAB '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+' <\n'
            fff.write(str(a))
     
            for ii in self.pdb_header: fff.write(str(ii))
            
            for ii in self.pdb_ss:     fff.write(str(ii))
     
            a='CRYST1'  #1-6
            a+="%9.3f" % float(self.traj[traj].frame[frame].cell[0,0]) #7-15
            a+="%9.3f" % float(self.traj[traj].frame[frame].cell[1,1]) #16-24
            a+="%9.3f" % float(self.traj[traj].frame[frame].cell[2,2]) #25-33
            a+="%7.2f" % float(self.traj[traj].frame[frame].cell[0,1]) #34-40
            a+="%7.2f" % float(self.traj[traj].frame[frame].cell[0,2]) #41-47
            a+="%7.2f" % float(self.traj[traj].frame[frame].cell[1,2]) #48-54
            a+="  "
            a+="%-10s" % "P 1"   #56-66
            a+="%-3s"  % "1"     #67-70
            a+='\n' 
            fff.write(str(a))

            dct_aux={'ATOM': 'ATOM  ', 'HETATM': 'HETATM'}
            
            new_index=0
            for jj in range(nlist):
                ii=select[jj]
                new_index+=1
                try:
                    a=dct_aux[self.atom[ii].type_pdb]          # 1-6
                except:
                    a='ATOM  '                                 # 1-6
                a+="%5d" % (new_index)                         # 7-11
                #a+="%5d" % self.atom[ii].pdb_index            # 7-11
                a+=' '                                         # 12
                if len(self.atom[ii].name)<4:
                    a+=' '+"%-3s" % self.atom[ii].name             # 13-16
                else:
                    a+=self.atom[ii].name[3]+self.atom[ii].name[:3] # 13-16
                a+=' '                                         # 17
                a+="%3s" % self.atom[ii].resid.name[:3]            # 18-20
                a+=' '                                         # 21
                a+="%1s" % self.atom[ii].chain.name            # 22
                a+="%4d" % (self.atom[ii].resid.pdb_index % 10000)       # 23-26
                a+=' '                                         # 27
                a+='   '                                       # 28-30
                a+="%8.3f" % float(self.traj[traj].frame[frame].coors[ii][0])   # 31-38
                a+="%8.3f" % float(self.traj[traj].frame[frame].coors[ii][1])   # 39-46
                a+="%8.3f" % float(self.traj[traj].frame[frame].coors[ii][2])   # 47-54
                a+="%6.2f" % self.atom[ii].occup               # 55-60
                a+="%6.2f" % self.atom[ii].bfactor             # 61-66
                a+='          '                                # 67-76
                a+="%2s" % self.atom[ii].elem_symb             # 77-78
                #a+="%2s" % self.atom[ii].charge                # 79-80
                a+='\n' 
                fff.write(str(a))         
                # if ii<(nlist-1):
                #     if self.atom[ii].type_pdb!=self.atom[ii+1].type_pdb or self.atom[ii].chain.name!=self.atom[ii+1].chain.name :
                #         new_index+=1
                #         a="TER   "
                #         a+="%5d" % (new_index)
                #         a+=' '
                #         a+='  '
                #         a+=' '                                         
                #         a+="%3s" % self.atom[ii].resid.name            
                #         a+='\n' 
                #         fff.write(str(a))
            a='END   '+'\n'
            fff.write(str(a))
            fff.close()

        elif filename.endswith('.gro'):

            if path.exists(filename): 
                print '# The file '+filename+' already exists.'
                return
     
            fff=open(filename,'w')
     
            a='FILE CREATED BY AQUALAB '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+'\n'
            fff.write(str(a))
     
            a=str(nlist)+'\n'
            fff.write(str(a))

            new_index=0
            for jj in range(nlist):
                ii=select[jj]
                new_index+=1
                a="%5d" % self.atom[ii].resid.pdb_index     # 1-5
                a+="%-5s" % self.atom[ii].resid.name        # 6-10
                a+="%5s"  % self.atom[ii].name              # 11-15
                a+="%5d" % new_index                         # 16-20
                a+="%8.3f" % float(self.traj[traj].frame[frame].coors[ii][0]/10.0)   # 
                a+="%8.3f" % float(self.traj[traj].frame[frame].coors[ii][1]/10.0)   # 
                a+="%8.3f" % float(self.traj[traj].frame[frame].coors[ii][2]/10.0)   # 
                a+="%8.4f" % float(0.0)   # 
                a+="%8.4f" % float(0.0)   # 
                a+="%8.4f" % float(0.0)   # 
                a+='\n' 
                fff.write(str(a))         
                
            a=str(float(self.traj[traj].frame[frame].box[0,0]/10.0))
            a+=' '+str(float(self.traj[traj].frame[frame].box[1,1]/10.0))
            a+=' '+str(float(self.traj[traj].frame[frame].box[2,2]/10.0))
            a+='\n'

            fff.write(str(a))
            fff.close()

        return None


    def write_set_to_file(self,name_of_file):
        fff=open(name_of_file,'w')
        pic.dump(self,fff)
        fff.close()

    def read_set_from_file(self,name_of_file):
        fff=open(name_of_file,'r')
        A=pic.load(fff)
        fff.close()
        return A

    def add_covalent_bonds(self,bonds=None,verbose=False):
        a = numpy.array(bonds)
        bad = False
        ashape=a.shape
        if len(ashape)==1:
            if ashape[0]==2:
                aa=a[0]
                bb=a[1]
                self.atom[aa].covalent_bonds.append(bb)
                self.atom[bb].covalent_bonds.append(aa)
            else:
                bad = True
        elif len(ashape)>1:
            if ashape[1]==2:
                for covab in bonds:
                    aa=covab[0]
                    bb=covab[1]
                    self.atom[aa].covalent_bonds.append(bb)
                    self.atom[bb].covalent_bonds.append(aa)
            else:
                bad = True
        else:
            bad = True

        if bad:
            print 'Error: list of atom pairs required.'
        pass

    def add_donors(self,select=None,verbose=False):
        setdon,nlist,numsys=__read_set_opt__(self,select)
        self.donors=list(self.donors)
        num_donors_added=0
        for ii in setdon:
            self.atom[ii].donor=True
            with_h=False
            for jj in self.atom[ii].covalent_bonds:
                if self.atom[jj].type=='H':
                    with_h=True
                    self.donors.append([ii,jj])
                    num_donors_added+=1
            if not with_h:
                print '# No H atom bonded to',self.atom[ii].info()
        self.donors=numpy.array(self.donors,dtype=int,order='Fortran')
        self.donors=self.donors[numpy.lexsort(numpy.fliplr(self.donors).T)]
        if verbose:
            print '#',num_donors_added,'donors added.'

    def add_acceptors(self,select=None,verbose=False):
        setacc,nlist,numsys=__read_set_opt__(self,select)
        self.acceptors=list(self.acceptors)
        for ii in setacc:
            if ii not in self.acceptors:
                self.atom[ii].acceptor=True
                self.acceptors.append(ii)
        self.acceptors=numpy.array(self.acceptors,dtype=int,order='Fortran')
        if verbose:
            print '#',nlist,'acceptors added.'

    # To handle coordinates

    def info_trajs(self):
        if len(self.traj):
            for aa in range(len(self.traj)):
                self.traj[aa].info(index=aa)
        else:
            print '# No coordinates'
        pass

    def coors2frame (self,traj=0,frame=0):

        self.which_traj=traj
        self.which_frame=frame
        self.box=self.traj[traj].frame[frame].box
        for ii in range(self.num_atoms):
            self.atom[ii].coors=self.traj[traj].frame[frame].coors[ii,:]

        pass

    def load_traj (self,file_input=None,frame=None,wrap=True,verbose=False):

        temp_traj=traj(file_input,frame,wrap,verbose=False)
        if verbose:
            temp_traj.info(index=len(self.traj))
        self.traj.append(temp_traj)
        del(temp_traj)
        if len(self.traj)==1:
            if len(self.traj[0].frame) :
                self.coors2frame()

    def delete_traj (self,index='ALL'):

        if index in ['all','All','ALL']:
            for ii in self.traj:
                if ii.io_opened:
                    ii.close()
            del(self.traj);self.traj=[]
            del(self.box);self.box=[]
            for ii in range(self.num_atoms):
                del(self.atom[ii].coors); self.atom[ii].coors=[]

        elif type(index) in [int]:
            self.traj.__delitem__(index)
            if self.which_traj==index:
                del(self.box);self.box=[]
                for ii in range(self.num_atoms):
                    del(self.atom[ii].coors)
                    self.atom[ii].coors=[]
            elif self.which_traj>index:
                self.which_traj-=1

        pass

    def selection (self,condition=None,traj=0,frame='ALL',pbc=True):

        if type(condition) in [tuple,list]:
            return condition
        else:
            list_condition=selection(self,condition,traj,frame,pbc)
            return list_condition

    def selection_covalent_chains(self,chain=None,select='protein'):

        salida=selection_covalent_chains(system=self,chain=chain,select=select)

        return salida

    def selection_hbonds(self,setA='ALL',verbose=False):
     
        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_system,diff_set=__read_sets_opt__(self,setA,None,None)

        don=[]
        don_H=[]
        don_start_H=[0]
        acc=[]
        acc_H=[]         # solo para skiner, topological, etc..
        acc_start_H=[0]
        all_wat=1

        gg_a=0
        gg_d=0
        for ii in setA:
            atom=self.atom[ii]
            if atom.resid.type!='Water': all_wat=0
            if atom.donor and atom.acceptor:
                don.append(ii)
                acc.append(ii)
                for jj in atom.covalent_bonds:
                    if self.atom[jj].type=='H':
                        don_H.append(jj)
                        gg_d+=1
                        acc_H.append(jj)
                        gg_a+=1
                don_start_H.append(gg_d)
                acc_start_H.append(gg_a)
            elif atom.donor:
                don.append(ii)
                for jj in atom.covalent_bonds:
                    if self.atom[jj].type=='H':
                        don_H.append(jj)
                        gg_d+=1
                don_start_H.append(gg_d)
            elif atom.acceptor:
                acc.append(ii)
                for jj in atom.covalent_bonds:
                    if self.atom[jj].type=='H':
                        acc_H.append(jj)
                        gg_a+=1
                acc_start_H.append(gg_a)                

        if verbose:
            print '# [ Donor, Hydrogen ]'
            for ii in range(len(don)):
                adon=self.atom[don[ii]]
                adonstr=adon.name+'-'+str(adon.pdb_index)+'/'+adon.resid.name+'-'+str(adon.resid.pdb_index)
                for jj in range(don_start_H[ii],don_start_H[ii+1]):
                    adon_H=self.atom[don_H[jj]]
                    adon_Hstr=adon_H.name+'-'+str(adon_H.pdb_index)+'/'+adon_H.resid.name+'-'+str(adon_H.resid.pdb_index)
                    print don[ii],'\t',adonstr,'\t',don_H[jj],'\t',adon_Hstr
            print ' '
            print '# [ Acceptor ]'
            for ii in acc:
                aacc=self.atom[ii]
                aaccstr=aacc.name+'-'+str(aacc.pdb_index)+'/'+aacc.resid.name+'-'+str(aacc.resid.pdb_index)
                print ii,'\t',aaccstr

        return [numpy.array(acc,order='F'),numpy.array(acc_start_H,order='F'),numpy.array(acc_H,order='F'), \
                numpy.array(don,order='F'),numpy.array(don_start_H,order='F'),numpy.array(don_H,order='F'),all_wat]

    def extract(self,select=None,verbose=False):

        if type(select) not in [list,tuple]:
            select=self.selection(select) 

        tmp_msystem=msystem()

        tmp_msystem.ff=self.ff
        tmp_msystem.file_topol=self.file_topol
        tmp_msystem.file_topol_type=self.file_topol_type
        tmp_msystem.parent=self
        tmp_msystem.parent_selection=select

        for ii in select:
            tmp_msystem.atom.append(ccopy.deepcopy(self.atom[ii]))

        # Global attributes
       
        tmp_msystem.name=self.name       # file=FOO.N.pdb -> name=FOO.N
        tmp_msystem.num_atoms=len(tmp_msystem.atom)
        tmp_msystem.dimensionality=tmp_msystem.num_atoms*3
        tmp_msystem.num_residues=len(tmp_msystem.resid)
        tmp_msystem.num_waters=len(tmp_msystem.water)
        tmp_msystem.num_chains=len(tmp_msystem.chains)
        tmp_msystem.num_ions=len(tmp_msystem.ion)
        tmp_msystem.num_proteins=len(tmp_msystem.protein)
        tmp_msystem.num_lipids=len(tmp_msystem.lipid)
        tmp_msystem.num_molecules=len(tmp_msystem.molecule)
        tmp_msystem.num_ligands=len(tmp_msystem.ligand)
        tmp_msystem.list_atoms=[ii for ii in range(tmp_msystem.num_atoms)]

        # Coordinates

        tmp_msystem.traj=ccopy.deepcopy(self.traj)
        for traj_id in range(len(tmp_msystem.traj)):
            for frame_id in range(len(self.traj[traj_id].frame)):
                tmp_msystem.traj[traj_id].frame[frame_id].coors=self.traj[traj_id].frame[frame_id].coors[select,:]

        if verbose:
            tmp_msystem.info()
            if tmp_msystem.traj:
                tmp_msystem.traj[0].info(index=0)

        return tmp_msystem 

    def rebuild_chains(self, list_atoms=None, names=None, pdb_names=None, types=None, verbose=False):

        self.num_chains=0
        self.chain=[]
        self.chains=[]
        for ii in self.atom:
            ii.chain=labels_unit()
        aa=0
        for ii in xrange(len(list_atoms)):
            temp_chain=labels_set()
            if names==None:
                bb=str(aa)
            else:
                bb=str(names[ii])
            temp_chain.index=aa
            temp_chain.name=bb
            temp_chain.list_atoms=list_atoms[ii]
            temp_chain.num_atoms=len(list_atoms[ii])
            list_resids={}
            for jj in temp_chain.list_atoms:
                self.atom[jj].chain.index=aa
                self.atom[jj].chain.name=bb
                list_resids[self.atom[jj].resid.index]=0
                cc=self.atom[jj].resid.type
            temp_chain.type=cc
            for jj in list_resids.keys():
                self.resid[jj].type=cc
                self.resid[jj].chain=temp_chain
            self.chain.append(temp_chain)
            aa+=1

        self.protein=[]
        for ii in self.chain:
            self.chains.append(ii.name)
            if ii.type=='Protein':
                self.protein.append(ii.list_atoms)

        if verbose:
            self.info()
                

###############################################################
###############################################################
    # To handle the set
    def center_of_mass(self,select='ALL',mass_weighted=False,traj=0,frame=0,pbc=True):

        pbc_opt=0
        if pbc:
            pbc_opt=1

        setcom,nlist_setcom,numsys=__read_set_opt__(self,select)

        num_frames=__length_frame_opt__(self,traj,frame)
        com=numpy.empty(shape=(num_frames,3),dtype=float,order='Fortran')

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            com[num_frames,:]=faux.center_of_mass(pbc_opt,setcom,iframe.coors,iframe.box,iframe.orthogonal,\
                                 nlist_setcom,numsys)
            num_frames+=1

        if num_frames==1:
            return com[0,:]
        else:
            return com


    def center(self,center_of=None,select='ALL',mass_weighted=False,traj=0,frame=0,pbc=True,wrap=True):

        ## I have to correct here the split molecules.

        if center_of==None:
            print 'Selection in input variable "center_of" required.'
            return

        if mass_weighted:
            print 'Mass_weighted option not implemented yet.'
            return

        pbc_opt=0
        if pbc:
            pbc_opt=1        

        wrap_opt=0
        if wrap:
            wrap_opt=1

        setcom,nlist_setcom,numsys=__read_set_opt__(self,center_of)
        setmov,nlist_setmov,numsys=__read_set_opt__(self,select)

        for iframe in __read_frame_opt__(self,traj,frame):
            iframe.coors=numpy.array(iframe.coors,order='Fotran')
            faux.center(pbc_opt,wrap_opt,setcom,setmov,iframe.coors,iframe.box,iframe.orthogonal,\
                                 nlist_setcom,nlist_setmov,numsys)
            
        pass

    


    def distance(self,sel1='ALL',sel2=None,points=None,traj=0,frame='ALL',legend=False,pbc=True):
        
        if pbc:
            #check_cell=self.traj[traj].frame[0].cell
            #if check_cell[0,1]!=90 or check_cell[0,2]!=90 or check_cell[1,2]!=90:
            #    print '# PBC not implemented for not orthorhombic boxes'
            #    return
            pbc=1

        num_frames=__length_frame_opt__(self,traj,frame)

        if points!=None:
            if type(points)==numpy.ndarray:
                setA,nlist_A,nsys_A=__read_set_opt__(self,sel1)
                points_aux=ccopy.deepcopy(points)
                if len(points.shape)==1:
                    num_points=1
                    points_aux.resize((1,3))
                elif points.shape[1]!=3:
                    print '# points needs to be a numpy.array with shape (3) or (N,3)'
                    return
                dists=numpy.empty(shape=(num_frames,nlist_A,num_points),dtype=float,order='Fortran')
                num_frames=0
                for iframe in __read_frame_opt__(self,traj,frame):
                    dists[num_frames,:,:]=faux.distance_p(pbc,setA,iframe.coors,iframe.box,iframe.orthogonal,points_aux,nlist_A,num_points,nsys_A)
                    num_frames+=1

                if num_frames==1:
                    if num_points==1:
                        return dists[0,:,0]
                    else:
                        return dists[0,:,0]
                else:
                    if num_points==1:
                        return dists[:,:,0]
                    else:
                        return dists
            else:
                print '# points needs to be a numpy.array with shape [3] or [N,3]'
        else:

            setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,sel1,None,sel2)
            
            dists=numpy.empty(shape=(num_frames,nlist_A,nlist_B),dtype=float,order='Fortran')
            
            num_frames=0
            for iframe in __read_frame_opt__(self,traj,frame):
                dists[num_frames,:,:]=faux.distance(diff_syst,diff_set,pbc,setA,iframe.coors,iframe.box,iframe.invbox,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
                num_frames+=1

            if legend:
                if num_frames==1:
                    return dists[0,:,:], setA, setB
                else:
                    return dists, setA, setB
            else:
                if num_frames==1:
                    return dists[0,:,:]
                else:
                    return dists


    def distance_image_pbc(self,setA='ALL',setB=None,traj=0,frame='ALL'):

        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,setA,None,setB)
        num_frames=__length_frame_opt__(self,traj,frame)
        min_dists=numpy.empty(shape=(nlist_A,num_frames),dtype=float,order='Fortran')
        ind_atoms_min=numpy.empty(shape=(nlist_A,num_frames),dtype=int,order='Fortran')
        min_image=numpy.empty(shape=(nlist_A,3,num_frames),dtype=int,order='Fortran')

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            min_dists[:,num_frames],ind_atoms_min[:,num_frames],min_image[:,:,num_frames]=faux.distance_images(diff_syst,diff_set,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
            num_frames+=1

        if num_frames==1:
            return min_dists[:,0],ind_atoms_min[:,0],min_image[:,:,0]
        else:
            return min_dists, ind_atoms_min,min_image



    def radius_gyration(self,setA='ALL',traj=0,frame='ALL',pbc=False):

        setA,nlist_A,nsys_A=__read_set_opt__(self,setA)
        num_frames=__length_frame_opt__(self,traj,frame)
        rgs=numpy.empty(shape=(num_frames),dtype=float,order='Fortran')

        pbc_opt=0
        if pbc:
            pbc_opt=1

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            rgs[num_frames]=faux.radius_gyration(pbc_opt,setA,iframe.coors,iframe.box,iframe.orthogonal,nlist_A,nsys_A)
            num_frames+=1

        if num_frames==1:
            return rgs[0]
        else:
            return rgs


    def principal_inertia_axis(self,setA='ALL',traj=0,frame='ALL'):

        setA,nlist_A,nsys_A=__read_set_opt__(self,setA)
        num_frames=__length_frame_opt__(self,traj,frame)
        piaxis=numpy.empty(shape=(3,3,num_frames),dtype=float,order='Fortran')

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            piaxis[:,:,num_frames]=faux.principal_inertia_axis(setA,iframe.coors,iframe.box,iframe.orthogonal,nlist_A,nsys_A)
            num_frames+=1

        if num_frames==1:
            return piaxis[:,:,0]
        else:
            return piaxis

    def principal_geometric_axis(self,setA='ALL',traj=0,frame='ALL'):

        setA,nlist_A,nsys_A=__read_set_opt__(self,setA)
        num_frames=__length_frame_opt__(self,traj,frame)
        pgaxis=numpy.empty(shape=(3,3,num_frames),dtype=float,order='Fortran')

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            pgaxis[:,:,num_frames]=faux.principal_geometric_axis(setA,iframe.coors,iframe.box,iframe.orthogonal,nlist_A,nsys_A)
            num_frames+=1

        if num_frames==1:
            return pgaxis[:,:,0]
        else:
            return pgaxis


    def ramachandran_map(self,resid='ALL',traj=0,frame='ALL',pdb_index=False,legend=False):

        if resid in ['ALL','All','all']:
            list_phi=selection_covalent_chains(system=self,chain=['C','N','CA','C'],select='protein')
            list_psi=selection_covalent_chains(system=self,chain=['N','CA','C','N'],select='protein')
            
        else:
            if pdb_index:
                sel_prefix='resid.pdb_index in '
            else:
                sel_prefix='resid.index in '
            if type(resid) in [int]:
                resid=[resid]
            aux_phi={}
            aux_psi={}
            for ii in resid:
                aux_phi[ii]=''
                aux_psi[ii]=''
                aux_phi[ii-1]=''
                aux_psi[ii+1]=''
            sel_phi=sel_prefix
            sel_psi=sel_prefix
            for ii in aux_phi.keys():
                sel_phi+=str(ii)+' '
            for ii in aux_psi.keys():
                sel_psi+=str(ii)+' '
            list_phi=selection_covalent_chains(system=self,chain=['C','N','CA','C'],select=sel_phi)
            list_psi=selection_covalent_chains(system=self,chain=['N','CA','C','N'],select=sel_psi)

        aux_tot={}
        aux_phi=[]
        aux_psi=[]
        for ii in list_phi:
            jj=self.atom[ii[2]].resid.index
            aux_tot[jj]=''
            aux_phi.append(jj)
        for ii in list_psi:
            jj=self.atom[ii[1]].resid.index
            aux_tot[jj]=''
            aux_psi.append(jj)

        list_resids=numpy.sort(aux_tot.keys())
        num_resid=list_resids.shape[0]
        num_phis=len(list_phi)
        num_psis=len(list_psi)
        aux_tot={}
        for ii in range(num_resid):
            aux_tot[list_resids[ii]]=ii
        for ii in range(num_phis):
            aux_phi[ii]=aux_tot[aux_phi[ii]]
        for ii in range(num_psis):
            aux_psi[ii]=aux_tot[aux_psi[ii]]

        
        if legend:
            if pdb_index:
                key_angs=[['Phi '+str(self.resid[ii].pdb_index),'Psi '+str(self.resid[ii].pdb_index)] for ii in list_resids]
            else:
                key_angs=[['Phi '+str(ii),'Psi '+str(ii)] for ii in list_resids]

        num_frames=__length_frame_opt__(self,traj,frame)
        dih_angs=numpy.zeros(shape=(num_frames,num_resid,2),dtype=float,order='Fortran')
        jj=-1
        for iframe in __read_frame_opt__(self,traj,frame):
            jj+=1
            if num_phis:
                phis=faux.dihedral_angles(iframe.coors,iframe.box,iframe.orthogonal,list_phi,num_phis,self.num_atoms)
                for kk in range(num_phis):
                    dih_angs[jj,aux_phi[kk],0]=phis[kk]
            if num_psis:
                psis=faux.dihedral_angles(iframe.coors,iframe.box,iframe.orthogonal,list_psi,num_psis,self.num_atoms)
                for kk in range(num_psis):
                    dih_angs[jj,aux_psi[kk],1]=psis[kk]

        if legend:
            if num_frames==1:
                if num_resid==1:
                    return dih_angs[0,0,:],key_angs[0]
                else:
                    return dih_angs[0,:,:],key_angs
            else:
                if num_resid==1:
                    return dih_angs[:,0,:],key_angs[0]
                else:
                    return dih_angs,key_angs
        else:
            if num_frames==1:
                if num_resid==1:
                    return dih_angs[0,0,:]
                else:
                    return dih_angs[0,:,:]
            else:
                if num_resid==1:
                    return dih_angs[:,0,:]
                else:
                    return dih_angs

    def dihedral_angle(self,covalent_chain=None,traj=0,frame='ALL',pbc=True):

        if not covalent_chain:
            print '# Error: Check method msystem.selection_covalent_chains()'
            return

        covalent_chain=numpy.array(covalent_chain,dtype=int,order='F')

        if covalent_chain.shape[-1]!=4:
            print '# Error: 4 atoms per covalent_chain'
            return

        if len(covalent_chain.shape)==1:
            covalent_chain.resize((1,4))

        num_dih_angs=covalent_chain.shape[0]

        num_frames=__length_frame_opt__(self,traj,frame)
        dih_angs=numpy.empty(shape=(num_frames,num_dih_angs),dtype=float,order='Fortran')
        
        jj=-1
        for iframe in __read_frame_opt__(self,traj,frame):
            jj+=1
            dih_angs[jj,:]=faux.dihedral_angles(pbc,iframe.coors,iframe.box,iframe.invbox,iframe.orthogonal,covalent_chain,num_dih_angs,self.num_atoms)

        if num_frames==1:
            if num_dih_angs==1:
                return dih_angs[0,0]
            else:
                return dih_angs[0,:]
        else:
            if num_dih_angs==1:
                return dih_angs[:,0]
            else:
                return dih_angs

    def rmsd(self,selection=None,traj_ref=None,frame_ref=None,traj=0,frame='ALL'):

         setA,n_A,natoms_A=__read_set_opt__(self,selection)

         num_frames=__length_frame_opt__(self,traj,frame)
         coors_reference=self.traj[traj_ref].frame[frame_ref].coors
         rmsd_vals=numpy.zeros((num_frames),dtype=float,order='F')

         jj=0
         for iframe in __read_frame_opt__(self,traj,frame):
             rmsd_val=faux.rmsd(coors_reference,iframe.coors,setA,n_A,natoms_A)
             rmsd_vals[jj]=rmsd_val
             jj+=1

         return rmsd_vals

    def least_rmsd(self,msystem_ref=None,selection_ref='ALL',traj_ref=0,frame_ref=0,selection='ALL',traj=0,frame='ALL',pbc=True):

        '''output should be the least rmsd and in addition and optionally, the translation and rotation.'''
        # Esto de wrap un wrap, es correcto pero se puede economizar mucho poniendolo como opcion de frame.
        # Para poder hacer iframe.wrap(select,new) y se solucionaria el tema.. y  no este lio
        if msystem_ref==None:
            msystem_ref=self
            setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system,diff_set=__read_sets_opt__(self,selection_ref,None,selection)
        else:
            setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system,diff_set=__read_sets_opt__(msystem_ref,selection_ref,self,selection)

        if n_A!=n_B :
            print '# Error: Different number of atoms'
            return

        if pbc:
            selfisunwrap=self.isunwrap(setB)
            refisunwrap =msystem_ref.isunwrap(setA)
            if not refisunwrap:
                msystem_ref.unwrap(selection=setA,traj=traj_ref,frame=frame_ref)
            if not selfisunwrap:
                self.unwrap(selection=setB,traj=traj,frame=frame)

        coors_reference=msystem_ref.traj[traj_ref].frame[frame_ref].coors
        
        num_frames=__length_frame_opt__(self,traj,frame)
        rmsd_traj=numpy.zeros((num_frames),dtype=float)
        center_ref_traj=numpy.zeros((num_frames,3),dtype=float)
        center_orig_traj=numpy.zeros((num_frames,3),dtype=float)
        rot_traj=numpy.zeros((num_frames,3,3),dtype=float,order='F')

        # Tal y como esta aqui... coors_reference no cambian en frame a frame (este no es el caso general, pero bueno)
        # Asi que deberia hacer un unwrap ya desde el principio.

        # A is the ref, B is self
        

        jj=0
        for iframe in __read_frame_opt__(self,traj,frame):
            rot,center_ref,center_orig,rmsd,g=faux.min_rmsd(coors_reference,iframe.coors,setA,setB,n_A,natoms_A,n_B,natoms_B)
            rmsd_traj[jj]=rmsd
            center_ref_traj[jj,:]=center_ref
            center_orig_traj[jj,:]=center_orig
            rot_traj[jj,:,:]=rot

        if pbc:
            if not selfisunwrap:
                self.wrap(selection=setB,traj=traj,frame=frame)
            if not refisunwrap:
                msystem_ref.wrap(selection=setA,traj=traj_ref,frame=frame_ref)
        
        if num_frames==1:
            return rmsd, center_orig, center_ref, rot
        else:
            return rmsd_traj, center_orig_traj, center_ref_traj, rot_traj

    def unwrap(self, selection='ALL',min_image_system=None,min_image_selection=None, traj=0,frame=0,new=False):

        min_image_switch=False
        if (min_image_system or min_image_selection)!=None :
            min_image_switch=True
            if not min_image_system:
               min_image_system=self
            if not min_image_selection:
                min_image_selection='ALL'

        # deberia hacerlo mirando los covalent chains.. para mas adelante oye. estas covalent chains se deberian hacer al principio
        setA,n_A,natoms_A=__read_set_opt__(self,selection)
        for iframe in __read_frame_opt__(self,traj,frame):
            if not numpy.isfortran(iframe.coors):
                iframe.coors=numpy.asfortranarray(iframe.coors)
            faux.unwrap(iframe.coors,iframe.box,iframe.orthogonal,setA,n_A,natoms_A)
            if min_image_switch:
                center_ref=min_image_system.center_of_mass(select=min_image_selection)
                center_unw=self.center_of_mass(select=selection)
                center_min=faux.min_image_point(center_unw,center_ref,iframe.box,iframe.invbox,iframe.orthogonal)
                trans_vect=center_min-center_unw
                if numpy.dot(trans_vect,trans_vect)>0.001: # Esto necesita ser corregido definiendo algo como frame.box.min_image_distance que se calcule al inicializar
                    coors_out=iframe.coors[setA,:]+trans_vect
                    iframe.coors[setA,:]=coors_out

        pass

    def iswrap(self, selection='ALL', traj=0,frame=0):
        
        # deberia hacerlo mirando los covalent chains.. para mas adelante oye. estas covalent chains se deberian hacer al principio
        setA,n_A,natoms_A=__read_set_opt__(self,selection)
        itis_frame=[]
        for iframe in __read_frame_opt__(self,traj,frame):
            # if not numpy.isfortran(iframe.coors):
            #     iframe.coors=numpy.asfortranarray(iframe.coors)
            itis=faux.iswrap(iframe.coors,iframe.box,iframe.orthogonal,setA,n_A,natoms_A)
            if itis==1:
                itis_frame.append(True)
            else:
                itis_frame.append(False)

        if len(itis_frame)==1:
            return itis_frame[0]
        else:
            return itis_frame

    def isunwrap(self, selection='ALL', traj=0,frame=0):
        
        # deberia hacerlo mirando los covalent chains.. para mas adelante oye. estas covalent chains se deberian hacer al principio
        setA,n_A,natoms_A=__read_set_opt__(self,selection)
        itis_frame=[]
        for iframe in __read_frame_opt__(self,traj,frame):
            # if not numpy.isfortran(iframe.coors):
            #     iframe.coors=numpy.asfortranarray(iframe.coors)
            itis=faux.isunwrap(iframe.coors,iframe.box,iframe.orthogonal,setA,n_A,natoms_A)
            if itis==1:
                itis_frame.append(True)
            else:
                itis_frame.append(False)

        if len(itis_frame)==1:
            return itis_frame[0]
        else:
            return itis_frame


    
    def wrap(self, selection='ALL', traj=0,frame='ALL',new=False):

        if selection in ['ALL','All','all']:
            for iframe in __read_frame_opt__(self,traj,frame):
                iframe.wrap()
        else:
            setA,n_A,natoms_A=__read_set_opt__(self,selection)
            for iframe in __read_frame_opt__(self,traj,frame):
                if not numpy.isfortran(iframe.coors):
                    iframe.coors=numpy.asfortranarray(iframe.coors)
                faux.wrap(iframe.coors,iframe.box,iframe.orthogonal,setA,n_A,natoms_A)
        pass

    def translation(self,translation=None,selection='ALL',traj=0,frame='ALL',new=False,wrap=False):

        # Esto deberia de ser sustituido por una funcion que sea como parse_set_arguments o algo asi
        # y lo mismo para el frame deberia ser...

        if new:
            print'Not implemented yet'
            return

        setA,n_A,natoms_A=__read_set_opt__(self,selection)
        for iframe in __read_frame_opt__(self,traj,frame):  ## Esto habria que hacerlo con mask... y ademas habria que hacerlo con una funcion sobre traj o frame o coors...
            coors_out=iframe.coors[setA,:]+translation
            iframe.coors[setA,:]=coors_out 

        pass

    def rotation(self,rotation_point=None,rotation_matrix=None,selection='ALL',traj=0,frame='ALL',new=False):

        # Esto deberia de ser sustituido por una funcion que sea como parse_set_arguments o algo asi
        # y lo mismo para el frame deberia ser...

        # if wrap:
        #     print'wrapping after rotation is risky... the box has to rotate also, an it is not the case here'

        if new:
            print'Not implemented yet'
            return
        #todo esto se puede economizar mucho mucho con operaciones definidas sobre coors... o sobre frame con una seleccion con lista de enteros y listo
        #ademas mucho mas rapido si uso transformaciones como:
        #http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
        #https://github.com/ahojnnes/numpy-snippets/blob/master/transformations.py
        setA,n_A,natoms_A=__read_set_opt__(self,selection)
        for iframe in __read_frame_opt__(self,traj,frame):  ## Esto habria que hacerlo con mask... y ademas habria que hacerlo con una funcion sobre traj o frame o coors...
            centered_coors=iframe.coors[setA,:]-rotation_point
            coors_out=numpy.zeros(centered_coors.shape,dtype=float)
            for ii in xrange(centered_coors.shape[0]):
                coors_out[ii,:]=numpy.matmul(rotation_matrix.T,centered_coors[ii,:]) # Funciona bien asi por como sale la matriz de giro de fortran... pero no se si vale para caso general
            iframe.coors[setA,:]=coors_out+rotation_point ## comparar que es mas eficiente, esto o trabajar con masks 

        pass


    def least_rmsd_fit(self,msystem_ref=None,selection_ref='ALL',traj_ref=0,frame_ref=0,selection='ALL',traj=0,frame='ALL',selection2move='ALL',new=False,pbc=True):

        # wrap after rotation makes no sense, it is not posible. The box show rotate also, but the box' axis have to be parallel to XYZ to be efficient.
        if new:
            print'not implemented yet'
            return

        #si las operaciones estuvieran definidas directamente sobre coors o frame... osea, sobre frame seguramente sea lo mejor, pues esto se podria economizar
        rmsd_traj, center_orig_traj, center_ref_traj, rot_traj=self.least_rmsd(msystem_ref,selection_ref,traj_ref,frame_ref,selection,traj,frame,pbc)
        self.rotation(center_orig_traj,rot_traj,selection2move,traj,frame,new)
        self.translation((center_ref_traj-center_orig_traj),selection2move,traj,frame,new,wrap=False)
        return rmsd_traj, center_orig_traj, center_ref_traj, rot_traj

#def min_distance(system,set_a,set_b=None,pbc=True,type_b='atoms'):
# 
#    if set_b==None:
# 
#        ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,True,system.frame[0].coors,system.frame[0].box,set_a,set_a,system.num_atoms,len(set_a),len(set_a))
# 
#    else:
#        if type_b=='atoms':
#            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,False,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),len(set_b))
#        elif type_b=='vectors':
#            l_vects=shape(set_b)
#            if len(l_vects)==1:
#                set_b=[set_b]
#                l_vects=shape(set_b)
#            numpy.array(set_b,order='Fortran')
#            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms_ref(pbc,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),l_vects[0])
#    
#    return ind_a1,ind_a2,min_dist
# 
# 


    def rdf(self,setA=None,setB=None,traj=0,frame='ALL',pbc=True,bins=100,segment=None):

        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,setA,None,setB)

        rdf_tot=numpy.zeros(shape=(bins),dtype=float,order='Fortran')
        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            dist_frame=faux.distance(1,1,pbc,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
            rdf_frame=faux.rdf_frame(dist_frame,iframe.box,segment[0],segment[1],bins,nlist_A,nlist_B)
            rdf_tot+=rdf_frame
            num_frames+=1.0

        rdf_tot=rdf_tot/(num_frames*1.0)

        return rdf_tot


    def neighbors(self,setA=None,setB=None,ranking=1,dist=None,asbonds=False,traj=0,frame=0,pbc=True):
     
        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,setA,None,setB)
        num_frames=__length_frame_opt__(self,traj,frame)

        if dist==None:
            neighbs=[]
            neighbs=numpy.empty(shape=(num_frames,nlist_A,ranking),dtype=int,order='Fortran')
            num_frames=0
            for iframe in __read_frame_opt__(self,traj,frame):
                neighbs[num_frames,:,:]=faux.neighbs_ranking(diff_syst,diff_set,pbc,ranking,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
                num_frames+=1
            if num_frames==1:
                return neighbs[0][:,:]
            else:
                return neighbs
            
        else:
            neighbs=[]
            sort_opt=0
            if ranking:
                sort_opt=1
            for iframe in __read_frame_opt__(self,traj,frame):
                contact_map,num_neighbs,dist_matrix=faux.neighbs_dist(diff_syst,diff_set,pbc,dist,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
                if asbonds:
                    tot_num_neighbs=sum(num_neighbs)
                    aux_neighbs,aux_dists=faux.translate2bonds(setA,setB,contact_map,dist_matrix,tot_num_neighbs,nlist_A,nlist_B)
                    neighbs.append([aux_neighbs,aux_dists])
                else:
                    aux_neighbs=[]
                    for ii in range(nlist_A):
                        if num_neighbs[ii]:
                            neighbs_A=faux.translate_list(sort_opt,setB,contact_map[ii,:],dist_matrix[ii,:],num_neighbs[ii],nlist_B)
                            aux_neighbs.append(neighbs_A)
                        else:
                            aux_neighbs.append([])
                    neighbs.append(aux_neighbs)
            if num_frames==1:
                return neighbs[0]
            else:
                return neighbs

#    def fast_contact_map (self,setA=None,setB=None,dist=None,traj=0,frame=0,pbc=True,verlet_list_init=None):
# 
#        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,setA,None,setB)
#        num_frames=__length_frame_opt__(self,traj,frame)
# 
#        c_map=numpy.empty(shape=(num_frames,nlist_A,nlist_B),dtype=int,order='Fortran')
#        num_frames=0
#        for iframe in __read_frame_opt__(self,traj,frame):
#            c_map[num_frames][:,:],aa,bb=faux.neighbs_dist_ns_list(diff_syst,diff_set,pbc,dist,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
#            num_frames+=1
#            del(aa); del(bb)
# 
#        if num_frames==1:
#            return c_map[0][:,:]
#        else:
#            return c_map

    def contact_map (self,setA=None,setB=None,dist=None,traj=0,frame=0,pbc=True):
     
        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,setA,None,setB)
        num_frames=__length_frame_opt__(self,traj,frame)

        c_map=numpy.empty(shape=(num_frames,nlist_A,nlist_B),dtype=int,order='Fortran')
        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            c_map[num_frames][:,:],aa,bb=faux.neighbs_dist(diff_syst,diff_set,pbc,dist,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
            num_frames+=1
            del(aa); del(bb)

        if num_frames==1:
            return c_map[0][:,:]
        else:
            return c_map
     
        #def plot_contact_map(contact_map):
        #    
        #    pylab.gray()
        #    pylab.imshow(contact_map==False,origin='lower',interpolation=None) # Would be better not to interpolate
        #    #pylab.matshow(contact_map==False)
        #    return pylab.show()



    def verlet_list_ns (self,r1=3.5,r2=6.0,traj=0,frame=0,pbc=True,update=False,verbose=False):

        print 'Not implemented as independent function yet'
        pass

    def grid_ns_list (self,rcell=6.0,setA='ALL',setB=None,traj=0,frame=0,pbc=True,verbose=False):
        
        print 'Not implemented as independent function yet'
        pass

    def make_cell_grid_ns(self,rcell=7.0,rcut=3.5,traj=0,frame=0):
        for iframe in __read_frame_opt__(self,traj,frame):
            faux.make_cell_ns(rcell,rcut,iframe.box,self.num_atoms)

    def verlet_list_grid_ns(self,r1=3.5,r2=7.0,rcell=7.0,traj=0,frame=0,iframe=None,pbc=True,update=False,verbose=False):

        pbc_opt=0
        if pbc:
            pbc_opt=1

        if iframe!=None:
            if update:
                faux.update_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
            else:
                #print 'aqui si'
                faux.make_cell_ns(rcell,r2,iframe.box,self.num_atoms)
                #print 'aqui tambien'
                faux.make_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
                #print 'aqui sale'
        else:
            if update:
                for iframe in __read_frame_opt__(self,traj,frame):
                    faux.update_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
            else:
                for iframe in __read_frame_opt__(self,traj,frame):
                    faux.make_cell_ns(rcell,r2,iframe.box,self.num_atoms)
                    faux.make_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
        
    def relative_water_position(self,pairs=None,pbc=True,verbose=False):

        traj=0
        frame=0

        pbc_opt=0
        if pbc:
            check_cell=self.traj[traj].frame[0].cell
            if check_cell[0,1]!=90 or check_cell[0,2]!=90 or check_cell[1,2]!=90:
                print '# PBC not implemented for not orthorhombic boxes'
                return
            pbc_opt=1

        numpairs=len(pairs)
        atsinds=numpy.empty((numpairs,6),dtype=int,order='Fortran')
        for ii in range(numpairs):
            w1=pairs[ii][0]; w2=pairs[ii][1]
            atsinds[ii,0]=self.water[w1].O.index
            atsinds[ii,1]=self.water[w1].H1.index
            atsinds[ii,2]=self.water[w1].H2.index
            atsinds[ii,3]=self.water[w2].O.index
            atsinds[ii,4]=self.water[w2].H1.index
            atsinds[ii,5]=self.water[w2].H2.index

        iframe=self.traj[0].frame[0]
        diff_syst=0
        diff_set=1
        carvars=faux.relative_water_position(diff_syst,diff_set,pbc_opt,atsinds,iframe.coors,iframe.box,iframe.orthogonal,numpairs,self.num_atoms)

        return carvars


    def hbonds (self,definition=None,set_A=None,set_B=None,acc_don_A=None,acc_don_B=None,traj=0,frame=0,sk_param=0.00850,roh_param=2.3000,roo_param=3.50,ang_param=30.0,optimize=False,pbc=True,infile=False,verbose=False):

        opt_effic=0
        opt_diff_syst=0
        opt_diff_set=1
        opt_pbc=0
        if pbc:
            opt_pbc=1


        if acc_don_A==None and acc_don_B==None:
            if set_A==None:
                print 'set_A and/or set_B needed'
                return
            else:
                acc_don_A=self.selection_hbonds(setA=set_A,verbose=False)
                if set_B==None:
                    acc_don_B=acc_don_A
                    opt_diff_set=0
                else:
                    acc_don_B=self.selection_hbonds(setA=set_B,verbose=False)
        else:
            if acc_don_B==None:
                acc_don_B=acc_don_A
                opt_diff_set=0


        nA_acc       = acc_don_A[0].shape[0]
        nA_acc_sH    = acc_don_A[1].shape[0] # Just for water and skinner, topological, etc...
        nA_acc_H     = acc_don_A[2].shape[0] # Just for water and skinner, topological, etc...
        nA_don       = acc_don_A[3].shape[0]
        nA_don_sH    = acc_don_A[4].shape[0]
        nA_don_H     = acc_don_A[5].shape[0]
        allwat_A     = acc_don_A[6]
        nB_acc       = acc_don_B[0].shape[0]
        nB_acc_sH    = acc_don_B[1].shape[0] # Just for water and skinner, topological, etc...
        nB_acc_H     = acc_don_B[2].shape[0] # Just for water and skinner, topological, etc...
        nB_don       = acc_don_B[3].shape[0]
        nB_don_sH    = acc_don_B[4].shape[0]
        nB_don_H     = acc_don_B[5].shape[0]
        allwat_B     = acc_don_B[6]


        num_frames=__length_frame_opt__(self,traj,frame)

        faux.hbdefinition=hbonds_type(definition,verbose=False)
        

        if faux.hbdefinition == 0 : 
            return
        

        # Skinner
        elif faux.hbdefinition == 1 : 
            faux.sk_param=sk_param
            cut_cell=-0.343*numpy.log(sk_param/7.10)+0.10
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
                return

            if infile:

                if type(frame) in [int]:
                    frame=[frame]

                if type(frame) not in [list,tuple]:
                    print '# "frame" must be a list'
                    return

                if optimize:
                    gg=0
                    hbout=[]
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=cut_cell,r2=cut_cell,rcell=cut_cell,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=cut_cell,r2=cut_cell,rcell=cut_cell,iframe=iframe,update=True)

                        faux.get_hbonds_skinner_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

                else:
                    hbout=[]
                    gg=0
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        faux.get_hbonds_skinner( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            else:
                if optimize:
                    gg=0
                    hbout=[]
                    for iframe in __read_frame_opt__(self,traj,frame):
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=cut_cell,r2=cut_cell,rcell=cut_cell,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=cut_cell,r2=cut_cell,rcell=cut_cell,iframe=iframe,update=True)

                        faux.get_hbonds_skinner_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.get_hbonds_skinner( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)
                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout


        # R(o,h)
        elif faux.hbdefinition == 2 : 
            faux.roh2_param= roh_param**2

            if infile:
                print 'Not implemented yet'
                pass

            else:

                if optimize:
                    gg=0
                    hbout=[]
                    for iframe in __read_frame_opt__(self,traj,frame):
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=roh_param,r2=roh_param,rcell=roh_param,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=roh_param,r2=roh_param,rcell=roh_param,iframe=iframe,update=True)

                        faux.get_hbonds_roh_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.get_hbonds_roh( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)
                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout

                
        # R(o,o)-Ang(o,o,h)
        elif faux.hbdefinition == 3 :
            faux.roo2_param, faux.cos_angooh_param= roo_param**2, numpy.cos(numpy.radians(ang_param))

            if infile:

                if type(frame) in [int]:
                    frame=[frame]

                if type(frame) not in [list,tuple]:
                    print '# "frame" must be a list'
                    return

                if optimize:
                    gg=0
                    hbout=[]
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe,update=True)

                        faux.get_hbonds_roo_angooh_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        faux.get_hbonds_roo_angooh( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            else:

                if optimize:
                    gg=0
                    hbout=[]
                    for iframe in __read_frame_opt__(self,traj,frame):
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe,update=True)

                        faux.get_hbonds_roo_angooh_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.get_hbonds_roo_angooh( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])
                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout

        # R(o,o)-Ang(o,h,o)
        elif faux.hbdefinition == 8 :
            faux.roo2_param, faux.cos_angoho_param= roo_param**2, numpy.cos(numpy.radians(ang_param))

            if infile:

                if type(frame) in [int]:
                    frame=[frame]

                if type(frame) not in [list,tuple]:
                    print '# "frame" must be a list'
                    return

                if optimize:
                    gg=0
                    hbout=[]
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe,update=True)

                        faux.get_hbonds_roo_angoho_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        faux.get_hbonds_roo_angoho( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            else:

                if optimize:
                    gg=0
                    hbout=[]
                    for iframe in __read_frame_opt__(self,traj,frame):
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe,update=True)

                        faux.get_hbonds_roo_angoho_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.get_hbonds_roo_angoho( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout


        # R(o,h)-Ang(o,h,o)
        elif faux.hbdefinition == 9 :
            faux.roh2_param, faux.cos_angoho_param= roh_param**2, numpy.cos(numpy.radians(ang_param))

            if infile:

                if type(frame) in [int]:
                    frame=[frame]

                if type(frame) not in [list,tuple]:
                    print '# "frame" must be a list'
                    return

                if optimize:
                    gg=0
                    hbout=[]
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=roh_param,r2=roh_param,rcell=roh_param,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=roh_param,r2=roh_param,rcell=roh_param,iframe=iframe,update=True)

                        faux.get_hbonds_roh_angoho_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        faux.get_hbonds_roh_angoho( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            else:

                if optimize:
                    gg=0
                    hbout=[]
                    for iframe in __read_frame_opt__(self,traj,frame):
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=roh_param,r2=roh_param,rcell=roh_param,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=roh_param,r2=roh_param,rcell=roh_param,iframe=iframe,update=True)

                        faux.get_hbonds_roh_angoho_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.get_hbonds_roh_angoho( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout


        # Donor-Acceptor-Number
        elif faux.hbdefinition == 4 : 
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
                return

            if infile:
                print 'Not implemented yet'
                pass

            else:

                if optimize:
                    gg=0
                    hbout=[]
                    for iframe in __read_frame_opt__(self,traj,frame):
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=6.0,r2=6.0,rcell=6.0,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=6.0,r2=6.0,rcell=6.0,iframe=iframe,update=True)

                        faux.get_hbonds_don_acc_num_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.get_hbonds_don_acc_num( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout


        # Topological
        elif faux.hbdefinition == 5 : 
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
                return

            if infile:
                print 'Not implemented yet'
                pass

            else:

                if optimize:
                    gg=0
                    hbout=[]
                    for iframe in __read_frame_opt__(self,traj,frame):
                        if (gg==0): 
                            self.verlet_list_grid_ns(r1=6.0,r2=6.0,rcell=6.0,iframe=iframe)
                        else:
                            self.verlet_list_grid_ns(r1=6.0,r2=6.0,rcell=6.0,iframe=iframe,update=True)

                        faux.get_hbonds_top_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.get_hbonds_top( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)

                        if faux.hbs_vals_out.shape[0]:
                            hbout.append([ccopy.deepcopy(faux.hbs_out),ccopy.deepcopy(faux.hbs_vals_out)])
                        else:
                            hbout.append([[],[]])

                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout


        # Donor-Number-Ang(o,o,h)
        elif faux.hbdefinition == 6 : 
            faux.cos_angooh_param= numpy.cos(numpy.radians(ang_param))
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
            print 'Not implemented yet'
            pass

        # Nearest-Neighbour
        elif faux.hbdefinition == 7 : 
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
            print 'Not implemented yet'
            pass

            #if verbose:
            #    print 'Listo'
        pass


#    def mss_hbonds_wat(self,definition=1,hbonds=None,bonds=None,verbose=True):
#        
#        mss_funcs.definition_hbs=faux.hbdefinition
# 
#        if hbonds==None:
#            print '# hbonds needed.'
#            return
# 
#        if definition==1:
# 
#            aux=numpy.zeros((self.num_atoms,3),dtype=int,order='F')
# 
#            for ii in range(len(self.water)):
#                jo=self.water[ii].O.index
#                jh1=self.water[ii].H1.index
#                jh2=self.water[ii].H2.index
#                aux[jo,0]=ii
#                aux[jo,1]=0
#                aux[jh1,0]=ii
#                aux[jh1,1]=1
#                aux[jh2,0]=ii
#                aux[jh2,1]=2
# 
# 
#            #for hbs in hbonds:
#            #    print hbs[1]
#            #print hbonds[0].shape,hbonds[1].shape
#            aux=numpy.array(aux,dtype=int,order='F')
# 
#            if type(hbonds[0][0][0]) in [numpy.ndarray]:
#                mss_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                mss_ind_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                for jj in range(len(hbonds)):
#                    num_hbs=hbonds[jj][0].shape[0]
#                    mss_ind=mss_funcs.ind_wat_limit_4_nosim(aux,hbonds[jj][0],hbonds[jj][1],self.num_waters,self.num_atoms,num_hbs)
#                    mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                    mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
#                    mss_tot[jj,:,:]=mss[:,:]
#                    mss_ind_tot[jj,:,:]=mss_ind[:,:]
#                return mss_tot,mss_ind_tot
# 
#            else:
#                num_hbs=hbonds[0].shape[0]
#                mss_ind=mss_funcs.ind_wat_limit_4_nosim(aux,hbonds[0],hbonds[1],self.num_waters,self.num_atoms,num_hbs)
#                mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
# 
#                return mss,mss_ind

#    def mss_hbonds_wat_prot(self,definition=1,hbonds=None,verbose=True):
#        
#        mss_funcs.definition_hbs=faux.hbdefinition
# 
#        if hbonds==None:
#            print '# hbonds needed.'
#            return
# 
#        if definition==1:
# 
#            aux=numpy.zeros((self.num_atoms,3),dtype=int,order='F')
#            filt_water=numpy.zeros(self.num_atoms,dtype=bool,order='F')
# 
#            for ii in range(len(self.water)):
#                jo=self.water[ii].O.index
#                jh1=self.water[ii].H1.index
#                jh2=self.water[ii].H2.index
#                aux[jo,0]=ii
#                aux[jo,1]=0
#                aux[jh1,0]=ii
#                aux[jh1,1]=1
#                aux[jh2,0]=ii
#                aux[jh2,1]=2
#                filt_water[jo]=True
#                filt_water[jh1]=True
#                filt_water[jh2]=True
# 
#            #for hbs in hbonds:
#            #    print hbs[1]
#            #print hbonds[0].shape,hbonds[1].shape
#            aux=numpy.array(aux,dtype=int,order='F')
# 
#            if type(hbonds[0][0][0]) in [numpy.ndarray]:
#                mss_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                mss_ind_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                for jj in range(len(hbonds)):
#                    num_hbs=hbonds[jj][0].shape[0]
#                    mss_ind=mss_funcs.ind_wat_limit_4_nosim_prot(aux,hbonds[jj][0],hbonds[jj][1],self.num_waters,self.num_atoms,num_hbs)
#                    mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                    mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
#                    mss_tot[jj,:,:]=mss[:,:]
#                    mss_ind_tot[jj,:,:]=mss_ind[:,:]
#                return mss_tot,mss_ind_tot
# 
#            else:
#                num_hbs=hbonds[0].shape[0]
#                mss_ind=mss_funcs.ind_wat_limit_4_nosim_prot(aux,filt_water,hbonds[0],hbonds[1],self.num_waters,self.num_atoms,num_hbs)
#                mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
# 
#                return mss,mss_ind


#    def mss_hbonds_wation(self,definition=1,hbonds=None,bonds=None,tipo=1,verbose=True):
#        
#        mss_funcs.definition_hbs=faux.hbdefinition
# 
#        if hbonds==None:
#            print '# hbonds needed.'
#            return
# 
#        if definition==1:
# 
#            aux=numpy.zeros((self.num_atoms,3),dtype=int,order='F')
# 
#            for ii in range(len(self.water)):
#                jo=self.water[ii].O.index
#                jh1=self.water[ii].H1.index
#                jh2=self.water[ii].H2.index
#                aux[jo,0]=ii
#                aux[jo,1]=0
#                aux[jh1,0]=ii
#                aux[jh1,1]=1
#                aux[jh2,0]=ii
#                aux[jh2,1]=2
# 
# 
#            #for hbs in hbonds:
#            #    print hbs[1]
#            #print hbonds[0].shape,hbonds[1].shape
#            aux=numpy.array(aux,dtype=int,order='F')
#            
#            if bonds==None:
#                if type(hbonds[0][0][0]) in [numpy.ndarray]:
#                    mss_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                    mss_ind_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                    for jj in range(len(hbonds)):
#                        num_hbs=hbonds[jj][0].shape[0]
#                        mss_ind=mss_funcs.ind_wat_limit_4_nosim(aux,hbonds[jj][0],hbonds[jj][1],self.num_waters,self.num_atoms,num_hbs)
#                        mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                        mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
#                        mss_tot[jj,:,:]=mss[:,:]
#                        mss_ind_tot[jj,:,:]=mss_ind[:,:]
#                    return mss_tot,mss_ind_tot
#                else:
#                    num_hbs=hbonds[0].shape[0]
#                    mss_ind=mss_funcs.ind_wat_limit_4_nosim(aux,hbonds[0],hbonds[1],self.num_waters,self.num_atoms,num_hbs)
#                    mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                    mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
#                    return mss,mss_ind
#            else:
#                num_bonds=len(bonds)
#                if type(hbonds[0][0][0]) in [numpy.ndarray]:
#                    mss_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                    mss_ind_tot=numpy.empty((len(hbonds),self.num_waters,17),dtype=int,order='F')
#                    for jj in range(len(hbonds)):
#                        num_hbs=hbonds[jj][0].shape[0]
#                        mss_ind=mss_funcs.ind_wat_limit_4_nosim(aux,hbonds[jj][0],hbonds[jj][1],self.num_waters,self.num_atoms,num_hbs)
#                        mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                        mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
#                        mss_tot[jj,:,:]=mss[:,:]
#                        mss_funcs.addbonds(tipo,mss,mss_ind,bonds,self.num_waters,num_bonds)
#                        mss_ind_tot[jj,:,:]=mss_ind[:,:]
#                    return mss_tot,mss_ind_tot
#                else:
#                    num_hbs=hbonds[0].shape[0]
#                    mss_ind=mss_funcs.ind_wat_limit_4_nosim(aux,hbonds[0],hbonds[1],self.num_waters,self.num_atoms,num_hbs)
#                    mss=mss_funcs.remove_index_mol(mss_ind,self.num_waters)
#                    mss_funcs.remove_permutations_limit_4_nosim(mss,mss_ind,self.num_waters)
#                    mss_funcs.addbonds(tipo,mss,mss_ind,bonds,self.num_waters,num_bonds)
#                    return mss,mss_ind



#    def mss_hbonds (self,definition=None,set_A=None,set_B=None,acc_don_A=None,acc_don_B=None,traj=0,frame=0,sk_param=0.00850,roh_param=2.3000,roo_param=3.50,angooh_param=30.0,optimize=False,pbc=True,verbose=False):
# 
#        opt_effic=0
#        opt_diff_syst=0
#        opt_diff_set=1
#        opt_pbc=0
#        if pbc:
#            opt_pbc=1
# 
# 
#        if acc_don_A==None and acc_don_B==None:
#            if set_A==None:
#                print 'set_A and/or set_B needed'
#                return
#            else:
#                acc_don_A=self.selection_hbonds(setA=set_A,verbose=False)
#                if set_B==None:
#                    acc_don_B=acc_don_A
#                    opt_diff_set=0
#                else:
#                    acc_don_B=self.selection_hbonds(setA=set_B,verbose=False)
#        else:
#            if acc_don_B==None:
#                acc_don_B=acc_don_A
#                opt_diff_set=0
# 
# 
#        nA_acc       = acc_don_A[0].shape[0]
#        nA_acc_sH    = acc_don_A[1].shape[0] # Just for water and skinner, topological, etc...
#        nA_acc_H     = acc_don_A[2].shape[0] # Just for water and skinner, topological, etc...
#        nA_don       = acc_don_A[3].shape[0]
#        nA_don_sH    = acc_don_A[4].shape[0]
#        nA_don_H     = acc_don_A[5].shape[0]
#        allwat_A     = acc_don_A[6]
#        nB_acc       = acc_don_B[0].shape[0]
#        nB_acc_sH    = acc_don_B[1].shape[0] # Just for water and skinner, topological, etc...
#        nB_acc_H     = acc_don_B[2].shape[0] # Just for water and skinner, topological, etc...
#        nB_don       = acc_don_B[3].shape[0]
#        nB_don_sH    = acc_don_B[4].shape[0]
#        nB_don_H     = acc_don_B[5].shape[0]
#        allwat_B     = acc_don_B[6]
# 
# 
#        num_frames=__length_frame_opt__(self,traj,frame)
# 
#        faux.hbdefinition=hbonds_type(definition,verbose=False)
#        if faux.hbdefinition == 0 : 
#            return
#        
#        elif faux.hbdefinition == 1 : 
#            faux.sk_param=sk_param
#            if not (allwat_A and allwat_B):
#                print '# This type of hbond only works for water molecules.'
#            print 'Not implemented yet'
#            pass
# 
#        elif faux.hbdefinition == 2 : 
#            faux.roh2_param= roh_param**2
#            print 'Not implemented yet'
#            pass
# 
#        elif faux.hbdefinition == 3 : # ROO_ANG
#            faux.roo2_param, faux.cos_angooh_param= roo_param**2, numpy.cos(numpy.radians(angooh_param))
# 
#            if optimize:
#                gg=0
#                hbout=[]
#                for iframe in __read_frame_opt__(self,traj,frame):
#                    if (gg==0): 
#                        self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe)
#                    else:
#                        self.verlet_list_grid_ns(r1=roo_param,r2=roo_param,rcell=roo_param,iframe=iframe,update=True)
# 
#                    faux.get_hbonds_roo_ang_ns_list( opt_diff_set, opt_pbc, \
#                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
#                                           iframe.coors,iframe.box,iframe.orthogonal, \
#                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
#                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
#                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
#                                           self.num_atoms)
# 
#                    hbout.append([faux.hbs_out,faux.hbs_vals_out])
#                    gg+=1
#            else:
#                hbout=[]
#                gg=0
#                for iframe in __read_frame_opt__(self,traj,frame):
#                    faux.get_hbonds_roo_ang( opt_diff_set, opt_pbc, \
#                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
#                                                        iframe.coors,iframe.box,iframe.orthogonal, \
#                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
#                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
#                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
#                                                        self.num_atoms)
#                    hbout.append([faux.hbs_out,faux.hbs_vals_out])
#                    gg+=1
# 
#            if gg==1:
#                return hbout[0]
#            else:
#                return hbout
# 
#        elif faux.hbdefinition == 4 : 
#            if not (allwat_A and allwat_B):
#                print '# This type of hbond only works for water molecules.'
#            print 'Not implemented yet'
#            pass
# 
#        elif faux.hbdefinition == 5 : 
#            if not (allwat_A and allwat_B):
#                print '# This type of hbond only works for water molecules.'
#            print 'Not implemented yet'
#            pass
# 
#        elif faux.hbdefinition == 6 : 
#            faux.cos_angooh_param= numpy.cos(numpy.radians(angooh_param))
#            if not (allwat_A and allwat_B):
#                print '# This type of hbond only works for water molecules.'
#            print 'Not implemented yet'
#            pass
# 
#        elif faux.hbdefinition == 7 : 
#            if not (allwat_A and allwat_B):
#                print '# This type of hbond only works for water molecules.'
#            print 'Not implemented yet'
#            pass
# 
#            #if verbose:
#            #    print 'Listo'
#        pass


#    def hbonds2 (self,definition=None,set_A=None,set_B=None,acc_don_A=None,acc_don_B=None,traj=0,frame=0,sk_param=0.00850,roh_param=2.3000,roo_param=3.5,angooh_param=30.0,optimize=False,pbc=True,verbose=False):
# 
#        opt_effic=0
#        opt_diff_syst=0
#        opt_diff_set=1
#        opt_pbc=0
#        if pbc:
#            opt_pbc=1
# 
#        if acc_don_A==None and acc_don_B==None:
#            if set_A==None:
#                print 'set_A and/or set_B needed'
#                return
#            else:
#                acc_don_A=self.selection_hbonds(setA=set_A,verbose=False)
#                if set_B==None:
#                    acc_don_B=acc_don_A
#                    opt_diff_set=0
#                else:
#                    acc_don_B=self.selection_hbonds(setA=set_B,verbose=False)
#        else:
#            if acc_don_B==None:
#                acc_don_B=acc_don_A
#                opt_diff_set=0
# 
# 
#        nA_acc       = acc_don_A[0].shape[0]
#        nA_acc_sH    = acc_don_A[1].shape[0] # Just for water and skinner, topological, etc...
#        nA_acc_H     = acc_don_A[2].shape[0] # Just for water and skinner, topological, etc...
#        nA_don       = acc_don_A[3].shape[0]
#        nA_don_sH    = acc_don_A[4].shape[0]
#        nA_don_H     = acc_don_A[5].shape[0]
#        allwat_A     = acc_don_A[6]
#        nB_acc       = acc_don_B[0].shape[0]
#        nB_acc_sH    = acc_don_B[1].shape[0] # Just for water and skinner, topological, etc...
#        nB_acc_H     = acc_don_B[2].shape[0] # Just for water and skinner, topological, etc...
#        nB_don       = acc_don_B[3].shape[0]
#        nB_don_sH    = acc_don_B[4].shape[0]
#        nB_don_H     = acc_don_B[5].shape[0]
#        allwat_B     = acc_don_B[6]
# 
#        natomA=self.num_atoms
#        natomB=self.num_atoms
#        num_frames=__length_frame_opt__(self,traj,frame)
# 
#        faux.hbdefinition=hbonds_type(definition,verbose=False)
# 
#        if faux.hbdefinition == 3 : # ROO_ANG
#            faux.roo2_param, faux.cos_angooh_param= roo_param**2, numpy.cos(numpy.radians(angooh_param))
#            
#            gg=0
#            for iframe in __read_frame_opt__(self,traj,frame):
#                if (gg==0): 
#                    self.verlet_list_grid_ns(r1=3.5,r2=3.5,rcell=3.5,iframe=iframe)
#                else:
#                    self.verlet_list_grid_ns(r1=3.5,r2=3.5,rcell=3.5,iframe=iframe,update=True)
#                    faux.get_hbonds_roo_ang_ns_list(opt_effic, opt_diff_syst, opt_diff_set, opt_pbc, \
#                                                       acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
#                                                       iframe.coors,iframe.box,iframe.orthogonal, \
#                                                       acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
#                                                       iframe.coors,nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
#                                                       nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
#                                                       natomA,natomB)
#                gg+=1
# 
# 
#        pass


    def water_bisector (self,water='ALL',traj=0,frame=0,pbc=True):

        opt_pbc=0
        if pbc:
            opt_pbc=1

        if water in ['ALL','All','all']:
            wat_ind=[ii for ii in range(self.num_waters)]
        elif type(water) in [int]:
            wat_ind=[water]
        elif type(water) in [list,tuple]:
            wat_ind=water
        else:
            print '# List of water molecules in input variable "water" needed.'
            return

        list_atoms=numpy.empty((len(wat_ind),3),dtype=int,order='F')
        for ii in wat_ind:
            list_atoms[ii,0]=self.water[ii].O.index
            list_atoms[ii,1]=self.water[ii].H1.index
            list_atoms[ii,2]=self.water[ii].H2.index

        ddd=[]
        for iframe in __read_frame_opt__(self,traj,frame):
            dd=faux.water_bisector(opt_pbc,list_atoms,iframe.coors,iframe.box,iframe.orthogonal,len(list_atoms),self.num_atoms)
            ddd.append(dd)
            
    def water_angle_bisector_atom (self,water='ALL',atoms=None,traj=0,frame=0,pbc=True):

        opt_pbc=0
        if pbc:
            opt_pbc=1

        if water in ['ALL','All','all']:
            wat_ind=[ii for ii in range(self.num_waters)]
        elif type(water) in [int]:
            wat_ind=[water]
        elif type(water) in [list,tuple]:
            wat_ind=water
        else:
            print '# List of water molecules in input variable "water" needed.'
            return

        list_atoms=numpy.empty((len(wat_ind),3),dtype=int,order='F')
        for ii in wat_ind:
            list_atoms[ii,0]=self.water[ii].O.index
            list_atoms[ii,1]=self.water[ii].H1.index
            list_atoms[ii,2]=self.water[ii].H2.index

        if type(atoms) in [int]:
            atoms=[atoms]
            natoms=1
        elif type(atoms) in [tuple,list]:
            natoms=len(atoms)
        else:
            print '# List of atoms in input variable "atoms" needed.'
            return

        num_frames=__length_frame_opt__(self,traj,frame)
        angls=numpy.empty(shape=(num_frames,len(wat_ind),natoms),dtype=float,order='Fortran')
        
        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            angls[num_frames,:,:]=faux.water_angle_bisector_atom(opt_pbc,atoms,list_atoms,iframe.coors,iframe.box,iframe.orthogonal,len(list_atoms),natoms,self.num_atoms)
            num_frames+=1

        if num_frames==1:
            return angls[0,:,:]
        else:
            return angls


    #def contact_map (self,setA=None,setB=None,dist=None,traj=0,frame=0,pbc=True):
    # 
    #    setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,setA,None,setB)
    #    num_frames=__length_frame_opt__(self,traj,frame)
    # 
    #    contact_map=numpy.empty(shape=(n_A,n_B,num_frames),dtype=int,order='Fortran')
    #    num_frames=0
    #    for iframe in __read_frame_opt__(self,traj,frame):
    #        contact_map[:,:][num_frames]=f.contact_map(diff_system,pbc,dist,setA,iframe.coors,iframe.box,setB,iframe.coors,n_A,n_B,natoms_A,natoms_B)
    #        num_frames+=1
    # 
    #    if num_frames==1:
    #        return contact_map[:,:][0]
    #    else:
    #        return contact_map
    # 
    #    #def plot_contact_map(contact_map):
    #    #    
    #    #    pylab.gray()
    #    #    pylab.imshow(contact_map==False,origin='lower',interpolation=None) # Would be better not to interpolate
    #    #    #pylab.matshow(contact_map==False)
    #    #    return pylab.show()


    #def displ_vector(self,set_reference=None):
    # 
    #    self.d_vector=set_reference.frame[0].coors - self.frame[0].coors


#### END CLASSES
#######################################################
#######################################################



#######################################################
#######################################################
#### EXTERNAL OBJECTS AND FUNCTIONS
    
####
#### Functions
####

def __read_frame_opt__ (syst=None,traj=0,frame=None):

    #Options: integer,list,tuple,'ALL'

    if frame in ['ALL','All','all']:
        return syst.traj[traj].frame
    elif type(frame) in [numpy.int32,int]:
        return [syst.traj[traj].frame[ii] for ii in [frame]]
    elif type(frame) in [list,tuple]:
        return [syst.traj[traj].frame[ii] for ii in frame]

def __length_frame_opt__ (syst=None,traj=0,frame=None):

    #Options: integer,list,tuple,'ALL'

    if frame in ['ALL','All','all']:
        return syst.traj[traj].num_frames
    elif type(frame) in [numpy.int32,int]:
        return 1
    elif type(frame) in [list,tuple]:
        return len(frame)

def __read_sets_opt__(systA=None,setA=None,systB=None,setB=None):

    if setA==None:
        print '# SetA needed.'
        pass

    setA,nlist_a,nsys_a=__read_set_opt__(systA,setA)

    diff_syst=1
    diff_set=1


    if systB==None:

        systB=systA
        diff_syst=0

        if setB in [None]:
            setB=setA
            diff_set=0

    setB,nlist_b,nsys_b=__read_set_opt__(systB,setB)

    return setA,nlist_a,nsys_a,setB,nlist_b,nsys_b,diff_syst,diff_set

          
def __read_set_opt__(systA=None,setA=None):

    if setA==None:
        print '# SetA needed.'
        pass

    nsys_a=systA.num_atoms

    if setA in ['ALL','All','all']:
        setA=[ii for ii in range(systA.num_atoms)]
        nlist_a=systA.num_atoms
    elif type(setA) in [numpy.int32,numpy.int64,int]:
        setA=[setA]
        nlist_a=1
    elif type(setA) in [list,tuple]:
        nlist_a=len(setA)
    else:
        setA=systA.selection(setA)
        nlist_a=len(setA)

    return setA,nlist_a,nsys_a


          


#def min_distance(system,set_a,set_b=None,pbc=True,type_b='atoms'):
# 
#    if set_b==None:
# 
#        ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,True,system.frame[0].coors,system.frame[0].box,set_a,set_a,system.num_atoms,len(set_a),len(set_a))
# 
#    else:
#        if type_b=='atoms':
#            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,False,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),len(set_b))
#        elif type_b=='vectors':
#            l_vects=shape(set_b)
#            if len(l_vects)==1:
#                set_b=[set_b]
#                l_vects=shape(set_b)
#            numpy.array(set_b,order='Fortran')
#            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms_ref(pbc,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),l_vects[0])
#    
#    return ind_a1,ind_a2,min_dist
# 
# 
#def xtc2bin(xtc_name,bin_name):
#    command=home_path+'xtc2bin %s %s'%(xtc_name,bin_name)
# 
#    if path.exists(bin_name):
#        command2='mv %s %s#'%(bin_name,bin_name)
#        print 'file',bin_name,' was moved to ', bin_name+'#'
#        system(command2)
# 
#    system(command)
# 
# 
#def dot_product_3d(vect1,vect2):
# 
#    return f.aux_funcs_general.proj3d(vect1,vect2,len(vect1))
# 
#def isothermal_compressibility(system,Temp,input_file,frame=None,begin=None,end=None):
# 
#    frame=[ii for ii in range(begin,end+1)]
# 
#    V2a=0.0
#    Va=0.0
#    Kt=0.0
#    for ii in frame :
#        system.delete_coors()
#        system.load_coors (input_file,frame=ii)
#        xx=0.0
#        xx=system.frame[0].box[0][0]*system.frame[0].box[1][1]*system.frame[0].box[2][2]
#        Va+=xx
#        V2a+=xx**2
#    
#    Va=Va/(len(frame)*1.0)
#    V2a=V2a/(len(frame)*1.0)
#    Kt=(V2a-Va**2)/(Temp*Va)
#    
#    return Kt*0.10,'(nm/Kb)'

def hbonds_type(option=None,verbose=True):

    hbs_type={}
    hbs_info={}
    hbs_type['Skinner']=1; hbs_info['Skinner']='R.Kumar, J.R. Schmidt and J.L. Skinner. J. Chem. Phys. 126, 204107 (2007)' 
    hbs_type['R(o,h)']=2;  hbs_info['R(o,h)']='V. J. Buch. J. Chem. Phys. 96, 3814-3823 (1992)'
    hbs_type['R(o,o)-Ang(o,o,h)']=3; hbs_info['R(o,o)-Ang(o,o,h)']='A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)'
    hbs_type['R(o,o)-Ang(o,h,o)']=8; hbs_info['R(o,o)-Ang(o,h,o)']=' '
    hbs_type['R(o,h)-Ang(o,h,o)']=9; hbs_info['R(o,h)-Ang(o,h,o)']=' '
    hbs_type['Donor-Acceptor-Number']=4; hbs_info['Donor-Acceptor-Number']='A. D. Hammerich, V. J. Buch. J. Chem. Phys. 128, 111101 (2008)'
    hbs_type['Topological']=5; hbs_info['Topological']='R. H. Henchman and S. J. Irudayam. J. Phys. Chem. B. 114, 16792-16810 (2010)'
    hbs_type['Donor-Number-Ang(o,o,h)']=6; hbs_info['Donor-Number-Ang(o,o,h)']='J. D. Smith, C. D. Cappa, et al. Proc. Natl. Acad. Sci. U.S.A. 102, 14171 (2005).'
    hbs_type['Nearest-Neighbour']=7; hbs_info['Nearest-Neighbour']='This is not a hydrogen bond definition but just a topological characterization.'



    #if verbose:
    #    if option not in hbs_type.keys():
    #        for ii in hbs_type.keys():
    #            if len(ii)<=12: tab='\t\t\t'
    #            if 12<len(ii)<=18: tab='\t\t'
    #            if 18<len(ii): tab='\t'
    #            print '  ',ii,tab+'[',hbs_info[ii],']'
    #    return

    if verbose:
        if option not in hbs_type.keys():
            for ii in ['R(o,h)','R(o,o)-Ang(o,o,h)','R(o,o)-Ang(o,h,o)','R(o,h)-Ang(o,h,o)',
                       'Skinner','Donor-Acceptor-Number','Donor-Number-Ang(o,o,h)','Topological',
                       'Nearest-Neighbour']:
                if len(ii)<=12: tab='\t\t\t'
                if 12<len(ii)<=18: tab='\t\t'
                if 18<len(ii): tab='\t'
                print '  ',ii,tab+'[',hbs_info[ii],']'
        return


    if option != None :
        if option not in hbs_type.keys():
            print option, ': Hbond type not defined.'
            print 'List of definitions:'
            for ii in hbs_type.keys():
                if len(ii)<=12: tab='\t\t\t'
                if 12<len(ii)<=18: tab='\t\t'
                if 18<len(ii): tab='\t'
                print '  ',ii,tab+'[',hbs_info[ii],']'
            return 0
        return hbs_type[option]


######################################################

def add_user_topol(file_topol=None,verbose=False):

    if file_topol.endswith('.py'):
        file_topol.replace('.py','')

    tp.add_topol(tp,file_topol,verbose)

######################################################


def selection_covalent_chains(system=None,chain=None,select='ALL'):
 
    if chain=='ALL':

        nn=len(system.atom)
        aux_dict={}
        aux_list=numpy.zeros((nn),dtype=int)
        aux_filt=numpy.zeros((nn),dtype=bool)

        gg=0
        for ii in range(nn):
            if aux_filt[ii]==False:
                gg=gg+1
                aux_dict[gg]=[]
                aux_dict[gg].append(ii)
                aux_filt[ii]=True
                aux_list[ii]=gg
                ggg=gg
            else:
                ggg=aux_list[ii]
            for jj in system.atom[ii].covalent_bonds:
                if aux_filt[jj] and aux_list[jj]!=ggg:
                    kk=aux_list[jj]
                    ll=aux_dict.pop(kk)
                    aux_dict[ggg].extend(ll)
                    for mm in ll:
                        aux_list[mm]=ggg
                elif aux_filt[jj]==False:
                    aux_filt[jj]=True
                    aux_list[jj]=ggg
                    aux_dict[ggg].append(jj)

        salida=[]
        for aa,bb in aux_dict.iteritems():
            salida.append(bb)

        del(aux_dict,aux_list,aux_filt)

        return salida
        pass

    else:

        setC,nlist_C,nsys_C=__read_set_opt__(system,select)
        
        for ii in range(len(chain)):
            if type(chain[ii]) is not list:
                chain[ii]=[chain[ii]]

        aux_dict={}
        for ii in chain:
            for jj in ii:
                aux_dict[jj]=[]
 
        aux_list=aux_dict.keys()

        for ii in setC:
            if system.atom[ii].name in aux_list:
                aux_dict[system.atom[ii].name].append(ii)

        aux_list=[]
        for ii in range(len(chain)):
            aux_list.append([])
            for jj in chain[ii]:
                aux_list[ii].extend(aux_dict[jj])
            aux_list[ii].sort()
 
        bonds=[]
        for ii in range(len(aux_list)-1):
            bonds.append({})
            for jj in aux_list[ii]:
                bonds[ii][jj]=[]
                for kk in system.atom[jj].covalent_bonds:
                    if kk in aux_list[ii+1]:
                        bonds[ii][jj].append(kk)

        salida=[[ii] for ii in aux_list[0]]
        for kk in range(len(aux_list)-1):
            salida2=[]
            for ii in range(len(salida)):
                desde=salida[ii][-1]
                for hacia in bonds[kk][desde]:
                    salida2.append(salida[ii]+[hacia])
            salida=salida2

        del(salida2,bonds,aux_list,aux_dict,setC,nlist_C,nsys_C)

        return salida


def selection(system=None,condition=None,traj=0,frame='ALL',pbc=True):

    icondition=' '+condition+' '

    # attributes syntaxis:

    dict_selects={
        'backbone':  '(atom.resid.type Protein and atom.name N CA C O CH3)',
        'sidechain': '(atom.resid.type Protein and not atom.name N CA C O H1 H2)',
        'protein':   '(atom.resid.type Protein)',
        'water':     '(atom.resid.type Water)',
        'ion':       '(atom.resid.type Ion)',
        'lipid':     '(atom.resid.type Lipid)'
        }

    for ii,jj in dict_selects.iteritems():
        icondition=icondition.replace(ii,jj)

    ### First block.
    icondition=icondition.replace(',',' ')
    icondition=icondition.replace('(',' ( ')
    icondition=icondition.replace(')',' ) ')
    icondition=icondition.replace('[',' ( ')
    icondition=icondition.replace(']',' ) ')
    icondition=icondition.replace(' AND ',' and ')
    icondition=icondition.replace(' And ',' and ')
    icondition=icondition.replace(' OR ',' or ')
    icondition=icondition.replace(' Or ',' or ')
    icondition=icondition.replace(' NOT ',' not ')
    icondition=icondition.replace(' Not ',' not ')
    icondition=icondition.replace(' Within ',' within ')
    icondition=icondition.replace(' WITHIN ',' within ')
    icondition=icondition.replace(' Of ',' of ')
    icondition=icondition.replace(' OF ',' of ')
    icondition=icondition.replace(' in ',' ')
    icondition=icondition.replace(' In ',' ')
    icondition=icondition.replace(' IN ',' ')
    icondition=icondition.replace(' : ',':')
    ### Second block.
    icondition=icondition.replace(' chain.',' atom.chain.')
    icondition=icondition.replace(' resid.',' atom.resid.')
    icondition=icondition.replace('  ',' ')

    # logic syntaxis

    icondition=icondition.split()
    aux_cond=[True for ii in range(len(icondition))]
    ocondition=[]

    for ii in range(len(icondition)):
        if icondition[ii] in ['(',')','and','or','not','within','of']:
            aux_cond[ii]=False
        if icondition[ii].startswith('atom'):
            aux_cond[ii]=False

    aux_cond.append(False)

    for ii in range(len(icondition)):
        part=icondition[ii]
        if part in ['(',')','and','or','not','of']:
            ocondition.append(part)
        elif part in ['within']:
            part2=icondition[ii+1]
            aux_cond[ii+1]=False
            ocondition.append(part)
            ocondition.append(float(part2))
        elif part.startswith('atom'):
            ocondition.append(part)
            if aux_cond[ii+1]:
                ocondition.append('in'); ocondition.append('[')
                jj=ii
                while True:
                    jj+=1
                    if not aux_cond[jj]:
                        ocondition.append(']')
                        break
                    else:
                        part2=icondition[jj]
                        if ':' in part2:
                            cc=part2.split(':')
                            from_cc=int(cc[0])
                            to_cc  =int(cc[1])
                            try:
                                inc_cc=int(cc[2])
                            except:
                                inc_cc=1
                            part2=str(range(from_cc,to_cc,inc_cc)).replace('[','')
                            part2=part2.replace(']','')
                            ocondition.append(part2)
                        else:
                            try:
                                kk=float(part2)
                                ocondition.append(part2+',')
                            except:
                                ocondition.append("'"+part2+"',")
                        aux_cond[jj]=False
        else:
            if aux_cond[ii]:
                ocondition.append(part)
                aux_cond[ii]=False

    # Solving the 'withins':
    for ii in range(len(ocondition)):
        if 'within' == ocondition[ii]:
            cutoff=float(ocondition[ii+1])
            if ')' not in ocondition[0:ii]:
                sel1=[]
                cond1=' '.join(ocondition[0:ii])
                for atom in system.atom:
                    if eval(cond1):
                        sel1.append(atom.index)
            else:
                sel1=[]
                cond1=' '.join(ocondition[0:ii])
                for atom in system.atom:
                    if eval(cond1):
                        sel1.append(atom.index)
            if '(' not in ocondition[(ii+3):]:
                sel2=[]
                cond2=' '.join(ocondition[(ii+3):])
                for atom in system.atom:
                    if eval(cond2):
                        sel2.append(atom.index)
            else:
                sel2=[]
                cond2=' '.join(ocondition[(ii+3):])
                for atom in system.atom:
                    if eval(cond2):
                        sel2.append(atom.index)
            dists_sels=system.distance(sel1,sel2,traj=traj,frame=frame,pbc=pbc)
            list_sel1=[]
            for gg in range(len(dists_sels)):
                list_sel_frame=[]
                for jj in range(len(sel1)):
                    if faux.within(dists_sels[gg][jj,:],cutoff,len(sel2)):
                        list_sel_frame.append(sel1[jj])
                list_sel1.append(list_sel_frame)
            #'(atom.resid.type Water and atom.type O) within 3.0 of atom.resid.type Protein'
            #atom.name OW within 3.0 of [3,4,5]
            #atom.name OW within 3.0 of sel1
            #atom.name OW within 3.0 of atom.name HW1
            if len(list_sel1)==1:
                return list_sel1[0]
            else:
                return list_sel1

    # Applying selection
    list_select=[]
    condition=' '.join(ocondition)    

    for atom in system.atom:
        if eval(condition):
            list_select.append(atom.index)
     
     
    return list_select










