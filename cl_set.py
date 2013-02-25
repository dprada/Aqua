####################################
# GENERAL COMMENTS

NAME_VERSION="Pynoramix 0.1"

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
from pyn_cl_coors import *
import top_par as tp
import pyn_fort_general as faux
import pyn_math as pyn_math

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
        self.polarizability=False       # True of falsel

####
#### Class residue (set of atoms)
####


class cl_residue(labels_set):           # Attributes of a residue (inherits from label_set)

    def __init__( self ):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms, .__int_type

        self.chain=labels_unit()         # Chain which this residue belongs to.
        self.__int_dict_atoms__={}         # Dictionary with internal atom names and indexes

        pass


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


class molecule(labels_set):               # The suptra-estructure: System (waters+cofactors+proteins...)

    
    def __init__(self,input_file=None,download=None,coors=True,verbose=True,with_bonds=True,missing_atoms=True):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        # > Instantation options:
        self.file_topol=input_file       # pdb,gro...
        self.file_hbonds=''             # still not useful -do not remove-
        self.file_mss=''                # still not useful -do not remove-
        self.file_shell=''              # still not useful -do not remove-
        #self.selection=False            # is the set a selection (a mirror -pointer- of a parent set)
        
        # > Topological properties
        self.atom=[]                    # list of atoms    (objects: cl_unit)
        self.resid=[]                   # list of residues (objects: molecule)
        self.chain=[]                   # list of chains   (objects: molecule)
        self.chains=[]                  # list of chain names (strings)
        self.ion=[]                     # list of ions (objects: molecule)
        self.water=[]                   # list of waters   (objects: cl_water)
        self.water_model=None           # water model
        self.parent=None                # Parent set if it comes from selection (object: labels_parent)

        # > Physical and Chemical properties
        self.acceptors=[]               # list of acceptors
        self.donors=[[],[]]             # list of [donors,donor_hydrogen]
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

        # A SET CAN BE BUILT FROM A FILE OR FROM A SELECTION

        # IF IT COMES FROM A FILE, DOES IT NEED TO BE DOWNLOADED?

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
                    temp_residue.type=tp.residue_type[atom.resid.name]
                    temp_residue.__int_name__=tp.residue[atom.resid.name]
                    xxx_resid_type=temp_residue.type
                    xxx_resid__int_name__=temp_residue.__int_name__
                    temp_residue.chain.name=atom.chain.name
                    temp_residue.chain.index=kk
                    temp_residue.__int_dict_atoms__={}
                    self.resid.append(temp_residue)
                ii+=1                                      #### Atom
                atom.index=ii                   
                atom.resid.index=jj             
                atom.chain.index=kk
                atom.hbonds=[]
                atom.__int_name__=tp.atom[atom.name]
                atom.type=tp.atom_type[atom.__int_name__]
                atom.resid.__int_name__=xxx_resid__int_name__
                atom.resid.type=xxx_resid_type
                self.resid[jj].list_atoms.append(ii)
                self.resid[jj].__int_dict_atoms__[atom.__int_name__]=ii
                self.chain[kk].list_atoms.append(ii)
                if atom.type=='H':
                        without_hs=False

            ### Setting up the subsets.

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
             
                if residue.type=='Ion':        ### Ions
                    temp_residue=cl_residue()
                    temp_residue.list_atoms=residue.list_atoms
                    temp_residue.index=residue.index
                    temp_residue.pdb_index=residue.pdb_index
                    temp_residue.name=residue.name
                    self.ion.append(temp_residue)


            ### Setting up the local attributes

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

                    # Peptide bond [C(i)->N(i+1)]
                    try:
                        next_residue=self.resid[residue.index+1]
                        if residue.type=='Protein' and next_residue.type=='Protein':
                            if residue.chain.name==next_residue.chain.name :
                                aa=residue.__int_dict_atoms__['atC']
                                bb=next_residue.__int_dict_atoms__['atN']
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
                if atom.__int_name__ in tp.donors_with_exceptions:
                    try:
                        exception=donors_exception[atom.__int_name__][atom.resid.__int_name__]
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
                if atom.__int_name__ in tp.acceptors_with_exceptions:
                    try:
                        exception=acceptors_exception[atom.__int_name__][tp.resid.__int_name__]
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
                            self.donors[0].append(atom.index)
                            self.donors[1].append(ii)

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
            self.list_atoms=[ii for ii in range(self.num_atoms)]


            ### Loading coordinates
            if coors:
                self.load_traj(self.file_topol,frame='ALL',verbose=False)


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
        print '#','System created from the file',self.file_topol,':'
        print '#',self.num_atoms,' atoms'
        print '#',self.num_residues,' residues'
        print '#',self.num_chains,' chains'
        print '#',self.num_waters,' waters'
        print '#',self.num_ions,' ions'

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

    def write_pdb (self,filename=None):
        
        if filename==None:
            print 'Enter filename: '
            print '      foo.write_pdb("foo.pdb")'
        else:
            if not filename.endswith('.pdb'): filename+='.pdb'
            if path.exists(filename): 
                print '# The file '+filename+' already exists.'
                return
     
            fff=open(filename,'w')
     
            a='HEADER    '+'> CREATED BY PYNORAMIX '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+' <\n'
            fff.write(str(a))
     
            for ii in self.pdb_header: fff.write(str(ii))
            
            for ii in self.pdb_ss:     fff.write(str(ii))
     
            dct_aux={'ATOM': 'ATOM  ', 'HETATM': 'HETATM'}
            
            new_index=0
            for ii in range(self.num_atoms):
                new_index+=1
                a=dct_aux[self.atom[ii].type_pdb]              # 1-6
                a+="%5d" % (new_index)                         # 7-11
                #a+="%5d" % self.atom[ii].pdb_index            # 7-11
                a+=' '                                         # 12
                a+=' '+"%-3s" % self.atom[ii].name             # 13-16
                a+=' '                                         # 17
                a+="%3s" % self.atom[ii].resid.name            # 18-20
                a+=' '                                         # 21
                a+="%1s" % self.atom[ii].chain.name            # 22
                a+="%4d" % self.atom[ii].resid.pdb_index       # 23-26
                a+=' '                                         # 27
                a+='   '                                       # 28-30
                a+="%8.3f" % float(self.frame[0].coors[ii][0])   # 31-38
                a+="%8.3f" % float(self.frame[0].coors[ii][1])   # 39-46
                a+="%8.3f" % float(self.frame[0].coors[ii][2])   # 47-54
                a+="%6.2f" % self.atom[ii].occup               # 55-60
                a+="%6.2f" % self.atom[ii].bfactor             # 61-66
                a+='          '                                # 67-76
                a+="%2s" % self.atom[ii].elem_symb             # 77-78
                a+="%2s" % self.atom[ii].charge                # 79-80
                a+='\n' 
                fff.write(str(a))         
                if ii<(self.num_atoms-1):
                    if self.atom[ii].type_pdb!=self.atom[ii+1].type_pdb or self.atom[ii].chain.name!=self.atom[ii+1].chain.name :
                        new_index+=1
                        a="TER   "
                        a+="%5d" % (new_index)
                        a+=' '
                        a+='  '
                        a+=' '                                         
                        a+="%3s" % self.atom[ii].resid.name            
                        a+='\n' 
                        fff.write(str(a))
            a='END   '+'\n'
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

    def load_traj (self,file_input=None,frame=None,begin=None,end=None,increment=1,units='frames',verbose=True):

        temp_traj=cl_traj(file_input,frame,begin,end,increment,units,verbose=False)
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

        list_condition=selection(self,condition,traj,frame,pbc)
        return list_condition

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

###############################################################
###############################################################
    # To handle the set

    def translation(self,select='ALL',vector=None,wrap=True):

        pass

    def center(self,select='ALL',center_of=None,mass_weighted=False,traj=0,frame=0,pbc=True,wrap=True):

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

        setcom,nlist_setcom,numsys=__read_set_opt__(self,center_of)
        setmov,nlist_setmov,numsys=__read_set_opt__(self,select)

        for iframe in __read_frame_opt__(self,traj,frame):
            faux.glob.center(pbc_opt,setcom,setmov,iframe.coors,iframe.box,iframe.orthogonal,\
                                 nlist_setcom,nlist_setmov,numsys)
            if wrap: 
                iframe.wrap()
            
        pass

    def rotation(self,select='ALL',rotor=None,wrap=True):

        pass


    def distance(self,setA='ALL',setB=None,traj=0,frame='ALL',pbc=True):
        
        if pbc:
            check_cell=self.traj[traj].frame[0].cell
            if check_cell[0,1]!=90 or check_cell[0,2]!=90 or check_cell[1,2]!=90:
                print '# PBC not implemented for not orthorhombic boxes'
                return
            pbc=1

        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,setA,None,setB)
        num_frames=__length_frame_opt__(self,traj,frame)
        dists=numpy.empty(shape=(num_frames,nlist_A,nlist_B),dtype=float,order='Fortran')

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            dists[num_frames,:,:]=faux.glob.distance(diff_syst,diff_set,pbc,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
            num_frames+=1

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
            min_dists[:,num_frames],ind_atoms_min[:,num_frames],min_image[:,:,num_frames]=faux.glob.distance_images(diff_syst,diff_set,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
            num_frames+=1

        if num_frames==1:
            return min_dists[:,0],ind_atoms_min[:,0],min_image[:,:,0]
        else:
            return min_dists, ind_atoms_min,min_image

    def radius_gyration(self,setA='ALL',traj=0,frame='ALL'):

        setA,nlist_A,nsys_A=__read_set_opt__(self,setA)
        num_frames=__length_frame_opt__(self,traj,frame)
        rgs=numpy.empty(shape=(num_frames),dtype=float,order='Fortran')

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            rgs[num_frames]=faux.glob.radius_gyration(setA,iframe.coors,iframe.box,iframe.orthogonal,nlist_A,nsys_A)
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
            piaxis[:,:,num_frames]=faux.glob.principal_inertia_axis(setA,iframe.coors,iframe.box,iframe.orthogonal,nlist_A,nsys_A)
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
            pgaxis[:,:,num_frames]=faux.glob.principal_geometric_axis(setA,iframe.coors,iframe.box,iframe.orthogonal,nlist_A,nsys_A)
            num_frames+=1

        if num_frames==1:
            return pgaxis[:,:,0]
        else:
            return pgaxis


    def dihedral_angles(self,select='protein',phi=True,psi=True,omega=True,xis=False,traj=0,frame='ALL',begin=0,end=None,legend=False,verbose=False):
 
        setA,nlist_A,nsys_A=__read_set_opt__(self,select)

        list_omega=[]
        list_phi=[]
        list_psi=[]
        list_angs=[]
        key_legend=[]

        if omega:
            list_omega=select_covalent_chain(system=self,select=setA,chain=[['CA','CH3'],'C','N',['CA','CH3']])
            for ii in list_omega:
                list_angs.append(ii)
                if legend:
                    key_legend.append(['Omega '+str(self.atom[ii[2]].resid.pdb_index),ii])
            
        if phi:
            list_phi=select_covalent_chain(system=self,select=setA,chain=['C','N','CA','C'])
            for ii in list_phi:
                list_angs.append(ii)
                if legend:
                    key_legend.append(['Phi '+str(self.atom[ii[2]].resid.pdb_index),ii])
        if psi:
            list_psi=select_covalent_chain(system=self,select=setA,chain=['N','CA','C','N'])
            for ii in list_psi:
                list_angs.append(ii)
                if legend:
                    key_legend.append(['Psi '+str(self.atom[ii[1]].resid.pdb_index),ii])
        if xis:
            print '# Option not implemented yet: xis=True'
            return

        list_angs=numpy.array(list_angs,dtype=int,order='F')
        num_dih_angs=list_angs.shape[0]

        if not end:
            num_frames=__length_frame_opt__(self,traj,frame)
            dih_angs=numpy.empty(shape=(num_frames,num_dih_angs),dtype=float,order='Fortran')
        
            for iframe in __read_frame_opt__(self,traj,frame):
                dih_angs[num_frames,:]=faux.glob.dihedral_angles(iframe.coors,iframe.box,iframe.orthogonal,list_angs,num_dih_angs,nsys_A)

        if end:
            num_frames=end-begin+1
            dih_angs=numpy.empty(shape=(num_frames,num_dih_angs),dtype=float,order='Fortran')
            for ii in range(end):
                self.traj[traj].reload_frame(frame=ii)
                iframe=self.traj[traj].frame[0]
                dih_angs[ii,:]=faux.glob.dihedral_angles(iframe.coors,iframe.box,iframe.orthogonal,list_angs,num_dih_angs,nsys_A)

         
         
        #legend=0
        if num_frames==1:
            if legend:
                return key_legend,dih_angs[0,:]
            else:
                return dih_angs[0,:]
        else:
            if legend:
                return key_legend,dih_angs
            else:
                return dih_angs


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

        xxx=pyn_math.binning(None,bins,segment,None,None)
        rdf_tot=numpy.zeros(shape=(bins),dtype=float,order='Fortran')
        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            dist_frame=faux.glob.distance(1,pbc,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
            rdf_frame=faux.rdf.rdf_frame(dist_frame,frame.box,segment[0],segment[1],bins,nlist_A,nlist_B)
            rdf_tot+=rdf_frame
            num_frames+=1.0

        rdf_tot=rdf_tot/(num_frames*1.0)
        return xxx,rdf_tot


    def neighbs_old(self,setA=None,setB=None,ranking=1,dist=None,traj=0,frame=0,pbc=True):
     
        setA,nlist_A,nsys_A,setB,nlist_B,nsys_B,diff_syst,diff_set=__read_sets_opt__(self,setA,None,setB)
        num_frames=__length_frame_opt__(self,traj,frame)

        if dist==None:
            neighbs=[]
            neighbs=numpy.empty(shape=(nlist_A,ranking,num_frames),dtype=int,order='Fortran')
            num_frames=0
            for iframe in __read_frame_opt__(self,traj,frame):
                neighbs[:,:][num_frames]=faux.glob.neighbs_ranking(diff_syst,diff_set,pbc,ranking,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
                num_frames+=1
            if num_frames==1:
                return neighbs[:,:][0]
            else:
                return neighbs
            
        else:
            neighbs=[]
            sort_opt=0
            if ranking:
                sort_opt=1
            for iframe in __read_frame_opt__(self,traj,frame):
                contact_map,num_neighbs,dist_matrix=faux.glob.neighbs_dist(diff_syst,diff_set,pbc,dist,setA,iframe.coors,iframe.box,iframe.orthogonal,setB,iframe.coors,nlist_A,nlist_B,nsys_A,nsys_B)
                aux_neighbs=[]
                for ii in range(nlist_A):
                    if num_neighbs[ii]:
                        neighbs_A=faux.glob.translate_list(sort_opt,setB,contact_map[ii,:],dist_matrix[ii,:],num_neighbs[ii],nlist_B)
                        aux_neighbs.append(neighbs_A)
                    else:
                        aux_neighbs.append([])
                neighbs.append(aux_neighbs)
            if num_frames==1:
                return neighbs[0]
            else:
                return neighbs



    def verlet_list_ns (self,r1=3.5,r2=6.0,traj=0,frame=0,pbc=True,update=False,verbose=False):

        print 'Not implemented as independent function yet'
        pass

    def grid_ns_list (self,rcell=6.0,setA='ALL',setB=None,traj=0,frame=0,pbc=True,verbose=False):
        
        print 'Not implemented as independent function yet'
        pass

    def make_cell_grid_ns(self,rcell=7.0,rcut=3.5,traj=0,frame=0):
        for iframe in __read_frame_opt__(self,traj,frame):
            faux.glob.make_cell_ns(rcell,rcut,iframe.box,self.num_atoms)

    def verlet_list_grid_ns(self,r1=3.5,r2=7.0,rcell=7.0,traj=0,frame=0,iframe=None,pbc=True,update=False,verbose=False):

        pbc_opt=0
        if pbc:
            pbc_opt=1

        if iframe!=None:
            if update:
                faux.glob.update_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
            else:
                #print 'aqui si'
                faux.glob.make_cell_ns(rcell,r2,iframe.box,self.num_atoms)
                #print 'aqui tambien'
                faux.glob.make_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
                #print 'aqui sale'
        else:
            if update:
                for iframe in __read_frame_opt__(self,traj,frame):
                    faux.glob.update_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
            else:
                for iframe in __read_frame_opt__(self,traj,frame):
                    faux.glob.make_cell_ns(rcell,r2,iframe.box,self.num_atoms)
                    faux.glob.make_verlet_list_grid_ns(r1,r2,pbc_opt,iframe.coors,iframe.box,iframe.volume,iframe.orthogonal,self.num_atoms)
        


    def hbonds (self,definition=None,set_A=None,set_B=None,acc_don_A=None,acc_don_B=None,traj=0,frame=0,sk_param=0.00850,roh_param=2.3000,roo_param=3.50,angooh_param=30.0,optimize=False,pbc=True,infile=False,verbose=False):

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

        faux.hbonds.definition=hbonds_type(definition,verbose=False)
        

        if faux.hbonds.definition == 0 : 
            return
        
        elif faux.hbonds.definition == 1 : 
            faux.hbonds.sk_param=sk_param
            cut_cell=-0.343*numpy.log(sk_param/7.10)+0.10
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
                return

            if infile:

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

                        faux.hbonds.get_hbonds_skinner_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
                        gg+=1

                else:
                    hbout=[]
                    gg=0
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        faux.hbonds.get_hbonds_skinner( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)
                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
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

                        faux.hbonds.get_hbonds_skinner_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.hbonds.get_hbonds_skinner( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)
                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout


        elif faux.hbonds.definition == 2 : 
            faux.hbonds.roh2_param= roh_param**2
            print 'Not implemented yet'
            pass

        elif faux.hbonds.definition == 3 : # ROO_ANG
            faux.hbonds.roo2_param, faux.hbonds.cos_angooh_param= roo_param**2, numpy.cos(numpy.radians(angooh_param))

            if infile:

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

                        faux.hbonds.get_hbonds_roo_ang_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for aa in frame:
                        self.traj[0].reload_frame(frame=aa)
                        iframe=self.traj[0].frame[0]
                        faux.hbonds.get_hbonds_roo_ang( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)
                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
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

                        faux.hbonds.get_hbonds_roo_ang_ns_list( opt_diff_set, opt_pbc, \
                                           acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                           iframe.coors,iframe.box,iframe.orthogonal, \
                                           acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                           nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                           nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                           self.num_atoms)

                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
                        gg+=1
                else:
                    hbout=[]
                    gg=0
                    for iframe in __read_frame_opt__(self,traj,frame):
                        faux.hbonds.get_hbonds_roo_ang( opt_diff_set, opt_pbc, \
                                                        acc_don_A[0],acc_don_A[1],acc_don_A[2],acc_don_A[3],acc_don_A[4],acc_don_A[5], \
                                                        iframe.coors,iframe.box,iframe.orthogonal, \
                                                        acc_don_B[0],acc_don_B[1],acc_don_B[2],acc_don_B[3],acc_don_B[4],acc_don_B[5], \
                                                        nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, \
                                                        nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H, \
                                                        self.num_atoms)
                        hbout.append([ccopy.deepcopy(faux.glob.hbs_out),ccopy.deepcopy(faux.glob.hbs_vals_out)])
                        gg+=1

            if gg==1:
                return hbout[0]
            else:
                return hbout

        elif faux.hbonds.definition == 4 : 
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
            print 'Not implemented yet'
            pass

        elif faux.hbonds.definition == 5 : 
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
            print 'Not implemented yet'
            pass

        elif faux.hbonds.definition == 6 : 
            faux.hbonds.cos_angooh_param= numpy.cos(numpy.radians(angooh_param))
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
            print 'Not implemented yet'
            pass

        elif faux.hbonds.definition == 7 : 
            if not (allwat_A and allwat_B):
                print '# This type of hbond only works for water molecules.'
            print 'Not implemented yet'
            pass

            #if verbose:
            #    print 'Listo'
        pass


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
            dd=faux.glob.water_bisector(opt_pbc,list_atoms,iframe.coors,iframe.box,iframe.orthogonal,len(list_atoms),self.num_atoms)
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
            angls[num_frames,:,:]=faux.glob.water_angle_bisector_atom(opt_pbc,atoms,list_atoms,iframe.coors,iframe.box,iframe.orthogonal,len(list_atoms),natoms,self.num_atoms)
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

    #def lrmsd_fit(self,set_ref=None,traj_ref=None,frame_ref=None,selection=None,traj=None,frame=None,new=False):
    # 
    #    setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,set_ref,None,selection)
    # 
    #    if n_A!=n_B :
    #        print '# Error: Different number of atoms'
    #        return
    # 
    #    rot,center_ref,center_orig,rmsd,g=f.aux_funcs_general.min_rmsd(coors_reference,coors_original,len(coors_original))
    #    coors_new=f.aux_funcs_general.rot_trans(coors_original,rot,center_orig,center_ref,len(coors_original))
    # 
    #    
    #    if new:
    #        fitted_set=copy.deepcopy(self)
    #        fitted_set.frame[0].coors=coors_new
    #         fitted_set.pdb_header="Mensaje del fitteo"
    #        fitted_set.rmsd=rmsd
    # 
    #        return fitted_set
    #    else:
    #        self.frame[0].coors=copy.deepcopy(coors_new)
    #        self.rmsd=rmsd
    # 
    #    print '# RMSD fitting:',rmsd
    #        # Use coors_new
    # 
    #    return

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
        return [syst.traj[0].frame[ii] for ii in [frame]]
    elif type(frame) in [list,tuple]:
        return [syst.traj[0].frame[ii] for ii in frame]

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

    diff_syst=1
    diff_set=1

    if systB==None:

        diff_syst=0
        nsys_a=systA.num_atoms
        nsys_b=systA.num_atoms

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

        if setB in [None]:
            setB=setA
            nlist_b=nlist_a
            diff_set=0
        elif setB in ['ALL','All','all']:
            setB=[ii for ii in range(systA.num_atoms)]
            nlist_b=systA.num_atoms
        elif type(setB) in [numpy.int64,numpy.int32,int]:
            setB=[setB]
            nlist_b=1
        elif type(setB) in [list,tuple]:
            nlist_b=len(setB)
        else:
            setB=systA.selection(setB)
            nlist_b=len(setB)

        if (setA==setB): diff_set=0

        return setA,nlist_a,nsys_a,setB,nlist_b,nsys_b,diff_syst,diff_set

    else:
        print 'Different systems: Not implemented yet'
        pass
          
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
    hbs_type['Donor-Acceptor-Number']=4; hbs_info['Donor-Acceptor-Number']='A. D. Hammerich, V. J. Buch. J. Chem. Phys. 128, 111101 (2008)'
    hbs_type['Topological']=5; hbs_info['Topological']='R. H. Henchman and S. J. Irudayam. J. Phys. Chem. B. 114, 16792-16810 (2010)'
    hbs_type['Donor-Number-Ang(o,o,h)']=6; hbs_info['Donor-Number-Ang(o,o,h)']='J. D. Smith, C. D. Cappa, et al. Proc. Natl. Acad. Sci. U.S.A. 102, 14171 (2005).'
    hbs_type['Nearest-Neighbour']=7; hbs_info['Nearest-Neighbour']='This is not a hydrogen bond definition but just a topological characterization.'


    if verbose:
        if option not in hbs_type.keys():
            for ii in hbs_type.keys():
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

def select_covalent_chain(system=None,select=None,chain=None):

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
        'backbone':  '(atom.name N CA C O)',
        'sidechain': '(atom.resid.type Protein and not atom.name N CA C O H1 H2)',
        'protein':   '(atom.resid.type Protein)',
        'water':     '(atom.resid.type Water)',
        }

    for ii,jj in dict_selects.iteritems():
        icondition=icondition.replace(ii,jj)

    ### First block.
    icondition=icondition.replace(',',' ')
    icondition=icondition.replace('(',' ( ')
    icondition=icondition.replace(')',' ) ')
    icondition=icondition.replace('[',' ( ')
    icondition=icondition.replace(']',' ) ')
    icondition=icondition.replace('AND','and')
    icondition=icondition.replace('And','and')
    icondition=icondition.replace('OR','or')
    icondition=icondition.replace('Or','or')
    icondition=icondition.replace('NOT','not')
    icondition=icondition.replace('Not','not')
    icondition=icondition.replace('Within','within')
    icondition=icondition.replace('WITHIN','within')
    icondition=icondition.replace('Of','of')
    icondition=icondition.replace('OF','of')
    icondition=icondition.replace(' in ',' ')
    icondition=icondition.replace(' In ',' ')
    icondition=icondition.replace(' IN ',' ')
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
            if '(' not in ocondition[(ii+3):]:
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
                    if faux.glob.within(dists_sels[gg][jj,:],cutoff,len(sel2)):
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










