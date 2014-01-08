 # tendria que pensar si incluyo los centros de los aromaticos para hacer hbonds
import numpy
import copy
from libmss import glob as mss_funcs

class atom():

    def __init__(self,index=None,donor=False,acceptor=False,nonpolar=False):

        self.index=index
        self.name=None
        self.donor=donor
        self.acceptor=acceptor
        self.nonpolar=nonpolar

        self.bond=[]
        self.bond_node=[]
        self.bond_value=[]
        self.num_bonds=0
        self.hbond=[]
        self.hbond_node=[]
        self.hbond_value=[]
        self.num_hbonds=0

    def info(self):

        print '# ', self.name
        if self.acceptor:
            print '# acceptor'
        if self.donor:
            print '# donor'
        if self.nonpolar:
            print '# nonpolar'
        print '# hbonds:', self.num_hbonds
        print '# bonds:', self.num_bonds

class shell():

    def __init__(self):

        self.mss=[]
        self.mss_ind_atoms=[]
        self.mss_ind_nodes=[]
        self.new_symm=[]

class node():

    def __init__(self,name=None,type=None):

        self.index=None
        self.name=name
        self.label=None
        self.type=type

        self.acceptors=[]
        self.donors=[]
        self.nonpolars=[]
        self.atoms=[]
        self.atom={}
        self.num_acceptors=0
        self.num_donors=0
        self.num_nonpolars=0
        self.num_atoms=0

        self.shell1st=shell()
        self.shell2nd=shell()

    def add_atom(self,index=None,donor=False,acceptor=False,nonpolar=False):

        self.atom[index]=atom(index=index,donor=donor,acceptor=acceptor,nonpolar=nonpolar)

    def merge_node(self,node=None,update_ivars=True):
     
        for ii,jj in node.atom.iteritems():
            self.atom[ii]=jj

    def info(self,atoms=False):

        if atoms:
            for ii in self.atoms:
                self.atom[ii].info()
        else:
            print '# ',self.name,self.type
            print '# atoms:', self.num_atoms
            print '# acceptors:', self.num_acceptors
            print '# donors:', self.num_donors
            print '# nonpolars:', self.num_nonpolars


class mss():

    def __init__(self,msystem=None,sets='chains',symm_ats=None,symm_nodes=None,symm_sets_nodes=None,verbose=True):

        self.type=sets   #'chains','residues','molecules'
        self.msystem=msystem
        self.node=[]
        self.num_nodes=0
        self.num_atoms=0
        self.num_acceptors=0
        self.num_donors=0
        self.num_nonpolars=0
        self.atoms=[]
        self.acceptors=[]
        self.donors=[]
        self.nonpolars=[]

        self.proteins=[]
        self.lipids=[]
        self.ions=[]
        self.waters=[]

        self.hbtype=None
        self.btype=None

        self.build_nodes()


    def info(self):

        print '#',self.num_nodes,'nodes:'
        print '#',len(self.proteins),'in proteins'
        print '#',len(self.lipids),'in lipids'
        print '#',len(self.ions),'in ions'
        print '#',len(self.waters),'in waters'
        print '#'
        print '#',self.num_atoms,'atoms:'
        print '#',self.num_acceptors,'acceptors'
        print '#',self.num_donors,'donors'
        print '#',self.num_nonpolars,'nonpolar'

    def __build_nodes_chains__(self):

        aux_dict={}

        for ii in self.msystem.donors:
            if not aux_dict.has_key(ii[0]):
                jj=self.msystem.atom[ii[0]].resid
                tmp_node=node(jj.name+'-'+str(jj.pdb_index),jj.type)
                aux_dict[ii[0]]=tmp_node
            aux_dict[ii[0]].add_atom(ii[1],donor=True)

        for ii in self.msystem.acceptors:
            if not aux_dict.has_key(ii):
                jj=self.msystem.atom[ii].resid
                tmp_node=node(jj.name+'-'+str(jj.pdb_index),jj.type)
                aux_dict[ii]=tmp_node
            aux_dict[ii].add_atom(ii,acceptor=True)

        sel_ions=self.msystem.selection('ion')
        for ii in sel_ions:
            jj=self.msystem.atom[ii].resid
            tmp_node=node(jj.name+'-'+str(jj.pdb_index),jj.type)
            aux_dict[ii]=tmp_node
            aux_dict[ii].add_atom(ii,nonpolar=True)

        return aux_dict


    def build_nodes(self):

        if self.type=='chains+ions':

            aux_dict=self.__build_nodes_chains__()
            aux_keys=aux_dict.keys()
            aux_keys.sort()
            for ii in aux_keys:
                self.node.append(aux_dict[ii])

            del(aux_dict,aux_keys)

        if self.type=='chains+XOn+ions':
            
            aux_dict=self.__build_nodes_chains__()

            # symm in ASP,GLU and Terminals
            con_ASP =self.msystem.selection_covalent_chains(['OD1','CG','OD2'],'protein')
            con_GLU =self.msystem.selection_covalent_chains(['OE1','CD','OE2'],'protein')
            con_Term=self.msystem.selection_covalent_chains(['OC1','C','OC2'],'protein')
            for ii in con_ASP:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb,False)

            for ii in con_GLU:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb,False)

            for ii in con_Term:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb,False)

            # head of lipid AOT

            con_head=self.msystem.selection_covalent_chains(['OS1','S','OS2'],'lipid')
            for ii in con_head:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb,False)

            con_head=self.msystem.selection_covalent_chains(['OS1','S','OS3'],'lipid')
            for ii in con_head:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb,False)

            aux_keys=aux_dict.keys()
            aux_keys.sort()
            for ii in aux_keys:
                self.node.append(aux_dict[ii])

            del(aux_dict,aux_keys,con_ASP,con_GLU,con_Term)


