 # tendria que pensar si incluyo los centros de los aromaticos para hacer hbonds
import numpy
import copy
from libmss_sets import glob as mss_funcs_sets
from libmss_nosets import glob as mss_funcs_nosets


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

    def add_bond(self,atom=None,node=None,value=None):

        self.bond.append(atom)
        self.bond_node.append(node)
        self.bond_value.append(value)
        self.num_bonds+=1

    def add_hbond(self,atom=None,node=None,value=None):

        self.hbond.append(atom)
        self.hbond_node.append(node)
        self.hbond_value.append(value)
        self.num_hbonds+=1

    def sort_bonds(self,rever=False):

        if self.num_bonds>1:
            tups = zip(self.bond_value,self.bond,self.bond_node)
            tups.sort(reverse=rever)
            [self.bond_value,self.bond,self.bond_node]=zip(*tups)
            self.bond_value = list(self.bond_value)
            self.bond       = list(self.bond)
            self.bond_node  = list(self.bond_node)

    def sort_hbonds(self,rever=False):

        if self.num_hbonds>1:
            tups = zip(self.hbond_value,self.hbond,self.hbond_node)
            tups.sort(reverse=rever)
            [self.hbond_value,self.hbond,self.hbond_node]=zip(*tups)
            self.hbond_value = list(self.hbond_value)
            self.hbond       = list(self.hbond)
            self.hbond_node  = list(self.hbond_node)

    def reset(self):

        self.bond=[]
        self.bond_value=[]
        self.bond_node=[]
        self.num_bonds=0
        self.hbond=[]
        self.hbond_value=[]
        self.hbond_node=[]
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

        self.mss           = numpy.array([],dtype='int32',order='Fortran') 
        self.mss_ind_atoms = numpy.array([],dtype='int32',order='Fortran') 
        self.mss_ind_nodes = numpy.array([],dtype='int32',order='Fortran') 
        self.new_symm      = numpy.array([],dtype='int32',order='Fortran') 
        self.mss_str       = []
        self.mss_str2      = []


class node():

    def __init__(self,name=None,type=None):

        self.index=None
        self.name=name
        self.label=None
        self.type=type
        self.category=[0,0,0] # code_node,code_set,code_inside_set

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

        self.symm_ats=[]

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

        self.__sets__=False
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
        self.symm_ats=symm_ats
        self.symm_nodes=symm_nodes
        self.symm_sets_nodes=symm_sets_nodes
        self.num_categories=0

        self.build_nodes()

        self.mss_translator={}
        self.atom2node={}
        self.trad2f_node={}
        self.trad2f_atom={}
        self.trad2py_node=[]
        self.trad2py_atom=[]
        ii_n=0
        ii_at=0

        self.num_nodes=len(self.node)
        self.x_node_run_ats=numpy.zeros((self.num_nodes+1),dtype=int,order='Fortran')
        self.x_atom2node=[]
        
        for ii in range(self.num_nodes):
            node=self.node[ii]
            ii_n+=1 ; self.trad2f_node[ii]=ii_n ; self.trad2py_node.append(ii)
            if self.node[ii].type=='Protein':
                self.proteins.append(ii)
            elif self.node[ii].type=='Water':
                self.waters.append(ii)
            elif self.node[ii].type=='Ion':
                self.ions.append(ii)
            elif self.node[ii].type=='Lipid':
                self.lipids.append(ii)
            for jj,atom in sorted(node.atom.iteritems()):
                ii_at+=1 ; self.trad2f_atom[jj]=ii_at ; self.trad2py_atom.append(jj)
                self.num_atoms+=1
                node.atoms.append(jj)
                atms=self.msystem.atom[jj]
                atom.name=atms.name+'-'+str(atms.pdb_index)+'/'+atms.resid.name+'-'+str(atms.resid.pdb_index)
                self.atom2node[jj]=ii
                self.x_atom2node.append(ii_n)
                if atom.acceptor:
                    node.acceptors.append(jj)
                if atom.donor:
                    node.donors.append(jj)
                if atom.nonpolar:
                    node.nonpolars.append(jj)
            node.index=ii
            node.acceptors.sort()
            node.donors.sort()
            node.nonpolars.sort()
            node.num_acceptors=len(node.acceptors)
            node.num_donors=len(node.donors)
            node.num_nonpolars=len(node.nonpolars)
            node.num_atoms=len(node.atoms)
            self.x_node_run_ats[ii+1]=node.num_atoms+self.x_node_run_ats[ii]

        ## symmetric atoms

        self.x_symm_ats=[]
        self.x_symm_ats_crits=numpy.zeros((self.num_nodes),dtype=int,order='Fortran')
        self.x_symm_ats_start=numpy.zeros((self.num_nodes),dtype=int,order='Fortran')

        if self.symm_ats:
            symm_ats_list=[]
            if type(self.symm_ats) in [list,tuple]:
                for sel in self.symm_ats:
                    symm_ats_list.append(self.msystem.selection(sel))
            else:
                symm_ats_list.append(self.msystem.selection(sel))
            gg=1
            for node in self.node:
                self.x_symm_ats_start[node.index]=gg
                hh=0
                for criterium in symm_ats_list:
                    aa=numpy.in1d(node.atoms,criterium)
                    bb=aa.sum()
                    if bb:
                        node.symm_ats.append(aa)
                        self.x_symm_ats.append(bb)
                        self.x_symm_ats.extend(aa.nonzero()[0]+1) #aqui
                        hh+=1
                        gg+=bb+1
                self.x_symm_ats_crits[node.index]=hh
            del(symm_ats_list)
            self.x_symm_ats      =numpy.array(self.x_symm_ats,dtype=int,order='Fortran')
        else:
            for node in self.node:
                self.x_symm_ats=numpy.zeros((1),dtype=int,order='Fortran')

        ## symmetric nodes and categories
        
        translator_counter=0

        self.x_symm_nods=[]
        self.x_symm_nods_num=0

        if self.symm_nodes:
            symm_nodes_list=[]
            if type(self.symm_nodes) in [list,tuple]:
                for sel in self.symm_nodes:
                    symm_nodes_list.append(self.msystem.selection(sel))
            else:
                symm_nodes_list.append(self.msystem.selection(sel))
            self.x_symm_nods_num=len(symm_nodes_list)
            aux2=range(self.x_symm_nods_num)
            aux=[[] for ii in aux2]
            pref=['' for ii in aux2]
            for node in self.node:
                for ii in aux2:
                    if True in numpy.in1d(node.atoms,symm_nodes_list[ii]):
                        aux[ii].append(node.index+1)
                        pref[ii]=node.type[0].lower()
                        break
            for ii in aux2:
                self.x_symm_nods.append(len(aux[ii]))
                self.x_symm_nods.extend(aux[ii])
                for jj in xrange(len(aux[ii])):
                    translator_counter+=1
                    self.mss_translator[translator_counter]=pref[ii]+str(jj)
            del(aux,aux2)
            self.x_symm_nods=numpy.array(self.x_symm_nods,dtype=int,order='Fortran')


        ## symmetric sets of nodes

        self.x_symm_sets_num=0
        self.x_symm_sets=[]

        if self.symm_sets_nodes:
            self.__sets__=True
            self.x_symm_sets_num=len(self.symm_sets_nodes)
            for criterium in self.symm_sets_nodes:
                if criterium == 'lipids':
                    symm_sets_nodes=self.msystem.lipid
                    pref='l'
                elif criterium == 'proteins':
                    symm_sets_nodes=self.msystem.protein
                    pref='p'
                aux=[item for sublist in symm_sets_nodes for item in sublist]
                num_sets=len(symm_sets_nodes)
                sets={ii:[] for ii in range(num_sets)}
                for node in self.node:
                    if True in numpy.in1d(node.atoms,aux):
                        for ii in range(num_sets):
                            if True in numpy.in1d(node.atoms,symm_sets_nodes[ii]):
                                sets[ii].append(node.index+1)
                                break
                self.x_symm_sets.append(num_sets)
                self.x_symm_sets.append(len(sets[0]))
                for ii in xrange(num_sets):
                    self.x_symm_sets.extend(sets[ii])
                    for jj in xrange(len(sets[ii])):
                        translator_counter+=1
                        self.mss_translator[translator_counter]=pref+str(ii)+'-'+str(jj)

        self.x_symm_sets=numpy.array(self.x_symm_sets,dtype=int,order='Fortran')


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

    def load_net(self,hbonds=None,bonds=None,hbtype='R(o,o)-Ang(o,o,h)',btype='dists'):

        for node in self.node:
            for ii,atom in node.atom.iteritems():
                atom.reset()

        self.hbtype=hbtype
        self.btype=btype

        if self.hbtype in ['R(o,o)-Ang(o,o,h)','R(o,h)']:
            rever_hb=False
        elif self.hbtype in ['Skinner']:
            rever_hb=True
 
        if self.btype in ['dists']:
            rever_b=False
        else:
            rever_b=True
        
        if hbonds:
 
            for hb_ind,hb_val in zip(hbonds[0],hbonds[1]):
                atdon=hb_ind[1]
                atacc=hb_ind[2]
                ndon=self.atom2node[atdon]
                nacc=self.atom2node[atacc]
                self.node[ndon].atom[atdon].add_hbond(atacc,nacc,hb_val)
                self.node[nacc].atom[atacc].add_hbond(atdon,ndon,hb_val)
 
            for node in self.node:
                for atom in node.atom.values():
                    atom.sort_hbonds(rever_hb)
 
        if bonds:
 
            for bond_ind,bond_val in zip(bonds[0],bonds[1]):
                ata=bond_ind[0]
                atb=bond_ind[1]
                na=self.atom2node[ata]
                nb=self.atom2node[atb]
                self.node[na].atom[ata].add_bond(atb,nb,bond_val)
                self.node[nb].atom[atb].add_bond(ata,na,bond_val)
 
            for node in self.node:
                for atom in node.atom.values():
                    atom.sort_bonds(rever_b)


        # Hago red:

        T_num_hbs_atom=numpy.zeros((self.num_atoms),dtype=int,order='Fortran')
        T_hbs=[]
        T_num_bs_atom=numpy.zeros((self.num_atoms),dtype=int,order='Fortran')
        T_bs=[]

        jj=0
        for node in self.node:
            for ii,atom in sorted(node.atom.iteritems()):
                T_hbs.extend(atom.hbond)
                T_bs.extend(atom.bond)
                T_num_hbs_atom[jj]=atom.num_hbonds
                T_num_bs_atom[jj]=atom.num_bonds
                jj+=1

        T_hbs =numpy.array(T_hbs,dtype=int,order='Fortran')
        T_bs  =numpy.array(T_bs,dtype=int,order='Fortran')
        T_num_hbs= T_hbs.shape[0]
        T_num_bs = T_bs.shape[0]
        for ii in xrange(T_num_hbs):
            T_hbs[ii]=self.trad2f_atom[T_hbs[ii]]
        for ii in xrange(T_num_bs):
            T_bs[ii]=self.trad2f_atom[T_bs[ii]]


        if self.__sets__ :
         
            mss_funcs_sets.load_net(T_hbs,T_bs,T_num_hbs_atom,T_num_bs_atom,
                                    T_num_hbs,T_num_bs,self.num_atoms)
         
        else:
         
            mss_funcs_nosets.load_net(T_hbs,T_bs,T_num_hbs_atom,T_num_bs_atom,
                                      T_num_hbs,T_num_bs,self.num_atoms)


        del(T_hbs,T_bs,T_num_hbs_atom,T_num_bs_atom,T_num_hbs,T_num_bs)

    def load_topol(self):

        if self.__sets__ :

            mss_funcs_sets.load_topol(self.x_atom2node,self.trad2py_atom,self.trad2py_node,
                                      self.x_symm_ats_start,self.x_symm_ats_crits,self.x_symm_ats,
                                      self.x_symm_nods,self.x_symm_nods_num,
                                      self.x_symm_sets,self.x_symm_sets_num,
                                      self.num_atoms,self.num_nodes,
                                      self.x_symm_ats.shape[0],self.x_symm_nods.shape[0],self.x_symm_sets.shape[0])
        else:

            mss_funcs_nosets.load_topol(self.x_atom2node,self.trad2py_atom,self.trad2py_node,
                                        self.x_symm_ats_start,self.x_symm_ats_crits,self.x_symm_ats,
                                        self.x_symm_nods,self.x_symm_nods_num,
                                        self.num_atoms,self.num_nodes,
                                        self.x_symm_ats.shape[0],self.x_symm_nods.shape[0])


    def build_shell1st(self,node='ALL',tostr=False):
     
        # x_node_run_ats,x_atom2node,trad2py_node,trad2py_atom
        # T_hbs_start,T_hbs_ind,T_bs_start,T_bs_ind,T_num_hbs,T_num_bs
        if self.__sets__ :
            if node=='ALL':
                for ii in range(self.num_nodes):
                    jj=self.trad2f_node[ii]
                    mss_funcs_sets.build_shell1st(jj)
                    self.node[ii].shell1st.mss           = numpy.copy(mss_funcs_sets.mss)
                    self.node[ii].shell1st.mss_ind_atoms = numpy.copy(mss_funcs_sets.mss_ind_ats)
                    self.node[ii].shell1st.mss_ind_nodes = numpy.copy(mss_funcs_sets.mss_ind_nods)
                    self.node[ii].shell1st.mss_symm      = numpy.copy(mss_funcs_sets.mss_symm)
                if tostr:
                    self.mss2mss_str(shell1st=True)
            elif type(node)==int:
                jj=self.trad2f_node[node]
                mss_funcs_sets.build_shell1st(jj)
                self.node[node].shell1st.mss           = numpy.copy(mss_funcs_sets.mss)
                self.node[node].shell1st.mss_ind_atoms = numpy.copy(mss_funcs_sets.mss_ind_ats)
                self.node[node].shell1st.mss_ind_nodes = numpy.copy(mss_funcs_sets.mss_ind_nods)
                self.node[node].shell1st.mss_symm      = numpy.copy(mss_funcs_sets.mss_symm)
                if tostr:
                    self.mss2mss_str(node=node,shell1st=True)
            pass
        else:
            if node=='ALL':
                for ii in range(self.num_nodes):
                    jj=self.trad2f_node[ii]
                    mss_funcs_nosets.build_shell1st(jj)
                    self.node[ii].shell1st.mss           = numpy.copy(mss_funcs_nosets.mss)
                    self.node[ii].shell1st.mss_ind_atoms = numpy.copy(mss_funcs_nosets.mss_ind_ats)
                    self.node[ii].shell1st.mss_ind_nodes = numpy.copy(mss_funcs_nosets.mss_ind_nods)
                    self.node[ii].shell1st.mss_symm      = numpy.copy(mss_funcs_nosets.mss_symm)
                if tostr:
                    self.mss2mss_str(shell1st=True)
            elif type(node)==int:
                jj=self.trad2f_node[node]
                mss_funcs_nosets.build_shell1st(jj)
                self.node[node].shell1st.mss           = numpy.copy(mss_funcs_nosets.mss)
                self.node[node].shell1st.mss_ind_atoms = numpy.copy(mss_funcs_nosets.mss_ind_ats)
                self.node[node].shell1st.mss_ind_nodes = numpy.copy(mss_funcs_nosets.mss_ind_nods)
                self.node[node].shell1st.mss_symm      = numpy.copy(mss_funcs_nosets.mss_symm)
                if tostr:
                    self.mss2mss_str(node=node,shell1st=True)
            pass
 


    def build_shell2nd(self,node='ALL',tostr=False):

        if self.__sets__ :    
            if node=='ALL':
                for ii in range(self.num_nodes):
                    jj=self.trad2f_node[ii]
                    mss_funcs_sets.build_shell2nd(jj)
                    self.node[ii].shell2nd.mss           = numpy.copy(mss_funcs_sets.mss)
                    self.node[ii].shell2nd.mss_ind_atoms = numpy.copy(mss_funcs_sets.mss_ind_ats)
                    self.node[ii].shell2nd.mss_ind_nodes = numpy.copy(mss_funcs_sets.mss_ind_nods)
                    self.node[ii].shell2nd.mss_symm      = numpy.copy(mss_funcs_sets.mss_symm)
                if tostr:
                    self.mss2mss_str(shell2nd=True)
            elif type(node)==int:
                jj=self.trad2f_node[node]
                mss_funcs_sets.build_shell2nd(jj)
                self.node[node].shell2nd.mss           = numpy.copy(mss_funcs_sets.mss)
                self.node[node].shell2nd.mss_ind_atoms = numpy.copy(mss_funcs_sets.mss_ind_ats)
                self.node[node].shell2nd.mss_ind_nodes = numpy.copy(mss_funcs_sets.mss_ind_nods)
                self.node[node].shell2nd.mss_symm      = numpy.copy(mss_funcs_sets.mss_symm)
                if tostr:
                    self.mss2mss_str(node=node,shell2nd=True)
        else:
            if node=='ALL':
                for ii in range(self.num_nodes):
                    jj=self.trad2f_node[ii]
                    mss_funcs_nosets.build_shell2nd(jj)
                    self.node[ii].shell2nd.mss           = numpy.copy(mss_funcs_nosets.mss)
                    self.node[ii].shell2nd.mss_ind_atoms = numpy.copy(mss_funcs_nosets.mss_ind_ats)
                    self.node[ii].shell2nd.mss_ind_nodes = numpy.copy(mss_funcs_nosets.mss_ind_nods)
                    self.node[ii].shell2nd.mss_symm      = numpy.copy(mss_funcs_nosets.mss_symm)
                if tostr:
                    self.mss2mss_str(shell2nd=True)
            elif type(node)==int:
                jj=self.trad2f_node[node]
                mss_funcs_nosets.build_shell2nd(jj)
                self.node[node].shell2nd.mss           = numpy.copy(mss_funcs_nosets.mss)
                self.node[node].shell2nd.mss_ind_atoms = numpy.copy(mss_funcs_nosets.mss_ind_ats)
                self.node[node].shell2nd.mss_ind_nodes = numpy.copy(mss_funcs_nosets.mss_ind_nods)
                self.node[node].shell2nd.mss_symm      = numpy.copy(mss_funcs_nosets.mss_symm)
                if tostr:
                    self.mss2mss_str(node=node,shell2nd=True)
        pass
 
    def filtrosymm(self,shell1st=False,shell2nd=False,node='ALL'):
 
        if node=='ALL':
            aa=self.node
        elif type(node)==int:
            aa=[self.node[node]]
 
        if shell2nd==True:
 
            for node in aa:
 
                gg=0
                nn1=node.shell2nd.mss_ind_nodes[gg]
                for ii in xrange(nn1):
                    gg+=1
                nn11=0
                for ii in xrange(nn1*2):
                    gg+=1
                    nn11+=node.shell2nd.mss_ind_nodes[gg]
                for ii in xrange(nn11):
                    gg+=1
                interruptor=False
                for ii in xrange(nn11):
                    gg+=1
                    bb=gg
                    nn2=node.shell2nd.mss_ind_nodes[gg]
                    veo=[]
                    for jj in xrange(nn2):
                        gg+=1
                        veo.append(node.shell2nd.mss_symm[gg])
                    interruptor=False
                    if 1 in veo:
                        ayuda=[]
                        for ii in xrange(nn2*2):
                            gg+=1
                            ayuda.append(node.shell2nd.mss_ind_nodes[gg])
                        ll=-1
                        for dd in xrange(nn2):
                            for kk in xrange(2):
                                ll+=1
                                if veo[dd]==1:
                                    tt=0
                                    for pp in xrange(ayuda[ll]):
                                        gg+=1
                                        if (node.shell2nd.mss_str.count(node.shell2nd.mss_str[gg])==1):
                                            tt+=1
                                    if ayuda[ll]>0:
                                        if ayuda[ll]==tt:
                                            node.shell2nd.mss_symm[dd+1+bb]=0
                                            interruptor=True
                                else:
                                    for jj in xrange(ayuda[ll]):
                                        gg+=1
                    else:
                        nn22=0
                        for jj in xrange(nn2*2):
                            gg+=1
                            nn22+=node.shell2nd.mss_ind_nodes[gg]
                        for jj in xrange(nn22):
                            gg+=1
                #if interruptor==True and sum(node.shell2nd.mss_symm)==0:
                #    print node.index
 
            for node in aa:
 
                gg=0
                nn1=node.shell2nd.mss_ind_nodes[gg]
                for ii in xrange(nn1):
                    gg+=1
                nn11=0
                for ii in xrange(nn1*2):
                    gg+=1
                    nn11+=node.shell2nd.mss_ind_nodes[gg]
                for ii in xrange(nn11):
                    gg+=1
                interruptor=False
                for ii in xrange(nn11):
                    gg+=1
                    bb=gg
                    nn2=node.shell2nd.mss_ind_nodes[gg]
                    for jj in xrange(nn2):
                        gg+=1
                    nn22=0
                    for jj in xrange(nn2*2):
                        gg+=1
                        nn22+=node.shell2nd.mss_ind_nodes[gg]
                    for jj in xrange(nn22):
                        gg+=1
                        if node.shell2nd.mss_symm[gg]==1:
                            if node.shell2nd.mss_str[gg].startswith('w'):
                                if (node.shell2nd.mss_str.count(node.shell2nd.mss_str[gg])==1):
                                    node.shell2nd.mss_symm[gg]=0
                                    interruptor=True
                            if node.shell2nd.mss_str[gg].startswith('i'):
                                if (node.shell2nd.mss_str.count(node.shell2nd.mss_str[gg])==1):
                                    node.shell2nd.mss_symm[gg]=0
                                    interruptor=True
                #if interruptor==True and (sum(node.shell2nd.mss_symm)==0):
                #    print node.index
 
 
            for node in aa:
 
                gg=0
                nn1=node.shell2nd.mss_ind_nodes[gg]
                for ii in xrange(nn1):
                    gg+=1
                nn11=0
                for ii in xrange(nn1*2):
                    gg+=1
                    nn11+=node.shell2nd.mss_ind_nodes[gg]
                veo=[]
                bb=gg
                for ii in xrange(nn11):
                    gg+=1
                    if (node.shell2nd.mss_symm[gg]==1):
                        veo.append(True)
                    else:
                        veo.append(False)
                interruptor2=False
                interruptor=True
                if True in veo:
                    #print 'Aqui puede haber', node.index
                    interruptor2=True
                for ii in xrange(nn11):
                    if veo[ii]==False:
                        gg+=1
                        nn2=node.shell2nd.mss_ind_nodes[gg]
                        for jj in xrange(nn2):
                            gg+=1
                        nn22=0
                        for jj in xrange(nn2*2):
                            gg+=1
                            nn22+=node.shell2nd.mss_ind_nodes[gg]
                        for jj in xrange(nn22):
                            gg+=1
                    else:
                        gg+=1
                        nn2=node.shell2nd.mss_ind_nodes[gg]
                        for jj in xrange(nn2):
                            gg+=1
                            if (node.shell2nd.mss_symm[gg]==1):
                                interruptor=False
                        nn22=0
                        for jj in xrange(nn2*2):
                            gg+=1
                            nn22+=node.shell2nd.mss_ind_nodes[gg]
                        for jj in xrange(nn22):
                            gg+=1
                            if (node.shell2nd.mss_symm[gg]==1):
                                interruptor=False
                if interruptor==True and interruptor2==True:
                    #print '......pues hay'
                    gg=bb
                    for ii in xrange(nn11):
                        gg+=1
                        if veo[ii]==True:
                            node.shell2nd.mss_symm[gg]=0
                    #if sum(node.shell2nd.mss_symm[:])==0:
                    #    print 'le di'
 
 
            for node in aa:
 
                gg=0
                nn1=node.shell2nd.mss_ind_nodes[gg]
                veo=[]
                bb=gg
                interruptor2=False
                interruptor=True
                for ii in xrange(nn1):
                    gg+=1
                    if (node.shell2nd.mss_symm[gg]==1):
                        veo.append(True)
                    else:
                        veo.append(False)
                if True in veo:
                    interruptor2=True
                nn11=0
                lista=[]
                for ii in xrange(nn1*2):
                    gg+=1
                    nn11+=node.shell2nd.mss_ind_nodes[gg]
                    lista.append(node.shell2nd.mss_ind_nodes[gg])
                ggg=0
                for ii in xrange(nn1):
                    for jj in range(lista[ggg]):
                        gg+=1
                        if veo[ii]:
                            if (node.shell2nd.mss_symm[gg]==1):
                                interruptor=False
                    ggg+=1
                    for jj in range(lista[ggg]):
                        gg+=1
                        if veo[ii]:
                            if (node.shell2nd.mss_symm[gg]==1):
                                interruptor=False
                    ggg+=1
                ggg=0
                for ii in xrange(nn1):
                    if veo[ii]==False:
                        for jj in range(lista[ggg]):
                            gg+=1
                            nn2=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn2):
                                gg+=1
                            nn22=0
                            for jjj in xrange(nn2*2):
                                gg+=1
                                nn22+=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn22):
                                gg+=1
                        ggg+=1
                        for jj in range(lista[ggg]):
                            gg+=1
                            nn2=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn2):
                                gg+=1
                            nn22=0
                            for jjj in xrange(nn2*2):
                                gg+=1
                                nn22+=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn22):
                                gg+=1
                        ggg+=1
                    else:
                        for jj in range(lista[ggg]):
                            gg+=1
                            nn2=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn2):
                                gg+=1
                                if (node.shell2nd.mss_symm[gg]==1):
                                    interruptor=False
                            nn22=0
                            for jjj in xrange(nn2*2):
                                gg+=1
                                nn22+=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn22):
                                gg+=1
                                if (node.shell2nd.mss_symm[gg]==1):
                                    interruptor=False
                        ggg+=1
                        for jj in range(lista[ggg]):
                            gg+=1
                            nn2=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn2):
                                gg+=1
                                if (node.shell2nd.mss_symm[gg]==1):
                                    interruptor=False
                            nn22=0
                            for jjj in xrange(nn2*2):
                                gg+=1
                                nn22+=node.shell2nd.mss_ind_nodes[gg]
                            for jjj in xrange(nn22):
                                gg+=1
                                if (node.shell2nd.mss_symm[gg]==1):
                                    interruptor=False
                        ggg+=1
                if interruptor==True and interruptor2==True:
                    #print '......pues hay'
                    gg=bb
                    for ii in xrange(nn1):
                        gg+=1
                        if veo[ii]==True:
                            node.shell2nd.mss_symm[gg]=0
                    #if sum(node.shell2nd.mss_symm[:])==0:
                    #    print 'le di'
  
        if shell1st==True:
 
            for node in aa:
                gg=0
                nn1=node.shell1st.mss_ind_nodes[gg]
                veo=[]
                for ii in xrange(nn1):
                    gg+=1
                    veo.append(node.shell1st.mss_symm[gg])
                interruptor=False
                if 1 in veo:
                    ayuda=[]
                    for ii in xrange(nn1*2):
                        gg+=1
                        ayuda.append(node.shell1st.mss_ind_nodes[gg])
                    ll=-1
                    for ii in xrange(nn1):
                        for kk in xrange(2):
                            ll+=1
                            if veo[ii]==1:
                                tt=0
                                for pp in xrange(ayuda[ll]):
                                    gg+=1
                                    if (node.shell1st.mss_str.count(node.shell1st.mss_str[gg])==1):
                                        tt+=1
                                if ayuda[ll]>0:
                                    if ayuda[ll]==tt:
                                        node.shell1st.mss_symm[ii+1]=0
                                        interruptor=True
                            else:
                                for jj in xrange(ayuda[ll]):
                                    gg+=1
                    #if interruptor==True and sum(node.shell1st.mss_symm)==0:
                    #    print node.index
 
            for node in aa:
                interruptor=False
                gg=0
                nn1=node.shell1st.mss_ind_nodes[gg]
                for ii in xrange(nn1):
                    gg+=1
                nn2=0
                for ii in xrange(nn1*2):
                    gg+=1
                    nn2+=node.shell1st.mss_ind_nodes[gg]
                for ii in xrange(nn2):
                        gg+=1
                        if node.shell1st.mss_symm[gg]==1:
                            if node.shell1st.mss_str[gg].startswith('w'):
                                if (node.shell1st.mss_str.count(node.shell1st.mss_str[gg])==1):
                                    node.shell1st.mss_symm[gg]=0
                                    interruptor=True
                            if node.shell1st.mss_str[gg].startswith('i'):
                                if (node.shell1st.mss_str.count(node.shell1st.mss_str[gg])==1):
                                    node.shell1st.mss_symm[gg]=0
                                    interruptor=True
                #if interruptor==True and (sum(node.shell1st.mss_symm)==0):
                #    print node.index
 
    def mss2mss_str(self,node='ALL',shell1st=False,shell2nd=False):
 
        if node=='ALL':
            aa=self.node
        elif type(node)==int:
            aa=[self.node[node]]
 
        if shell1st:
 
            for node in aa:
                node.shell1st.mss_str=node.shell1st.mss.tolist()
                n_ats=node.shell1st.mss[0]
                aux_str=self.mss_translator[node.shell1st.mss[1]]
                ii=1
                jj=ii+n_ats
                for iii in xrange(ii,jj):
                    node.shell1st.mss_str[iii]=aux_str
                ii=jj
                jj=jj+n_ats*2
                n_bs=node.shell1st.mss[ii:jj].sum()
                ii=jj
                jj=jj+n_bs
                for iii in xrange(ii,jj):
                    node.shell1st.mss_str[iii]=self.mss_translator[node.shell1st.mss[iii]]
                node.shell1st.mss_str2=mss1trans(node.shell1st.mss_str)
 
        if shell2nd:
 
            for node in aa:
                node.shell2nd.mss_str=node.shell2nd.mss.tolist()
                n_ats=node.shell2nd.mss[0]
                aux_str=self.mss_translator[node.shell2nd.mss[1]]
                ii=1
                jj=ii+n_ats
                for iii in xrange(ii,jj):
                    node.shell2nd.mss_str[iii]=aux_str
                ii=jj
                jj=jj+n_ats*2
                n_bs=node.shell2nd.mss[ii:jj].sum()
                ii=jj
                jj=jj+n_bs
                for iii in xrange(ii,jj):
                    node.shell2nd.mss_str[iii]=self.mss_translator[node.shell2nd.mss[iii]]
                for jjj in xrange(n_bs):
                    n_ats2=node.shell2nd.mss[jj]
                    aux_str=self.mss_translator[node.shell2nd.mss[jj+1]]
                    jj+=1
                    ii=jj
                    jj=jj+n_ats2
                    for iii in xrange(ii,jj):
                        node.shell2nd.mss_str[iii]=aux_str
                    ii=jj
                    jj=jj+n_ats2*2
                    n_bs2=node.shell2nd.mss[ii:jj].sum()
                    ii=jj
                    jj=jj+n_bs2
                    for iii in xrange(ii,jj):
                        node.shell2nd.mss_str[iii]=self.mss_translator[node.shell2nd.mss[iii]]
                node.shell2nd.mss_str2=mss2trans(node.shell2nd.mss_str)
        pass
 
def mss1trans(mss):
 
    gg=0
    mss_out='<'
    ii=mss[gg]
    mss_out+=mss[gg+1]
    gg=gg+ii
    aa=0
    mss_out+='|'
    joe=[]
    for jj in xrange(ii*2):
        gg+=1
        joe.append(str(mss[gg]))
        aa+=mss[gg]
    mss_out+=str.join(',',joe)
    mss_out+='|'
    joe=[]
    for jj in xrange(aa):
        gg+=1
        joe.append(mss[gg])
    mss_out+=str.join(',',joe)
    mss_out+='>'
    return mss_out
 
def mss2trans(mss):
 
    gg=0
    mss_out='<'
    ii=mss[gg]
    mss_out+=mss[gg+1]
    gg=gg+ii
    aa=0
    mss_out+='|'
    joe=[]
    for jj in xrange(ii*2):
        gg+=1
        joe.append(str(mss[gg]))
        aa+=mss[gg]
    mss_out+=str.join(',',joe)
    mss_out+='|'
    joe=[]
    for jj in xrange(aa):
        gg+=1
        joe.append(mss[gg])
    mss_out+=str.join(',',joe)
    for jj in xrange(aa):
        mss_out+='> <'
        gg+=1
        ii=mss[gg]
        mss_out+=mss[gg+1]
        gg=gg+ii
        bb=0
        mss_out+='|'
        joe=[]
        for jj in xrange(ii*2):
            gg+=1
            joe.append(str(mss[gg]))
            bb+=mss[gg]
        mss_out+=str.join(',',joe)
        mss_out+='|'
        joe=[]
        for kk in xrange(bb):
            gg+=1
            joe.append(mss[gg])
        mss_out+=str.join(',',joe)
    mss_out+='>'
    return mss_out
        
