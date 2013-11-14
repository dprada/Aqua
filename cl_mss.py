# tendria que pensar si incluyo los centros de los aromaticos para hacer hbonds
import numpy

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

        self.mss=[]
        self.mss_ind=[]

class node():

    def __init__(self,name=None,type=None):

        self.type=type
        self.name=name
        self.index=None

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

    def merge_node(self,node=None):

        for ii,jj in node.atom.iteritems():
            self.atom[ii]=jj


    def info(self):

        print '# ',self.name,self.type
        print '# atoms:', self.num_atoms
        print '# acceptors:', self.num_acceptors
        print '# donors:', self.num_donors
        print '# nonpolars:', self.num_nonpolars

class mss():

    def __init__(self,msystem=None,sets='chains',symm_ats=None,symm_nodes=None,verbose=True):

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
        self.atom2node={}
        self.acceptor2node={}
        self.donor2node={}
        self.nonpolar2node={}

        self.proteins=[]
        self.lipids=[]
        self.ions=[]
        self.waters=[]

        self.hbtype=None
        self.btype=None
        self.symm_ats=symm_ats
        self.symm_nodes=symm_nodes

        self.build_nodes()

        self.num_nodes=len(self.node)
        for ii in range(self.num_nodes):
            node=self.node[ii]
            if self.node[ii].type=='Protein':
                self.proteins.append(ii)
            elif self.node[ii].type=='Water':
                self.waters.append(ii)
            elif self.node[ii].type=='Ion':
                self.ions.append(ii)
            elif self.node[ii].type=='Lipid':
                self.lipids.append(ii)
            for jj,atom in node.atom.iteritems():
                node.atoms.append(jj)
                atms=self.msystem.atom[jj]
                atom.name=atms.name+'-'+str(atms.pdb_index)+'/'+atms.resid.name+'-'+str(atms.resid.pdb_index)
                self.atom2node[jj]=ii
                if atom.acceptor:
                    node.acceptors.append(jj)
                    self.acceptor2node[jj]=ii
                if atom.donor:
                    node.donors.append(jj)
                    self.donor2node[jj]=ii
                if atom.nonpolar:
                    node.nonpolars.append(jj)
                    self.nonpolar2node[jj]=ii
            node.index=ii
            node.acceptors.sort()
            node.donors.sort()
            node.nonpolars.sort()
            node.atoms.sort()
            node.num_acceptors=len(node.acceptors)
            node.num_donors=len(node.donors)
            node.num_nonpolars=len(node.nonpolars)
            node.num_atoms=len(node.atoms)
        self.atoms=self.atom2node.keys()
        self.acceptors=self.acceptor2node.keys()
        self.donors=self.donor2node.keys()
        self.nonpolars=self.nonpolar2node.keys()
        self.num_atoms=len(self.atoms)
        self.num_acceptors=len(self.acceptors)
        self.num_donors=len(self.donors)
        self.num_nonpolars=len(self.nonpolars)

        ## symmetric atoms
         
        if self.symm_ats:
            symm_ats_list=[]
            if type(self.symm_ats) in [list,tuple]:
                for sel in self.symm_ats:
                    symm_ats_list.append(self.msystem.selection(sel))
            else:
                symm_ats_list.append(self.msystem.selection(sel))
            for node in self.node:
                for criterium in symm_ats_list:
                    #aa=numpy.in1d(node.atoms,criterium)
                    aa=numpy.intersect1d(node.atoms,criterium)
                    if aa.sum():
                        node.symm_ats.append(aa)
            del(symm_ats_list)

        if verbose:
            self.info()

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

        if self.type=='chains':

            aux_dict=self.__build_nodes_chains__()
            aux_keys=aux_dict.keys()
            aux_keys.sort()
            for ii in aux_keys:
                self.node.append(aux_dict[ii])

            del(aux_dict,aux_keys)

        if self.type=='chains+XOn':
            
            aux_dict=self.__build_nodes_chains__()

            # symm in ASP,GLU and Terminals
            con_ASP =self.msystem.selection_covalent_chains(['OD1','CG','OD2'],'protein')
            con_GLU =self.msystem.selection_covalent_chains(['OE1','CD','OE2'],'protein')
            con_Term=self.msystem.selection_covalent_chains(['OC1','C','OC2'],'protein')
            for ii in con_ASP:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb)

            for ii in con_GLU:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb)

            for ii in con_Term:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb)

            # head of lipid AOT

            con_head=self.msystem.selection_covalent_chains(['OS1','S','OS2'],'lipid')
            for ii in con_head:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb)

            con_head=self.msystem.selection_covalent_chains(['OS1','S','OS3'],'lipid')
            for ii in con_head:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].merge_node(bb)

            aux_keys=aux_dict.keys()
            aux_keys.sort()
            for ii in aux_keys:
                self.node.append(aux_dict[ii])

            del(aux_dict,aux_keys,con_ASP,con_GLU,con_Term)

    def reset(self):

        for node in self.node:  
            for atom in node.atom.values():
                atom.reset()

    def build_shell1st(self,hbonds=None,bonds=None,hbtype='R(o,o)-Ang(o,o,h)',btype='dists'):

        self.hbtype=hbtype
        self.btype=btype

        self.reset()

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
                ndon=self.donor2node[atdon]
                nacc=self.acceptor2node[atacc]
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


    #def build_mss_ind_shell1st(self):
    # 
    #    for node in self.node:
    #        mss_ind=[]
    #        mss_ind.extend([node.don_num,node.acc_num,node.at_num])
    #        mss_ind.extend(node.shell1st.don_num)
    #        mss_ind.extend(node.shell1st.acc_num)
    #        mss_ind.extend(node.shell1st.bond_num)
    #        for ii in range(node.don_num):
    #            for jj in range(node.shell1st.don_num[ii]):
    #                mss_ind.append(node.shell1st.don_node[ii][jj][0])
    #        for ii in range(node.acc_num):
    #            for jj in range(node.shell1st.acc_num[ii]):
    #                mss_ind.append(node.shell1st.acc_node[ii][jj][0])
    #        for ii in range(node.at_num):
    #            for jj in range(node.shell1st.bond_num[ii]):
    #                mss_ind.append(node.shell1st.bond_node[ii][jj][0])
    #        node.shell1st.mss_ind=mss_ind
    #        
    #        node.shell1st.mss_ind_symm=mss_ind
    #        for eqats in node.filt_ats_symm:
    #            pass
    # 
    #    for node in self.node:
    #        mss_ind=[]
    #        mss_ind.extend([node.don_num,node.acc_num,node.at_num])
    # 
    #        
    # 
    #    
    # 
    #def symmetrize_mss_ind_shell1st(self,symm_type=None):
    # 
    #    if symm_type==1: # microstate only for water
    #        
    #        # H1 and H2 distinguishable in first shell
    #        # H1 and H2 indistinguishable in second shell
    #        # Nodes water indistinguishable
    #        # Nodes ion indistinguishable
    #        # Nodes lipid indistinguishable with different peptides ??
    # 
    #        for ii in self.list_water:
    #            ishell1st=self.node[ii].shell1st
    #            ishell1st.mss_ind_sym=ishell1st.mss_ind[3:7]
    #            ishell1st.mss_sym=ishell1st.mss_ind[3:7]
    #            numhb=sum(ishell1st.mss_ind_sym[0:3])
    #            numb=ishell1st.mss_ind_sym[3]
    #            wind=0
    #            iind=0
    #            aux_dict={}
    #            for kk in ishell1st.mss_ind[9:(9+numhb)]:
    #                if not aux_dict.has_key(kk):
    #                    if self.node[kk].type=='Water':
    #                        aux_dict[kk]='w'+str(wind)
    #                        wind+=1
    #                    else:
    #                        aux_dict[kk]=kk
    #                ishell1st.mss_ind_sym.append(kk)
    #                ishell1st.mss_sym.append(aux_dict[kk])
    #            for kk in ishell1st.mss_ind[(9+numhb):(9+numhb+numb)]:
    #                if not aux_dict.has_key(kk):
    #                    if self.node[kk].type=='Ion':
    #                        aux_dict[kk]='i'+str(iind)
    #                        iind+=1
    #                    else:
    #                        aux_dict[kk]=kk
    #                ishell1st.mss_ind_sym.append(kk)
    #                ishell1st.mss_sym.append(aux_dict[kk])
    #            print aux_dict
    #        del(aux_dict)
    # 
    #def build_mss_ind(self):
    # 
    #    self.build_mss_ind_shell1st()
    # 
    #    for node in self.node:
    #        node.mss_ind.append(node.index)
    #        node.mss_ind.extend(node.shell1st.mss_ind)
    #        ii=sum(node.shell1st.mss_ind[0:3])+3
    #        jj=sum(node.shell1st.mss_ind[3:ii])+ii
    #        for kk in node.shell1st.mss_ind[ii:jj]:
    #            node.mss_ind.extend(self.node[kk].shell1st.mss_ind)
            

