# tendria que pensar si incluyo los centros de los aromaticos para hacer hbonds
import numpy

class shell1st():

    def __init__(self):

        self.acc=[]
        self.don=[]
        self.bond=[]
        self.acc_node=[]
        self.don_node=[]
        self.bond_node=[]
        self.acc_val=[]
        self.don_val=[]
        self.bond_val=[]
        self.acc_num=[]
        self.don_num=[]
        self.bond_num=[]
        self.mss=[]
        self.mss_ind=[]
        self.mss_ind_symm=[]
        self.mss_symm=[]


class adnode():

    def __init__(self):

        self.type=None
        self.name=None
        self.index=None
        self.acceptors=[]
        self.donors=[]
        self.atoms=[]
        self.acc_num=0
        self.don_num=0
        self.at_num=0
        self.shell1st=shell1st()
        self.mss=[]
        self.mss_ind=[]
        self.symm_ats=None
        self.symm_ats_list=None
        self.symm_nodes=None
        self.symm_nodes_list=None
        self.filt_symm_ats=[]

    def info(self):

        print '# acceptors:', self.acc_num
        print '# donors:', self.don_num
        print '# atoms:', self.at_num

class mss():

    def __init__(self,msystem=None,sets='chains',symm_ats=None,symm_nodes=None,verbose=True):

        self.type=sets   #'chains','residues','molecules'
        self.msystem=msystem
        self.node=[]
        self.num_nodes=0
        self.hbtype=None
        self.btype=None
        self.symm_ats=symm_ats
        self.symm_nodes=symm_nodes

        self.build_nodes()
        self.num_nodes=len(self.node)
        
        self.list_protein=[]
        self.list_water=[]
        self.list_ion=[]
        self.list_lipid=[]
        for ii in range(self.num_nodes):
            if self.node[ii].type=='Protein':
                self.list_protein.append(ii)
            elif self.node[ii].type=='Water':
                self.list_water.append(ii)
            elif self.node[ii].type=='Ion':
                self.list_ion.append(ii)
            elif self.node[ii].type=='Lipid':
                self.list_lipid.append(ii)


        self.at2node={}
        self.acc2node={}
        self.don2node={}

        for ii in range(self.num_nodes):
            node=self.node[ii]
            node.index=ii
            node.don_num=len(node.donors)
            node.acc_num=len(node.acceptors)
            node.at_num=len(node.atoms)
            for jj in range(node.don_num):
                self.don2node[node.donors[jj]]=[ii,jj]
                node.shell1st.don.append([])
                node.shell1st.don_val.append([])
                node.shell1st.don_num.append(0)
            for jj in range(node.acc_num):
                self.acc2node[node.acceptors[jj]]=[ii,jj]
                node.shell1st.acc.append([])
                node.shell1st.acc_val.append([])
                node.shell1st.acc_num.append(0)
            for jj in range(node.at_num):
                self.at2node[node.atoms[jj]]=[ii,jj]
                node.shell1st.bond.append([])
                node.shell1st.bond_val.append([])
                node.shell1st.bond_num.append(0)

        self.at_list=self.at2node.keys()
        self.acc_list=self.acc2node.keys()
        self.don_list=self.don2node.keys()
        self.at_num=len(self.at_list)
        self.acc_num=len(self.acc_list)
        self.don_num=len(self.don_list)

        self.__dict_aux_don__={}
        self.__dict_aux_acc__={}
        self.__dict_aux_at__={}

        jj=0
        for ii in self.don_list:
            self.__dict_aux_don__[ii]=jj
            jj+=1

        jj=0
        for ii in self.acc_list:
            self.__dict_aux_acc__[ii]=jj
            jj+=1

        jj=0
        for ii in self.at_list:
            self.__dict_aux_at__[ii]=jj
            jj+=1


        # symmetric atoms

        if self.symm_ats:

            self.symm_ats_list=[]
            if type(self.symm_ats) in [list,tuple]:
                for sel in self.symm_ats:
                    self.symm_ats_list.append(self.msystem.selection(sel))
            else:
                self.symm_ats_list.append(self.msystem.selection(sel))

        for inode in self.node:
            for criterium in self.symm_ats_list:
                #aa=numpy.in1d(inode.atoms,criterium)
                aa=numpy.intersect1d(inode.atoms,criterium)
                if aa:
                    inode.filt_symm_ats.append(aa)


        if verbose:
            self.info()

    def info(self):

        print '#',self.num_nodes,'nodes:'
        print '#',len(self.list_protein),'in proteins'
        print '#',len(self.list_lipid),'in lipids'
        print '#',len(self.list_ion),'in ions'
        print '#',len(self.list_water),'in waters'
        print '#'
        print '#',self.at_num,'atoms:'
        print '#',self.acc_num,'acceptors'
        print '#',self.don_num,'donors'
        print '#',self.at_num-self.acc_num-self.don_num,'nonpolar'

    def build_nodes(self):

        if self.type=='chains':
            
            aux_dict={}
            aa=0
            for ii in self.msystem.donors:
                if not aux_dict.has_key(ii[0]):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii[0]].resid.type
                    tmp_adnode.name=self.msystem.atom[ii[0]].resid.name
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii[0]].resid.water
                    aux_dict[ii[0]]=tmp_adnode
                aux_dict[ii[0]].donors.append(ii[1])

            for ii in self.msystem.acceptors:
                if not aux_dict.has_key(ii):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii].resid.type
                    tmp_adnode.name=self.msystem.atom[ii].resid.name
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii].resid.water
                    aux_dict[ii]=tmp_adnode
                aux_dict[ii].acceptors.append(ii)

            aux_keys=aux_dict.keys()
            aux_keys.sort()

            for ii in aux_keys:
                self.node.append(aux_dict[ii])
            
            del(aux_dict,aux_keys)

        if self.type=='chains+XOn':
            
            aux_dict={}
            aa=0
            for ii in self.msystem.donors:
                if not aux_dict.has_key(ii[0]):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii[0]].resid.type
                    tmp_adnode.name=self.msystem.atom[ii[0]].resid.name
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii[0]].resid.water
                    aux_dict[ii[0]]=tmp_adnode
                aux_dict[ii[0]].donors.append(ii[1])
                aux_dict[ii[0]].atoms.append(ii[1])

            for ii in self.msystem.acceptors:
                if not aux_dict.has_key(ii):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii].resid.type
                    tmp_adnode.name=self.msystem.atom[ii].resid.name
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii].resid.water
                    aux_dict[ii]=tmp_adnode
                aux_dict[ii].acceptors.append(ii)
                aux_dict[ii].atoms.append(ii)

            sel_ions=self.msystem.selection('ion')
            for ii in sel_ions:
                tmp_adnode=adnode()
                tmp_adnode.type=self.msystem.atom[ii].resid.type
                tmp_adnode.name=self.msystem.atom[ii].resid.name
                aux_dict[ii]=tmp_adnode
                aux_dict[ii].atoms.append(ii)


            # symm in ASP,GLU and Terminals
            con_ASP =self.msystem.selection_covalent_chains(['OD1','CG','OD2'],'protein')
            con_GLU =self.msystem.selection_covalent_chains(['OE1','CD','OE2'],'protein')
            con_Term=self.msystem.selection_covalent_chains(['OC1','C','OC2'],'protein')
            for ii in con_ASP:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])
                aux_dict[ii[0]].atoms.append(bb.atoms[0])
            for ii in con_GLU:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])
                aux_dict[ii[0]].atoms.append(bb.atoms[0])
            for ii in con_Term:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])
                aux_dict[ii[0]].atoms.append(bb.atoms[0])

            # head of lipid AOT

            con_head=self.msystem.selection_covalent_chains(['OS1','S','OS2'],'lipid')
            for ii in con_head:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])
                aux_dict[ii[0]].atoms.append(bb.atoms[0])
            con_head=self.msystem.selection_covalent_chains(['OS1','S','OS3'],'lipid')
            for ii in con_head:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])
                aux_dict[ii[0]].atoms.append(bb.atoms[0])

            aux_keys=aux_dict.keys()
            aux_keys.sort()

            for ii in aux_keys:
                aux_node=aux_dict[ii]
                aux_node.atoms.sort()
                self.node.append(aux_node)
            
            del(aux_dict,aux_keys,con_ASP,con_GLU,con_Term)


    def reset(self):

        for node in self.node:
            node.shell1st.don       =[ [] for ii in range(node.don_num)]
            node.shell1st.acc       =[ [] for ii in range(node.acc_num)]
            node.shell1st.bond      =[ [] for ii in range(node.at_num)]
            node.shell1st.don_val   =[ [] for ii in range(node.don_num)]
            node.shell1st.acc_val   =[ [] for ii in range(node.acc_num)]
            node.shell1st.bond_val  =[ [] for ii in range(node.at_num)]
            node.shell1st.don_node  =[ [] for ii in range(node.don_num)]
            node.shell1st.acc_node  =[ [] for ii in range(node.acc_num)]
            node.shell1st.bond_node =[ [] for ii in range(node.at_num)]
            node.shell1st.don_num   =[ 0 for ii in range(node.don_num)]
            node.shell1st.acc_num   =[ 0 for ii in range(node.acc_num)]
            node.shell1st.bond_num  =[ 0 for ii in range(node.at_num)]
            node.mss                =[]
            node.mss_ind            =[]


    def build_shell1st(self,hbonds=None,bonds=None,hbtype='R(o,o)-Ang(o,o,h)',btype='dists'):

        self.hbtype=hbtype
        self.btype=btype

        self.reset()

        if self.hbtype in ['R(o,o)-Ang(o,o,h)','R(o,h)']:
            rever=False
        elif self.hbtype in ['Skinner']:
            rever=True

        filt_aux_don=numpy.zeros((self.don_num),dtype=int)
        filt_aux_acc=numpy.zeros((self.acc_num),dtype=int)
        filt_aux_at=numpy.zeros((self.at_num),dtype=int)

        if hbonds:

            for hb_ind,hb_val in zip(hbonds[0],hbonds[1]):
                atdon=hb_ind[1]
                atacc=hb_ind[2]
                idon=self.don2node[atdon]
                iacc=self.acc2node[atacc]
                self.node[idon[0]].shell1st.don[idon[1]].append(atacc)
                self.node[iacc[0]].shell1st.acc[iacc[1]].append(atdon)
                self.node[idon[0]].shell1st.don_node[idon[1]].append(iacc)
                self.node[iacc[0]].shell1st.acc_node[iacc[1]].append(idon)
                self.node[idon[0]].shell1st.don_val[idon[1]].append(hb_val)
                self.node[iacc[0]].shell1st.acc_val[iacc[1]].append(hb_val)
                self.node[idon[0]].shell1st.don_num[idon[1]]+=1
                self.node[iacc[0]].shell1st.acc_num[iacc[1]]+=1
                filt_aux_don[self.__dict_aux_don__[atdon]]+=1
                filt_aux_acc[self.__dict_aux_acc__[atacc]]+=1

            for ii in numpy.nonzero(filt_aux_don>1)[0]:
                idon=self.don2node[self.don_list[ii]]
                ishell1st=self.node[idon[0]].shell1st
                tups = zip(ishell1st.don_val[idon[1]], ishell1st.don[idon[1]], ishell1st.don_node[idon[1]])
                tups.sort(reverse=rever) 
                [ishell1st.don_val[idon[1]], ishell1st.don[idon[1]], ishell1st.don_node[idon[1]]]=zip(*tups)
                ishell1st.don_val[idon[1]] = list(ishell1st.don_val[idon[1]])
                ishell1st.don[idon[1]]     = list(ishell1st.don[idon[1]])
                ishell1st.don_node[idon[1]]= list(ishell1st.don_node[idon[1]])

            for ii in numpy.nonzero(filt_aux_acc>1)[0]:
                iacc=self.acc2node[self.acc_list[ii]]
                ishell1st=self.node[iacc[0]].shell1st
                tups = zip(ishell1st.acc_val[iacc[1]], ishell1st.acc[iacc[1]], ishell1st.acc_node[iacc[1]])
                tups.sort(reverse=rever)
                [ishell1st.acc_val[iacc[1]], ishell1st.acc[iacc[1]], ishell1st.acc_node[iacc[1]]]=zip(*tups)
                ishell1st.acc_val[iacc[1]] = list(ishell1st.acc_val[iacc[1]])
                ishell1st.acc[iacc[1]]     = list(ishell1st.acc[iacc[1]])
                ishell1st.acc_node[iacc[1]]= list(ishell1st.acc_node[iacc[1]])

        if bonds:

            if self.btype in ['dists']:
                rever=False

            for bond_ind,bond_val in zip(bonds[0],bonds[1]):
                ata=bond_ind[0]
                atb=bond_ind[1]
                ia=self.at2node[ata]
                ib=self.at2node[atb]
                self.node[ia[0]].shell1st.bond[ia[1]].append(atb)
                self.node[ib[0]].shell1st.bond[ib[1]].append(ata)
                self.node[ia[0]].shell1st.bond_node[ia[1]].append(ib)
                self.node[ib[0]].shell1st.bond_node[ib[1]].append(ia)
                self.node[ia[0]].shell1st.bond_val[ia[1]].append(bond_val)
                self.node[ib[0]].shell1st.bond_val[ib[1]].append(bond_val)
                self.node[ia[0]].shell1st.bond_num[ia[1]]+=1
                self.node[ib[0]].shell1st.bond_num[ib[1]]+=1
                filt_aux_at[self.__dict_aux_at__[ata]]+=1
                filt_aux_at[self.__dict_aux_at__[atb]]+=1

            for ii in numpy.nonzero(filt_aux_at>1)[0]:
                iat=self.at2node[self.at_list[ii]]
                ishell1st=self.node[iat[0]].shell1st
                tups = zip(ishell1st.bond_val[iat[1]], ishell1st.bond[iat[1]], ishell1st.bond_node[iat[1]])
                tups.sort(reverse=rever)
                [ishell1st.bond_val[iat[1]], ishell1st.bond[iat[1]], ishell1st.bond_node[iat[1]]]=zip(*tups)
                ishell1st.bond_val[iat[1]] = list(ishell1st.bond_val[iat[1]])
                ishell1st.bond[iat[1]]     = list(ishell1st.bond[iat[1]])
                ishell1st.bond_node[iat[1]]= list(ishell1st.bond_node[iat[1]])
        
        del(filt_aux_don,filt_aux_acc,filt_aux_at)
        

    def build_mss_ind_shell1st(self):

        for node in self.node:
            mss_ind=[]
            inddon=[]
            indacc=[]
            mss_ind.extend([node.don_num,node.acc_num,node.at_num])
            mss_ind.extend(node.shell1st.don_num)
            mss_ind.extend(node.shell1st.acc_num)
            mss_ind.extend(node.shell1st.bond_num)
            for ii in range(node.don_num):
                for jj in range(node.shell1st.don_num[ii]):
                    mss_ind.append(node.shell1st.don_node[ii][jj][0])
            for ii in range(node.acc_num):
                for jj in range(node.shell1st.acc_num[ii]):
                    mss_ind.append(node.shell1st.acc_node[ii][jj][0])
            for ii in range(node.at_num):
                for jj in range(node.shell1st.bond_num[ii]):
                    mss_ind.append(node.shell1st.bond_node[ii][jj][0])
            node.shell1st.mss_ind=mss_ind
            
            node.shell1st.mss_ind_symm=mss_ind
            for eqats in node.filt_ats_symm:
                
        

    def symmetrize_mss_ind_shell1st(self,symm_type=None):

        if symm_type==1: # microstate only for water
            
            # H1 and H2 distinguishable in first shell
            # H1 and H2 indistinguishable in second shell
            # Nodes water indistinguishable
            # Nodes ion indistinguishable
            # Nodes lipid indistinguishable with different peptides ??

            for ii in self.list_water:
                ishell1st=self.node[ii].shell1st
                ishell1st.mss_ind_sym=ishell1st.mss_ind[3:7]
                ishell1st.mss_sym=ishell1st.mss_ind[3:7]
                numhb=sum(ishell1st.mss_ind_sym[0:3])
                numb=ishell1st.mss_ind_sym[3]
                wind=0
                iind=0
                aux_dict={}
                for kk in ishell1st.mss_ind[9:(9+numhb)]:
                    if not aux_dict.has_key(kk):
                        if self.node[kk].type=='Water':
                            aux_dict[kk]='w'+str(wind)
                            wind+=1
                        else:
                            aux_dict[kk]=kk
                    ishell1st.mss_ind_sym.append(kk)
                    ishell1st.mss_sym.append(aux_dict[kk])
                for kk in ishell1st.mss_ind[(9+numhb):(9+numhb+numb)]:
                    if not aux_dict.has_key(kk):
                        if self.node[kk].type=='Ion':
                            aux_dict[kk]='i'+str(iind)
                            iind+=1
                        else:
                            aux_dict[kk]=kk
                    ishell1st.mss_ind_sym.append(kk)
                    ishell1st.mss_sym.append(aux_dict[kk])
                print aux_dict
            del(aux_dict)

    def build_mss_ind(self):

        self.build_mss_ind_shell1st()

        for node in self.node:
            node.mss_ind.append(node.index)
            node.mss_ind.extend(node.shell1st.mss_ind)
            ii=sum(node.shell1st.mss_ind[0:3])+3
            jj=sum(node.shell1st.mss_ind[3:ii])+ii
            for kk in node.shell1st.mss_ind[ii:jj]:
                node.mss_ind.extend(self.node[kk].shell1st.mss_ind)
            

