# tendria que pensar si incluyo los centros de los aromaticos para hacer hbonds
import numpy

class shell1st():

    def __init__(self):

        self.acc=[]
        self.don=[]
        self.acc_val=[]
        self.don_val=[]
        self.acc_num=[]
        self.don_num=[]

class adnode():

    def __init__(self):

        self.type=None
        self.name=None
        self.index=None
        self.acceptors=[]
        self.donors=[]
        self.num_acc=0
        self.num_don=0
        self.shell1st=shell1st()
        self.mss=[]
        self.mss_ind=[]

class mss():

    def __init__(self,msystem=None,sets='chains',symm_acc=True,symm_don=True,verbose=True):

        self.type=sets   #'chains','residues','molecules'
        self.msystem=msystem
        self.node=[]
        self.num_nodes=0
        self.hbtype=None
        self.btype=None

        self.build_nodes()
        self.num_nodes=len(self.node)
        
        self.list_protein=[]
        self.list_water=[]
        for ii in range(self.num_nodes):
            if self.node[ii].type=='Protein':
                self.list_protein.append(ii)
            elif self.node[ii].type=='Water':
                self.list_water.append(ii)

        self.at2node={}
        self.acc2node={}
        self.don2node={}

        for ii in range(self.num_nodes):
            node=self.node[ii]
            node.num_don=len(node.donors)
            node.num_acc=len(node.acceptors)
            for jj in range(node.num_don):
                self.at2node[node.donors[jj]]=[ii,'donor',jj]
                self.don2node[node.donors[jj]]=[ii,jj]
                node.shell1st.don.append([])
                node.shell1st.don_val.append([])
                node.shell1st.don_num.append(0)
            for jj in range(node.num_acc):
                self.at2node[node.acceptors[jj]]=[ii,'acceptor',jj]
                self.acc2node[node.acceptors[jj]]=[ii,jj]
                node.shell1st.acc.append([])
                node.shell1st.acc_num.append(0)

        self.at_list=self.at2node.keys()
        self.acc_list=self.acc2node.keys()
        self.don_list=self.don2node.keys()
        self.acc_num=len(self.acc_list)
        self.don_num=len(self.don_list)

        self.__dict_aux_don__={}
        self.__dict_aux_acc__={}

        jj=0
        for ii in self.don_list:
            self.__dict_aux_don__[ii]=jj
            jj+=1

        jj=0
        for ii in self.acc_list:
            self.__dict_aux_acc__[ii]=jj
            jj+=1


        if verbose:
            self.info()

    def info(self):

        print '#',self.num_nodes,'nodes'

    def build_nodes(self):

        if self.type=='chains':
            
            aux_dict={}
            aa=0
            for ii in self.msystem.donors:
                if not aux_dict.has_key(ii[0]):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii[0]].resid.type
                    tmp_adnode.name=self.msystem.atom[ii[0]].resid.name
                    tmp_adnode.index=self.msystem.atom[ii[0]].resid.index
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii[0]].resid.water
                    aux_dict[ii[0]]=tmp_adnode
                aux_dict[ii[0]].donors.append(ii[1])

            for ii in self.msystem.acceptors:
                if not aux_dict.has_key(ii):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii].resid.type
                    tmp_adnode.name=self.msystem.atom[ii].resid.name
                    tmp_adnode.index=self.msystem.atom[ii].resid.index
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii].resid.water
                    aux_dict[ii]=tmp_adnode
                aux_dict[ii].acceptors.append(ii)

            aux_keys=aux_dict.keys()
            aux_keys.sort()

            for ii in aux_keys:
                self.node.append(aux_dict[ii])
            
            del(aux_dict,aux_keys)

        if self.type=='chains+COn':
            
            aux_dict={}
            aa=0
            for ii in self.msystem.donors:
                if not aux_dict.has_key(ii[0]):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii[0]].resid.type
                    tmp_adnode.name=self.msystem.atom[ii[0]].resid.name
                    tmp_adnode.index=self.msystem.atom[ii[0]].resid.index
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii[0]].resid.water
                    aux_dict[ii[0]]=tmp_adnode
                aux_dict[ii[0]].donors.append(ii[1])

            for ii in self.msystem.acceptors:
                if not aux_dict.has_key(ii):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii].resid.type
                    tmp_adnode.name=self.msystem.atom[ii].resid.name
                    tmp_adnode.index=self.msystem.atom[ii].resid.index
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii].resid.water
                    aux_dict[ii]=tmp_adnode
                aux_dict[ii].acceptors.append(ii)

            # symm in ASP,GLU and Terminals
            con_ASP =self.msystem.selection_covalent_chains(['OD1','CG','OD2'])
            con_GLU =self.msystem.selection_covalent_chains(['OE1','CD','OE2'])
            con_Term=self.msystem.selection_covalent_chains(['OC1','C','OC2'])
            for ii in con_ASP:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])
            for ii in con_GLU:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])
            for ii in con_Term:
                bb=aux_dict.pop(ii[2])
                aux_dict[ii[0]].acceptors.append(bb.acceptors[0])

            aux_keys=aux_dict.keys()
            aux_keys.sort()

            for ii in aux_keys:
                self.node.append(aux_dict[ii])
            
            del(aux_dict,aux_keys,con_ASP,con_GLU,con_Term)

    def reset(self):

        for node in self.node:
            node.shell1st.don=[ [] for ii in range(node.num_don)]
            node.shell1st.acc=[ [] for ii in range(node.num_acc)]
            node.shell1st.don_val=[ [] for ii in range(node.num_don)]
            node.shell1st.acc_val=[ [] for ii in range(node.num_acc)]
            node.shell1st.don_num=[ 0 for ii in range(node.num_don)]
            node.shell1st.acc_num=[ 0 for ii in range(node.num_acc)]
            node.mss=[]
            node.mss_ind=[]


    def build_shell1st(self,hbonds=None,bonds=None,hbtype='R(o,o)-Ang(o,o,h)',btype=None):

        self.hbtype=hbtype
        self.btype=btype

        self.reset()

        filt_aux_don=numpy.zeros((self.don_num),dtype=int)
        filt_aux_acc=numpy.zeros((self.acc_num),dtype=int)

        for hb_ind,hb_val in zip(hbonds[0],hbonds[1]):
            atdon=hb_ind[1]
            atacc=hb_ind[2]
            idon=self.don2node[atdon]
            iacc=self.acc2node[atacc]
            self.node[idon[0]].shell1st.don[idon[1]].append(atacc)
            self.node[iacc[0]].shell1st.acc[iacc[1]].append(atdon)
            self.node[idon[0]].shell1st.don_val[idon[1]].append(hb_val)
            self.node[iacc[0]].shell1st.acc_val[iacc[1]].append(hb_val)
            self.node[idon[0]].shell1st.don_num[idon[1]]+=1
            self.node[iacc[0]].shell1st.acc_num[iacc[1]]+=1
            filt_aux_don[self.__dict_aux_don__[atdon]]+=1
            filt_aux_acc[self.__dict_aux_acc__[atacc]]+=1

        for ii in numpy.nonzero(filt_aux_don>1)[0]:
            idon=self.don2node[self.don_list[ii]]
            tups = zip(self.node[idon[0]].shell1st.don_val[idon[1]], self.node[idon[0]].shell1st.don[idon[1]])
            tups.sort() 
            [self.node[idon[0]].shell1st.don_val[idon[1]], self.node[idon[0]].shell1st.don[idon[1]]]=zip(*tups)
         
        for ii in numpy.nonzero(filt_aux_acc>1)[0]:
            iacc=self.acc2node[self.acc_list[ii]]
            tups = zip(self.node[iacc[0]].shell1st.acc_val[iacc[1]], self.node[iacc[0]].shell1st.acc[iacc[1]])
            tups.sort()
            [self.node[iacc[0]].shell1st.acc_val[iacc[1]], self.node[iacc[0]].shell1st.acc[iacc[1]]]=zip(*tups)

        del(filt_aux_don,filt_aux_acc)

        # is this faster?
        #for node in self.node:
        #    for ii in range(node.num_acc):
        #        if node.shell1st.acc_num[ii]>1:
        #            tups = zip(node.shell1st.acc_val[ii], node.shell1st.acc[ii])
        #            tups.sort()
        #            [node.shell1st.acc_val[ii], node.shell1st.acc[ii]]=zip(*tups)
        #    for ii in range(node.num_don):
        #        if node.shell1st.don_num[ii]>1:
        #            tups = zip(node.shell1st.don_val[ii], node.shell1st.don[ii])
        #            tups.sort()
        #            [node.shell1st.don_val[ii], node.shell1st.don[ii]]=zip(*tups)


    def build_mss(self):

        if self.hbtype in ['R(o,o)-Ang(o,o,h)']:
            reverse=False
        elif self.hbtype in ['Skinner']:
            reverse=True

        for node in self.node:
            if node.type=='Water':   # Puedo hacer una lista previa con las aguas y otra con la proteina
                node.mss_ind=numpy.zeros((16),dtype=int,order='Fortran')
                #if 
        
