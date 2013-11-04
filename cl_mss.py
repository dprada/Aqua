# tendria que pensar si incluyo los centros de los aromaticos para hacer hbonds
import numpy

class shell1st():

    def __init__(self):

        self.acc=[]
        self.don=[]
        self.acc_node=[]
        self.don_node=[]
        self.acc_val=[]
        self.don_val=[]
        self.acc_num=[]
        self.don_num=[]
        self.mss=[]
        self.mss_ind=[]


class adnode():

    def __init__(self):

        self.type=None
        self.name=None
        self.index=None
        self.acceptors=[]
        self.donors=[]
        self.acc_num=0
        self.don_num=0
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
            node.index=ii
            node.don_num=len(node.donors)
            node.acc_num=len(node.acceptors)
            for jj in range(node.don_num):
                self.at2node[node.donors[jj]]=[ii,'donor',jj]
                self.don2node[node.donors[jj]]=[ii,jj]
                node.shell1st.don.append([])
                node.shell1st.don_val.append([])
                node.shell1st.don_num.append(0)
            for jj in range(node.acc_num):
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
                    #tmp_adnode.index=self.msystem.atom[ii[0]].resid.index
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii[0]].resid.water
                    aux_dict[ii[0]]=tmp_adnode
                aux_dict[ii[0]].donors.append(ii[1])

            for ii in self.msystem.acceptors:
                if not aux_dict.has_key(ii):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii].resid.type
                    tmp_adnode.name=self.msystem.atom[ii].resid.name
                    #tmp_adnode.index=self.msystem.atom[ii].resid.index
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
                    #tmp_adnode.index=self.msystem.atom[ii[0]].resid.index
                    if tmp_adnode.type=='Water':
                        tmp_adnode.water=self.msystem.atom[ii[0]].resid.water
                    aux_dict[ii[0]]=tmp_adnode
                aux_dict[ii[0]].donors.append(ii[1])

            for ii in self.msystem.acceptors:
                if not aux_dict.has_key(ii):
                    tmp_adnode=adnode()
                    tmp_adnode.type=self.msystem.atom[ii].resid.type
                    tmp_adnode.name=self.msystem.atom[ii].resid.name
                    #tmp_adnode.index=self.msystem.atom[ii].resid.index
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
            node.shell1st.don=[ [] for ii in range(node.don_num)]
            node.shell1st.acc=[ [] for ii in range(node.acc_num)]
            node.shell1st.don_val=[ [] for ii in range(node.don_num)]
            node.shell1st.acc_val=[ [] for ii in range(node.acc_num)]
            node.shell1st.don_node=[ [] for ii in range(node.don_num)]
            node.shell1st.acc_node=[ [] for ii in range(node.acc_num)]
            node.shell1st.don_num=[ 0 for ii in range(node.don_num)]
            node.shell1st.acc_num=[ 0 for ii in range(node.acc_num)]
            node.mss=[]
            node.mss_ind=[]


    def build_shell1st(self,hbonds=None,bonds=None,hbtype='R(o,o)-Ang(o,o,h)',btype=None):

        self.hbtype=hbtype
        self.btype=btype

        self.reset()

        if self.hbtype in ['R(o,o)-Ang(o,o,h)','R(o,h)']:
            rever=False
        elif self.hbtype in ['Skinner']:
            rever=True

        filt_aux_don=numpy.zeros((self.don_num),dtype=int)
        filt_aux_acc=numpy.zeros((self.acc_num),dtype=int)

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
         
        for ii in numpy.nonzero(filt_aux_acc>1)[0]:
            iacc=self.acc2node[self.acc_list[ii]]
            ishell1st=self.node[iacc[0]].shell1st
            tups = zip(ishell1st.acc_val[iacc[1]], ishell1st.acc[iacc[1]], ishell1st.acc_node[iacc[1]])
            tups.sort(reverse=rever)
            [ishell1st.acc_val[iacc[1]], ishell1st.acc[iacc[1]], ishell1st.acc_node[iacc[1]]]=zip(*tups)

        del(filt_aux_don,filt_aux_acc)

        # is this faster?
        #for node in self.node:
        #    for ii in range(node.acc_num):
        #        if node.shell1st.acc_num[ii]>1:
        #            tups = zip(node.shell1st.acc_val[ii], node.shell1st.acc[ii])
        #            tups.sort()
        #            [node.shell1st.acc_val[ii], node.shell1st.acc[ii]]=zip(*tups)
        #    for ii in range(node.don_num):
        #        if node.shell1st.don_num[ii]>1:
        #            tups = zip(node.shell1st.don_val[ii], node.shell1st.don[ii])
        #            tups.sort()
        #            [node.shell1st.don_val[ii], node.shell1st.don[ii]]=zip(*tups)


    def build_mss(self):

        for node in self.node:
            mss=[]
            inddon=[]
            indacc=[]
            mss.extend([node.don_num,node.acc_num])
            mss.extend(node.shell1st.don_num)
            mss.extend(node.shell1st.acc_num)
            for ii in range(node.don_num):
                for jj in range(node.shell1st.don_num[ii]):
                    mss.append(node.shell1st.don_node[ii][jj][0])
            for ii in range(node.acc_num):
                for jj in range(node.shell1st.acc_num[ii]):
                    mss.append(node.shell1st.acc_node[ii][jj][0])            
            node.shell1st.mss=mss

        for node in self.node:
            node.mss.append(node.index)
            node.mss.extend(node.shell1st.mss)
            ii=sum(node.shell1st.mss[0:2])+2
            jj=sum(node.shell1st.mss[2:ii])+ii
            for kk in node.shell1st.mss[ii:jj]:
                node.mss.extend(self.node[kk].shell1st.mss)
            

    def build_mss_slow_antes(self):

        #for node in self.node:
        #    if node.type=='Water':   # Puedo hacer una lista previa con las aguas y otra con la proteina
        #        mss=numpy.zeros((4),dtype=int,order='Fortran')
        #        mss_ind=numpy.zeros((4),dtype=int,order='Fortran')
        #        if node.shell1st.don_num[0]:
        #            mss_ind[0]=node.shell1st.don[0][0]; mss[0]=2
        #        if node.shell1st.don_num[1]:
        #            mss_ind[1]=node.shell1st.don[1][0]; mss[1]=3
        #        if node.shell1st.acc_num[0]==1:
        #            mss_ind[2]=node.shell1st.acc[0][0]; mss[2]=4
        #        elif node.shell1st.acc_num[0]==2:
        #            mss_ind[2]=node.shell1st.acc[0][0]; mss[2]=4
        #            mss_ind[3]=node.shell1st.acc[0][1]; mss[3]=5
        #        node.shell1st.mss=mss; node.shell1st.mss_ind=mss_ind 
                
        #for node in self.node:
        #    if node.type=='Water':
        #        mss_ind=numpy.zeros((16),dtype=int)
        #        mss_filt=numpy.zeros((16),dtype=bool)
        #        if node.shell1st.don_num[0]:
        #            ii=self.acc2node[node.shell1st.don[0][0]]
        #            nodeh1=self.node[ii[0]]
        #            mss_filt[0]=True
        #            mss_ind[0]=ii[0]
        #            if nodeh1.type=='Water':
        #                if nodeh1.shell1st.don_num[0]:
        #                    jj=self.acc2node[nodeh1.shell1st.don[0][0]]
        #                    mss_filt[4]=True
        #                    mss_ind[4]=jj[0]
        #                if nodeh1.shell1st.don_num[1]:
        #                    jj=self.acc2node[nodeh1.shell1st.don[1][0]]
        #                    mss_filt[5]=True
        #                    mss_ind[5]=jj[0]
        #                if nodeh1.shell1st.acc_num[0]>1:
        #                    kk=self.don2node[nodeh1.shell1st.acc[0][0]]
        #                    ll=self.don2node[nodeh1.shell1st.acc[0][1]]
        #                    if kk==node.index or ll==node.index:
        #                        if kk==node.index:
        #                            mss_filt[6]=True
        #                            mss_ind[6]=ll
        #                        else:
        #                            mss_filt[6]=True
        #                            mss_ind[6]=kk
        #                    else:
        #                        mss_filt[6]=True
        #                        mss_ind[6]=kk
        #        if node.shell1st.don_num[1]:
        #            ii=self.acc2node[node.shell1st.don[1][0]]
        #            nodeh2=self.node[ii[0]]
        #            mss_filt[0]=True
        #            mss_ind[0]=ii[0]
        #            if nodeh2.type=='Water':
        #                if nodeh2.shell1st.don_num[0]:
        #                    jj=self.acc2node[nodeh2.shell1st.don[0][0]]
        #                    mss_filt[4]=True
        #                    mss_ind[4]=jj[0]
        #                if nodeh2.shell1st.don_num[1]:
        #                    jj=self.acc2node[nodeh2.shell1st.don[1][0]]
        #                    mss_filt[5]=True
        #                    mss_ind[5]=jj[0]
        #                if nodeh2.shell1st.acc_num[0]>1:
        #                    kk=self.don2node[nodeh2.shell1st.acc[0][0]]
        #                    ll=self.don2node[nodeh2.shell1st.acc[0][1]]
        #                    if kk==node.index or ll==node.index:
        #                        if kk==node.index:
        #                            mss_filt[6]=True
        #                            mss_ind[6]=ll
        #                        else:
        #                            mss_filt[6]=True
        #                            mss_ind[6]=kk
        #                    else:
        #                        mss_filt[6]=True
        #                        mss_ind[6]=kk
        #        if node.shell1st
        #                        
        #                    
        #                
        # 
        # 
        #        node.mss[0:4]=node.shell1st.mss[:]
        #        node.mss_ind[0:4]=node.shell1st.mss_ind[:]
        #        if node.mss[0]:
        #            jj.node.acceptor[0]
        #            ii=self.at2node[node.mss_ind[0]]
        #            node.mss_ind[4:6]=
        pass
