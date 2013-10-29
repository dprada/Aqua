# tendria que pensar si incluyo los centros de los aromaticos para hacer hbonds

class adnode():

    def __init__(self):

        self.type=None
        self.name=None
        self.index=None
        self.acceptors=[]
        self.donors=[]
        self.num_acc=0
        self.num_don=0
        self.shell1st_acc=[]
        self.shell1st_don=[]

class mss():

    def __init__(self,msystem=None,sets='chains',symm_acc=True,symm_don=True,verbose=True):

        self.type=sets   #'chains','residues','molecules'
        self.msystem=msystem
        self.node=[]
        self.num_nodes=0

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
        for ii in range(self.num_nodes):
            node=self.node[ii]
            node.num_don=len(node.donors)
            node.num_acc=len(node.acceptors)
            for jj in range(node.num_don):
                self.at2node[node.donors[jj]]=[ii,'donor',jj]
                node.shell1st_don.append([])
            for jj in range(node.num_acc):
                self.at2node[node.acceptors[jj]]=[ii,'acceptor',jj]
                node.shell1st_acc.append([])

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

    def build_shell1st(self,hbonds=None,bonds=None):

        for hb in hbonds[0]:
            atdon=hb[1]
            atacc=hb[2]
            idon=self.at2node[atdon]
            iacc=self.at2node[atacc]
            self.node[idon[0]].shell1st_don[idon[2]].append(atacc)
            self.node[iacc[0]].shell1st_acc[iacc[2]].append(atdon)
