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
        self.mss_ind_atoms=[]
        self.mss_ind_nodes=[]
        self.new_symm=[]

class node():

    def __init__(self,name=None,type=None):

        self.index=None
        self.name=name
        self.label=None
        self.type=type
        self.category=[0,0,0]  # code_node, code_set, code_inside_set

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
        self.symm_node=[]
        #self.new_order=None
        #self.symm_broken_ats=True
        #self.new_symm_ats=[]


    def add_atom(self,index=None,donor=False,acceptor=False,nonpolar=False):

        self.atom[index]=atom(index=index,donor=donor,acceptor=acceptor,nonpolar=nonpolar)

    def merge_node(self,node=None,update_ivars=True):
     
        for ii,jj in node.atom.iteritems():
            self.atom[ii]=jj

        if update_ivars:
            pass
            #Cuidado que hay que arreglar las variables internas

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

class internal_vars():

    def __init__(self,num_atoms=None,num_nodes=None):

        self.atom2id={}
        self.node2id={}
        self.categories=numpy.zeros((num_nodes,3),dtype=int,order="Fortran")
        self.atomid2nodeid=numpy.zeros((num_atoms),dtype=int,order="Fortran")
        self.max_atoms_node=0
        self.sets=[]
        self.num_sets=0
        self.nodes_per_set=[]
        self.num_sets_sets=[]
        self.num_categories=0

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
        self.symm_sets_nodes=symm_sets_nodes
        self.symm_sets_nodes_aux=[]

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
                    aa=numpy.in1d(node.atoms,criterium)
                    if aa.sum():
                        node.symm_ats.append(aa)
            del(symm_ats_list)

        ## symmetric nodes
        count_category=0
        if self.symm_nodes:
            symm_nodes_list=[]
            if type(self.symm_nodes) in [list,tuple]:
                for sel in self.symm_nodes:
                    symm_nodes_list.append(self.msystem.selection(sel))
            else:
                symm_nodes_list.append(self.msystem.selection(sel))
            for node in self.node:
                for ii in range(len(symm_nodes_list)):
                    if True in numpy.in1d(node.atoms,symm_nodes_list[ii]):
                        node.category[0]=count_category+ii
            count_category+=len(symm_nodes_list)

        ## symmetric sets of nodes

        aux_sets=[]
        aux_category_sets=[]
        if self.symm_sets_nodes:
            for criterium in self.symm_sets_nodes:
                if criterium == 'lipids':
                    symm_sets_nodes=self.msystem.lipid
                elif criterium == 'proteins':
                    symm_sets_nodes=self.msystem.protein
                aux_category_sets.append(count_category)
                aux=[item for sublist in symm_sets_nodes for item in sublist]
                num_sets=len(symm_sets_nodes)
                sets={ii:[] for ii in range(num_sets)}
                for node in self.node:
                    if True in numpy.in1d(node.atoms,aux):
                        node.category[0]=count_category
                        for ii in range(num_sets):
                            if True in numpy.in1d(node.atoms,symm_sets_nodes[ii]):
                                sets[ii].append(node.index)
                count_category+=1
                aux_sets.append(sets.values())
                for ii in range(num_sets):
                    aux=sets[ii]
                    for jj in range(len(aux)):
                        self.node[aux[jj]].category[1]=ii
                        self.node[aux[jj]].category[2]=jj


        ### building support

        self.ivars=internal_vars(self.num_atoms,self.num_nodes)

        self.ivars.num_categories=count_category
        self.ivars.sets=aux_sets
        self.ivars.category_sets=numpy.array(aux_category_sets)
        self.ivars.num_sets=len(self.ivars.sets)
        for ii in range(self.ivars.num_sets):
            self.ivars.nodes_per_set.append(len(self.ivars.sets[ii][0]))
            self.ivars.num_sets_sets.append(len(self.ivars.sets[ii]))

        ii_n=1
        ii_a=1
        for node in self.node:
            self.ivars.node2id[node.index]=ii_n
            self.ivars.categories[ii_n-1,:]=node.category
            for atom in node.atom.values():
                self.ivars.atom2id[atom.index]=ii_a
                self.ivars.atomid2nodeid[ii_a-1]=ii_n
                ii_a+=1
            ii_n+=1
            if self.ivars.max_atoms_node<node.num_atoms:
                self.ivars.max_atoms_node=node.num_atoms

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
 
 
        ## build support

        tot_num_hbs=len(hbonds[0])*2
        tot_num_bs=len(bonds[0])*2

        mss_funcs.support_up(tot_num_hbs,tot_num_bs)

        for node in self.node:
            for atom in node.atom.values():
                jj=self.ivars.atom2id[atom.index]
                for xx in atom.hbond:
                    ii=self.ivars.atom2id[xx]
                    mss_funcs.add_hb_support(jj,ii)
                for xx in atom.bond:
                    ii=self.ivars.atom2id[xx]
                    mss_funcs.add_b_support(jj,ii)

        codes_atom,codes_node=mss_funcs.support_down(self.ivars.atomid2nodeid,self.ivars.categories,self.ivars.num_sets_sets,
                                                     self.ivars.nodes_per_set,self.ivars.num_categories,self.num_atoms,
                                                     self.num_nodes,self.ivars.num_sets)

        return codes_atom,codes_node


#    def build_mss_shell1st(self):
# 
# 
#        # Hasta aqui el soporte. 
# 
#        otro_auxilio_hb=numpy.zeros((self.num_nodes,self.max_ats_node),dtype=int)
#        otro_auxilio_b=numpy.zeros((self.num_nodes,self.max_ats_node),dtype=int)
#        for node in self.node: # Para el orden de los mismos atomos del nodo. P.e.: H1 y H2
#            order=copy.copy(node.atoms)
#            node.shell1st.new_symm.extend(['X'])
#            if node.symm_ats:
#                support=numpy.zeros((node.num_atoms,40),dtype=int,order='Fortran')
#                for ii in range(node.num_atoms):
#                    jj=order[ii]
#                    atom=node.atom[jj]
#                    support[ii,0:10]=atom.suphb
#                    support[ii,10:20]=atom.supb
#                    support[ii,20:30]=atom.suphb_2sh
#                    support[ii,30:40]=atom.supb_2sh
#                broken=0
#                for criterium in node.symm_ats:
#                    order,new_symm=mss_funcs.breaking_symmetry_1st(criterium,order,support,node.num_atoms,40)
#                    node.new_symm_ats.append(new_symm)
#                    broken+=sum(new_symm)
#                node.shell1st.new_symm.extend(new_symm)
#            else:
#                node.shell1st.new_symm.extend([0 for ii in range(node.num_atoms)])
#            node.shell1st.new_symm.extend(['X' for ii in range(2*node.num_atoms)])
#            mss_ind_atoms=[]
#            mss_ind_nodes=[]
#            mss_ind_atoms.append(node.num_atoms)
#            mss_ind_atoms.extend(order)
#            node.new_order=order
#            mss_ind_nodes.append(node.num_atoms)
#            mss_ind_nodes.extend([node.index for ii in order])
#            kkk=0
#            for ii in order:                                            
#                cc=[node.atom[ii].num_hbonds,node.atom[ii].num_bonds]   
#                mss_ind_atoms.extend(cc)                                
#                mss_ind_nodes.extend(cc)
#                otro_auxilio_hb[node.index,kkk]=node.atom[ii].num_hbonds
#                otro_auxilio_b[node.index,kkk]=node.atom[ii].num_bonds
#                kkk+=1
#            node.shell1st.mss_ind_atoms=mss_ind_atoms
#            node.shell1st.mss_ind_nodes=mss_ind_nodes
# 
#        for node in self.node: # para quitar cierta simetria
#            if 1 in node.shell1st.new_symm:
#                support=numpy.zeros((node.num_atoms,2*self.max_ats_node),dtype=int,order='Fortran')
#                order=node.new_order
#                criterium=node.shell1st.new_symm
#                node.shell1st.new_symm=[]
#                node.shell1st.new_symm.extend(['X'])
#                for ii in range(node.num_atoms):
#                    atom=order[ii]
#                    for jj in node.atom[atom].hbond_node:
#                        support[ii,0:3]+=otro_auxilio_hb[jj,:]
#                        support[ii,3:6]+=otro_auxilio_b[jj,:]
#                    for jj in node.atom[atom].bond_node:
#                        support[ii,0:3]+=otro_auxilio_hb[jj,:]
#                        support[ii,3:6]+=otro_auxilio_b[jj,:]
#                broken=0
#                for criterium in node.symm_ats:
#                    order,new_symm=mss_funcs.breaking_symmetry_1st(criterium,order,support,node.num_atoms,6)
#                    node.new_symm_ats.append(new_symm)
#                    broken+=sum(new_symm)
#                node.shell1st.new_symm.extend(new_symm)
#                node.shell1st.new_symm.extend(['X' for ii in range(2*node.num_atoms)])
#                mss_ind_atoms=[]
#                mss_ind_nodes=[]
#                mss_ind_atoms.append(node.num_atoms)
#                mss_ind_atoms.extend(order)
#                node.new_order=order
#                mss_ind_nodes.append(node.num_atoms)
#                mss_ind_nodes.extend([node.index for ii in order])
#                kkk=0
#                for ii in order:                                            
#                    cc=[node.atom[ii].num_hbonds,node.atom[ii].num_bonds]   
#                    mss_ind_atoms.extend(cc)                                
#                    mss_ind_nodes.extend(cc)
#                    kkk+=1
#                node.shell1st.mss_ind_atoms=mss_ind_atoms
#                node.shell1st.mss_ind_nodes=mss_ind_nodes
#                
# 
#        for node in self.node: # Para el orden de los nodos enlazados.
# 
#            aa=[]
#            bb=[]
#            order=node.new_order
#            if self.symm_nodes:
#                for ii in order:
#                    cc=[node.atom[ii].num_hbonds,node.atom[ii].num_bonds]
#                    order_hb_atoms=node.atom[ii].hbond
#                    order_hb_nodes=node.atom[ii].hbond_node
#                    order_b_atoms =node.atom[ii].bond
#                    order_b_nodes =node.atom[ii].bond_node
#                    broken0=1
#                    broken1=1
#                    if cc[0]>1:
#                        aux_set={}
#                        support=numpy.zeros((cc[0],1+5+10+10+2*self.max_ats_node+1),dtype=int,order='Fortran')
#                        for jj in range(cc[0]):
#                            support[jj,0]=self.node[order_hb_nodes[jj]].codigo
#                            support[jj,6:16]=self.node[order_hb_nodes[jj]].suphb
#                            support[jj,16:26]=self.node[order_hb_nodes[jj]].supb
#                            support[jj,26:29]=otro_auxilio_hb[order_hb_nodes[jj],:]
#                            support[jj,29:32]=otro_auxilio_b[order_hb_nodes[jj],:]
#                            if self.node[order_hb_nodes[jj]].codigo==2:
#                                support[jj,32]=self.node[order_hb_nodes[jj]].codigo_sets[1]
#                                try:
#                                    aux_set[self.node[order_hb_nodes[jj]].codigo_sets[0]].append(jj)
#                                except:
#                                    aux_set[self.node[order_hb_nodes[jj]].codigo_sets[0]]=[jj]
#                        for same_set in aux_set.values():
#                            aaa=numpy.zeros((5),dtype=int,order='Fortran')
#                            bbb=numpy.zeros((26),dtype=int,order='Fortran')
#                            for ww in same_set:
#                                aaa[support[ww,32]]+=1
#                                bbb[:]+=support[ww,6:32]
#                            for ww in same_set:
#                                support[ww,1:6]=aaa[:]
#                                support[ww,6:32]=bbb[:]
#                        order_hb_atoms,order_hb_nodes,new_symm_hb=mss_funcs.breaking_symmetry_2nd(order_hb_atoms,order_hb_nodes,support,cc[0],33)
#                        node.shell1st.new_symm.extend(2*new_symm_hb)
#                    if cc[0]==1:
#                        node.shell1st.new_symm.extend([0])
#                    if cc[1]>1:
#                        aux_set={}
#                        support=numpy.zeros((cc[1],1+5+10+10+2*self.max_ats_node+1),dtype=int,order='Fortran')
#                        for jj in range(cc[1]):
#                            support[jj,0]=self.node[order_b_nodes[jj]].codigo
#                            support[jj,6:16]=self.node[order_b_nodes[jj]].suphb
#                            support[jj,16:26]=self.node[order_b_nodes[jj]].supb
#                            support[jj,26:29]=otro_auxilio_hb[order_b_nodes[jj],:]
#                            support[jj,29:32]=otro_auxilio_b[order_b_nodes[jj],:]
#                            if self.node[order_b_nodes[jj]].codigo==2:
#                                support[jj,32]=self.node[order_b_nodes[jj]].codigo_sets[1]
#                                try:
#                                    aux_set[self.node[order_b_nodes[jj]].codigo_sets[0]].append(jj)
#                                except:
#                                    aux_set[self.node[order_b_nodes[jj]].codigo_sets[0]]=[jj]
#                        for same_set in aux_set.values():
#                            aaa=numpy.zeros((5),dtype=int,order='Fortran')
#                            bbb=numpy.zeros((26),dtype=int,order='Fortran')
#                            for ww in same_set:
#                                aaa[support[ww,32]]+=1
#                                bbb[:]+=support[ww,6:32]
#                            for ww in same_set:
#                                support[ww,1:6]=aaa[:]
#                                support[ww,6:32]=bbb[:]
#                        order_b_atoms,order_b_nodes,new_symm_b=mss_funcs.breaking_symmetry_2nd(order_b_atoms,order_b_nodes,support,cc[1],33)
#                        node.shell1st.new_symm.extend(2*new_symm_b)
#                    if cc[1]==1:
#                        node.shell1st.new_symm.extend([0])
#                    aa.extend(order_hb_atoms)
#                    aa.extend(order_b_atoms)
#                    bb.extend(order_hb_nodes)
#                    bb.extend(order_b_nodes)
#            else:
#                for ii in order:                                           
#                    aa.extend(node.atom[ii].hbond)
#                    aa.extend(node.atom[ii].bond)
#                    bb.extend(node.atom[ii].hbond_node)
#                    bb.extend(node.atom[ii].bond_node)
#                    node.shell1st.new_symm.extend([0 for ii in range(node.atom[ii].num_hbonds)])
#                    node.shell1st.new_symm.extend([0 for ii in range(node.atom[ii].num_bonds)])
#            node.shell1st.mss_ind_atoms.extend(aa)
#            node.shell1st.mss_ind_nodes.extend(bb)
# 
#            pacambiar=[]
#            n_ats=node.shell1st.mss_ind_nodes[0]
#            n_bs=sum(node.shell1st.mss_ind_nodes[(1+n_ats):(1+n_ats*3)])
#            pacambiar.extend(range(1,(1+n_ats)))
#            pacambiar.extend(range((1+n_ats*3),(1+n_ats*3+n_bs)))
# 
#            mss=copy.copy(node.shell1st.mss_ind_nodes)
#            aux2_dict={}
#            aux2_set={}
#            cc_water=0
#            cc_lipid=0
#            cc_ion=0
#            for ii in pacambiar:
#                jj=mss[ii]
#                if aux2_dict.has_key(jj)==False:
#                    if jj in self.waters:
#                        aux2_dict[jj]='w'+str(cc_water)
#                        cc_water+=1
#                    elif jj in self.lipids:
#                        if aux2_set.has_key(self.node[jj].codigo_sets[0])==False:
#                            aux2_set[self.node[jj].codigo_sets[0]]=cc_lipid
#                            cc_lipid+=1
#                        aux2_dict[jj]='l'+str(aux2_set[self.node[jj].codigo_sets[0]])+'-'+str(self.node[jj].codigo_sets[1])
#                    elif jj in self.ions:
#                        aux2_dict[jj]='i'+str(cc_ion)
#                        cc_ion+=1
#                mss[ii]=aux2_dict[jj]
#            node.shell1st.mss=mss
#            
# 
#    def build_mss_shell2nd(self):
# 
#        # Pulo el hecho de que el nodo en el centro rompe asimetrias
#        for node in self.node:
#            node.shell2nd.mss_ind_atoms.extend(node.shell1st.mss_ind_atoms)
#            node.shell2nd.mss_ind_nodes.extend(node.shell1st.mss_ind_nodes)
#            node.shell2nd.new_symm.extend(node.shell1st.new_symm)
#            pacambiar=[]
#            n_ats=node.shell1st.mss_ind_nodes[0]
#            n_bs=sum(node.shell1st.mss_ind_nodes[(1+n_ats):(1+n_ats*3)])
#            pacambiar.extend(range(1,(1+n_ats)))
#            pacambiar.extend(range((1+n_ats*3),(1+n_ats*3+n_bs)))
#            ll=len(node.shell1st.mss_ind_nodes)
#            for ii in range((1+n_ats*3),(1+n_ats*3+n_bs)):
#                jj=node.shell1st.mss_ind_nodes[ii]
#                nn_ats=self.node[jj].shell1st.mss_ind_nodes[0]
#                nn_bs=sum(self.node[jj].shell1st.mss_ind_nodes[(1+nn_ats):(1+nn_ats*3)])
#                pacambiar.extend(range(ll+1,(ll+1+nn_ats)))
#                pacambiar.extend(range((ll+1+nn_ats*3),(ll+1+nn_ats*3+nn_bs)))
#                ll+=len(self.node[jj].shell1st.mss_ind_nodes)
#                node.shell2nd.mss_ind_atoms.extend(self.node[jj].shell1st.mss_ind_atoms)
#                node.shell2nd.mss_ind_nodes.extend(self.node[jj].shell1st.mss_ind_nodes)
#                node.shell2nd.new_symm.extend(self.node[jj].shell1st.new_symm)
#            
#            mss=copy.copy(node.shell2nd.mss_ind_nodes)
#            aux2_dict={}
#            aux2_set={}
#            cc_water=0
#            cc_lipid=0
#            cc_ion=0
#            for ii in pacambiar:
#                jj=mss[ii]
#                if aux2_dict.has_key(jj)==False:
#                    if jj in self.waters:
#                        aux2_dict[jj]='w'+str(cc_water)
#                        cc_water+=1
#                    elif jj in self.lipids:
#                        if aux2_set.has_key(self.node[jj].codigo_sets[0])==False:
#                            aux2_set[self.node[jj].codigo_sets[0]]=cc_lipid
#                            cc_lipid+=1
#                        aux2_dict[jj]='l'+str(aux2_set[self.node[jj].codigo_sets[0]])+'-'+str(self.node[jj].codigo_sets[1])
#                    elif jj in self.ions:
#                        aux2_dict[jj]='i'+str(cc_ion)
#                        cc_ion+=1
#                mss[ii]=aux2_dict[jj]
#            node.shell2nd.mss=mss
# 
# 
#    def breaking_symmetry_centrality(self,node=None,center=None):
# 
#        new_mss_ind_atoms=copy.copy(self.node[jj].shell1st.mss_ind_atoms)
#        new_mss_ind_nodes=copy.copy(self.node[jj].shell1st.mss_ind_atoms)
#        new_new_symm     =copy.copy(self.node[jj].shell1st.new_symm)

        
