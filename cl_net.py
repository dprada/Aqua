from pyn_fort_net import glob as f_net
from pyn_fort_kin_anal import glob as f_kin_anal
from pyn_fort_anal_trajs import glob as f_trajs
import pyn_math as pyn_math
import numpy
import copy
from os import system

#####################################################################################
##### Networks
#####################################################################################

class cl_io():

    def __init__(self):

        self.v={}
        self.reset_H1()
        self.reset_H2()
        pass
    
    def reset_H1(self):        
        self.v['num_nodes']=0
        self.v['with_index']=False
        self.v['directed']=False
        self.v['kinetic']=False

    def reset_H2(self):

        self.l=[]
        self.v['num_fields']=0
        self.v['with_node1']=False
        self.v['with_links']=False
        self.v['with_weight']=False
        self.v['with_coors']=False
        self.v['with_cluster']=False
        self.v['with_color']=False
        self.v['with_size']=False

        self.v['node1']=0
        self.v['node2']=0
        self.v['weight']=0
        self.v['cluster']=0
        self.v['coorx']=0.00
        self.v['coory']=0.00
        self.v['coorz']=0.00
        self.v['color']=''
        self.v['size']=0.00


    def read_H1(self,line):
        line=line.replace('@',''); line=line.replace(',',' ')
        line=line.split()
        for field in line:
            att,opt=field.split('=')
            if att in ['directed','with_index','kinetic']:
                if opt in ['True','true']:
                    self.v[att]=True
                elif opt in ['False','false']:
                    self.v[att]=False
                else:
                    raise Exception('#Error in attribute of network: '+opt)
            elif att in ['num_nodes']:
                self.v[att]=int(opt)
            else:
                raise Exception('# Option in file not recognized: '+att)


    def read_H2(self,line):
        self.reset_H2()
        line=line.replace('#','')
        line=line.split()
        self.v['num_fields']=len(line)
        for ind in range(self.v['num_fields']):
            if line[ind] == 'node':
                if self.v['with_node1']:
                    self.l.append('node2')
                    self.v['node2']=ind
                    self.v['with_links']=True
                else:
                    self.l.append('node1')
                    self.v['node1']=ind
                    self.v['with_node1']=True
            elif line[ind] in self.v.keys():
                self.v[line[ind]]=ind
                self.l.append(line[ind])
                if line[ind] in ['coorx','coory','coorz']: self.v['with_coors']=True
                if line[ind] in ['cluster']: self.v['with_cluster']=True
                if line[ind] in ['color']: self.v['with_color']=True
                if line[ind] in ['size']: self.v['with_size']=True
                if line[ind] in ['weight']: self.v['with_weight']=True

    def read_line(self,line):
        
        line=line.split()
        if len(line)==self.v['num_fields']:
            control=True
            for ind in range(self.v['num_fields']):
                self.v[self.l[ind]]=line[ind]
            if self.v['with_index']: 
                self.v['node1']=int(self.v['node1'])
                self.v['node2']=int(self.v['node2'])
            if self.v['with_weight']:
                self.v['weight']=float(self.v['weight'])
            if self.v['with_cluster']: self.v['cluster']=int(self.v['cluster'])
            if self.v['with_size']: self.v['size']=float(self.v['size'])
            if self.v['with_coors']:
                self.v['coorx']=float(self.v['coorx'])
                self.v['coory']=float(self.v['coory'])
                self.v['coorz']=float(self.v['coorz'])
        else:
            control=False

        return control

class cl_node():
    """ Fundamental unit to constitute a network together with links.

        Attr:
               label   [string] :    label or key
               weight  [int]    :    weight
               link    [dict]   :    

    """
    def __init__(self):
#        """Attributes:
#        
#           label[string]: label or key
#           weight[int]: weight
#           
#        """
        self.label=''
        self.link={}
        self.alt_link={}
        self.k_out=0
        self.k_in=0
        self.k=0
        self.weight=0
        self.cluster=0
        self.component=0
        self.coors=[]
        self.color=''
        self.size=''

    def most_weighted_links(self,length=1):
        """ Indexes **of** the ranked N most weighted links.
        
        :param length: number of links
        :type length: integer
        :return: array of node indexes
        :rtype: list

        .. **note** :: Esto es una nota
        """
        aux_bak=[[self.link[x],x] for x in self.link.keys()]
        aux_bak.sort(reverse=True)
        most_w_destin=[]
        for ii in range(length):
            most_w_destin.append(aux_bak[ii][1])
        return most_w_destin

class cl_cluster():

    def __init__(self):
        """Attributes of a cluster or community
        Args:
            Name blabla

        Returns:
           Bla bla bla
        """
        self.label=''
        self.link={}
        self.alt_link={}
        self.nodes=[]
        self.num_nodes=0
        self.weight=0
        self.k_out=0
        self.k_in=0
        self.k=0

class network():

    def __init_att__(self):

        self.num_nodes=0
        self.num_clusters=0
        self.num_components=0
        self.node=[]
        self.cluster=[]
        self.component=[]
        self.k_total=0
        self.k_max=0
        self.k_out=0
        self.k_in=0
        self.weight=0
        self.labels={}
        self.clustering_method=' '
        self.directed=True
        self.kinetic=False
 
        self.file_net=None
        self.file_labels=None
 
        self.__symmetric__=False

        pass

    def __init_Ts__(self):

        self.Ts=False
        self.T_ind=[]
        self.T_start=[]
        self.T_wl=[]
        self.T_wn=[]        

        pass

    def __init__(self,file_net=None,file_labels=None,net_format='text',labels_format='text',directed=True,kinetic=False,verbose=True):

        self.__init_att__()
        self.__init_Ts__()
        self.directed=directed
        self.kinetic=kinetic

        if file_net!=None:
            self.load_net(file_net,net_format,verbose)
            if file_labels!=None:
                self.load_labels(file_labels,labels_format)
        else:
            if verbose:
                self.info()

        pass

        
    def info(self,update=True,verbose=True):

        if update:

            self.num_nodes=len(self.node)
            self.weight=0
            self.k_total=0
            self.k_max=0
            if self.kinetic:
                for ii in range(self.num_nodes):
                    self.node[ii].weight=sum(self.node[ii].link.values())
            for ii in range(self.num_nodes):
                self.weight+=self.node[ii].weight
                k=len(self.node[ii].link)
                self.node[ii].k_out=k
                self.node[ii].k=k
                self.k_total+=k
                if (self.k_max<k): self.k_max=k

        if verbose:

            print '# Network:'
            print '#', self.num_nodes, 'nodes'
            print '#', self.k_total, 'links out'
            print '#', self.weight, 'total weight nodes'

        pass

    def add_node(self, label=None, weight=0,iout=False):
        """Description add node"""

        node=cl_node()
        node.weight=weight
        self.weight+=weight
        no=self.num_nodes
        self.num_nodes+=1
            
        if label:
            node.label=str(label)
            self.labels[node.label]=no

        self.node.append(node)

        if iout:
            return no
        pass

    def add_link(self,node_origin,node_final,weight=0,index_origin=False,index_final=False,iout=False):

        if index_origin:
            no=node_origin
        else:
            try:
                no=self.labels[str(node_origin)]
            except:
                no=self.add_node(node_origin,iout=True)

        if index_final:
            nf=node_final
        else:
            try:
                nf=self.labels[str(node_final)]
            except:
                nf=self.add_node(node_final,iout=True)

        try:
            self.node[no].link[nf]+=weight
        except:
            self.node[no].link[nf]=weight

        if not self.directed:
            try:
                self.node[nf].link[no]+=weight
            except:
                self.node[nf].link[no]=weight
        
        if iout:
            return no, nf

        pass

    def build_Ts(self,alt_links=False):

        if alt_links:

            alt_k_total=0
            for ii in self.node:
                alt_k_total+=len(ii.alt_link[ii])

            alt_T_ind=numpy.zeros(shape=(alt_k_total),dtype=int,order='Fortran')
            alt_T_start=numpy.zeros(shape=(self.num_nodes+1),dtype=int,order='Fortran')
            alt_T_wl=numpy.zeros(shape=(alt_k_total),dtype=float,order='Fortran')

            kk=0
            for ii in range(self.num_nodes):
                alt_T_start[ii]=kk
                aux_links=self.node[ii].alt_link.keys()
                if ii in aux_links:
                    alt_T_ind[kk]=ii+1
                    alt_T_wl[kk]=self.node[ii].alt_link[ii]
                    aux_links.remove(ii)
                    kk+=1
                for jj in aux_links:
                    alt_T_ind[kk]=jj+1
                    alt_T_wl[kk]=self.node[ii].alt_link[jj]
                    kk+=1
            alt_T_start[self.num_nodes]=kk

            return alt_k_total, alt_T_start, alt_T_ind, alt_T_wl

        else:

            self.info(verbose=False)
            
            self.T_ind=numpy.zeros(shape=(self.k_total),dtype=int,order='Fortran')
            self.T_start=numpy.zeros(shape=(self.num_nodes+1),dtype=int,order='Fortran')
            self.T_wl=numpy.zeros(shape=(self.k_total),dtype=float,order='Fortran')
            self.T_wn=numpy.zeros(shape=(self.num_nodes),dtype=float,order='Fortran')
            
            kk=0
            for ii in range(self.num_nodes):
                self.T_wn[ii]=self.node[ii].weight
                self.T_start[ii]=kk
                aux_links=self.node[ii].link.keys()
                if ii in aux_links:
                    self.T_ind[kk]=ii+1
                    self.T_wl[kk]=self.node[ii].link[ii]
                    aux_links.remove(ii)
                    kk+=1
                for jj in aux_links:
                    self.T_ind[kk]=jj+1
                    self.T_wl[kk]=self.node[ii].link[jj]
                    kk+=1
            self.T_start[self.num_nodes]=kk
            self.Ts=True

        pass

    def build_from_Ts(self):

        if self.Ts:
            self.__init_att__()
            for ii in range(len(self.T_start)-1):
                self.add_node()
                for jj in range(self.T_start[ii],self.T_start[ii+1]):
                    self.node[ii].link[self.T_ind[jj]-1]=self.T_wl[jj]

            self.kinetic=True
            self.info(update=True,verbose=False)

        else:
            print '# Ts arrays needed and self.Ts==True.'
            pass

    def remove_Ts(self):

        self.__init_Ts__()
        pass

    def merge_net(self,net=None,verbose=True):
         
        # merging the labels and weights of nodes
         
        net_to_total=[]
        labels_aux=copy.deepcopy(self.labels)
        for ii in range(net.num_nodes):
            try :
                no=labels_aux[net.node[ii].label]
            except:
                no=self.add_node(net.node[ii].label,iout=True)
            self.node[no].weight+=net.node[ii].weight
            net_to_total.append(no)

        # merging the links
         
        for no in range(net.num_nodes):
            no2=net_to_total[no]
            for nf,wf in net.node[no].link.iteritems():
                nf2=net_to_total[nf]
                self.add_link(no2,nf2,weight=wf,index_origin=True,index_final=True)

        self.info(update=True,verbose=verbose)
        del(net_to_total); del(labels_aux)

        pass

    def extract_net(self,nodes=None,verbose=True):
         
        # extracting the labels and weights of nodes
        
        aux=[-2 for ii in range(self.num_nodes)]
        temp_net=network(directed=self.directed,kinetic=self.kinetic,verbose=False)
        for ii in nodes:
            aux[ii]=temp_net.add_node(self.node[ii].label,self.node[ii].weight,iout=True)
        for ii in nodes:
            aa=aux[ii]
            for jj,kk in self.node[ii].link.iteritems():
                bb=aux[jj]
                if bb>-1:
                    temp_net.add_link(aa,bb,weight=kk,index_origin=True,index_final=True)

        del(aux)
        if temp_net.kinetic:
            for ii in temp_net.node:
                ii.weight=sum(ii.link.values())

        temp_net.info(update=True,verbose=verbose)
        return temp_net
         
    def load_net(self,name_file,format='text',verbose=True):
        """format:['text','native']"""

        self.file_net=name_file

        fff=open(name_file,'r')

        if format=='text':
            io=cl_io()
            for line in fff.readlines():
                if line[0]=='@':
                    io.read_H1(line)
                    if io.v['with_index']:
                        for ii in range(io.v['num_nodes']):
                            self.add_node('')
                        self.labels={}

                if line[0]=='#':
                    io.read_H2(line)
                else:
                    to_read=io.read_line(line)
                    if to_read:
                        if io.v['with_links']:
                            self.add_link(io.v['node1'],io.v['node2'],weight=io.v['weight'],index_origin=io.v['with_index'],index_final=io.v['with_index'])
                            if not io.v['directed']:
                                self.add_link(io.v['node2'],io.v['node1'],weight=io.v['weight'],index_origin=io.v['with_index'],index_final=io.v['with_index'])
                        else:
                            if io.v['with_index']:
                                self.node[io.v['node1']].weight=io.v['weight']
                            else:
                                io.v['node1']=self.add_node(io.v['node1'],weight=io.v['weight'],iout=True)
                            if io.v['with_coors']: self.node[io.v['node1']].coors=[io.v['coorx'],io.v['coory'],io.v['coorz']]
                            if io.v['with_cluster']: self.node[io.v['node1']].cluster=io.v['cluster']
                            if io.v['with_color']: self.node[io.v['node1']].color=io.v['color']
                            if io.v['with_size']: self.node[io.v['node1']].size=io.v['size']

            self.directed=False
            if io.v['directed']:
                self.directed=True

            self.kinetic=False
            if io.v['kinetic']:
                self.kinetic=True
                for ii in self.node:
                    ii.weight=sum(ii.link.values())

            if io.v['with_cluster']:
                jj=0
                # todo

            del(io)

        if format=='native':

            line=fff.readline()
            self.num_nodes=int(line.split()[0])
            
            k_max=int(line.split()[1])
            k_total=int(line.split()[2])
            
            self.T_ind=numpy.zeros(shape=(k_total),dtype=int,order='Fortran')
            self.T_start=numpy.zeros(shape=(self.num_nodes+1),dtype=int,order='Fortran')
            self.T_wl=numpy.zeros(shape=(k_total),dtype=float,order='Fortran')
            
            data=fff.read()

            data2=[]
            for aa in data.split():
                data2.append(float(aa))
                

            jj=-1
            sumk=0
            for ii in range(self.num_nodes):
                jj+=1

                node=cl_node()
                node_ind=int(data2[jj])
                k_out=int(data2[jj+1])
                weight=data2[jj+2]
                self.T_start[ii]=sumk
                
                jj=jj+2
                for kk in range(k_out):
                    jj+=1
                    neigh=int(data2[jj])
                    jj+=1
                    flux=data2[jj]
                    self.T_ind[sumk]=neigh
                    self.T_wl[sumk]=flux
                    sumk+=1
                    node.link[neigh-1]=flux

                node.k_out=len(node.link)
                node.weight=sum(node.link.values())
                self.node.append(node)

            self.T_start[self.num_nodes]=sumk

            self.k_max=k_max
            self.num_links=0
            for ii in self.node:
                self.num_links+=ii.k_out
            self.k_total=k_total
            self.Ts=True
        

        self.weight=0
        for ii in self.node:
            self.weight+=ii.weight

        self.info(verbose=verbose)
        fff.close()
        pass

    def load_labels(self,name_file,format='text'):

        """format=[text,water]"""

        self.file_labels=name_file

        fff=open(name_file,'r')

        if format == 'water':
            for ii in range(self.num_nodes):
                line=ff.readline()
                mss=line.split()[1]+' |'
                for jj in range(2,6):
                    mss=mss+' '+line.split()[jj]
                    mss=mss+' |' 
                for jj in range(6,9):
                    mss=mss+' '+line.split()[jj]
                    mss=mss+' |'
                for jj in range(9,12):
                    mss=mss+' '+line.split()[jj]
                    mss=mss+' |'
                for jj in range(12,15):
                    mss=mss+' '+line.split()[jj]
                    mss=mss+' |'
                for jj in range(15,18):
                    mss=mss+' '+line.split()[jj]
                index=int(line.split()[0])-1
                self.labels[mss]=index
                self.node[index].label=mss

        if format == 'text':
            line=fff.readline()
            line=line.replace('#','')
            line=line.split()
            num_fields=len(line)
            for ind in range(num_fields):
                if line[ind]=='index':
                    nind=ind
                if line[ind]=='label':
                    nlab=ind
            for line in fff.readlines():
                line=line.split()
                ind=int(line.pop(nind))
                lab=str(' ').join(line)
                self.node[ind].label=lab
                self.labels[lab]=ind

        if len(self.labels)>self.num_nodes:
            print '# Some labels have no node'
        if len(self.labels)<self.num_nodes:
            print '# Some nodes have no label'


        fff.close()

    def write_labels (self,name_file=None,format='text'):

        if name_file==None:
            print '# Error: name_file required'
            pass

        fff=open(name_file,'w')
        print >> fff,'#','index','label'
        for ii in range(self.num_nodes):
            print >> fff, ii,self.node[ii].label

        fff.close()
        

    def write_net (self,name_file=None,format='text',pymol=False,with_index=True,with_cluster=False):

        if name_file==None:
            print '# Error: name_file required'
            pass

        #todo: check if the file exists
        
        if format=='native':

            fff=open(name_file,'w')

            print >> fff, self.num_nodes, self.k_max, self.k_total

            for ii in range(self.num_nodes):
                aux=[]
                aux.append(ii+1)
                aux.append(self.node[ii].k_out)
                aux.append(self.node[ii].weight)
                if ii in self.node[ii].link.keys():
                    aux.append(ii+1)
                    aux.append(self.node[ii].link[ii])
                for jj in self.node[ii].link.keys():
                    if jj!=ii:
                        aux.append(jj+1)
                        aux.append(self.node[ii].link[jj])
                aux=str(aux).replace(',','')
                aux=aux.replace('[','')
                aux=aux.replace(']','')
                print >> fff,aux

            fff.close()

        elif format=='text':

            io=cl_io()
            io.v['num_nodes']=self.num_nodes
            io.v['with_index']=with_index
            io.v['directed']=self.directed
            io.v['kinetic']=self.kinetic

            fff=open(name_file,'w')

            header=[]
            if io.v['with_index']:
                header.append('num_nodes='+str(self.num_nodes))
                header.append('with_index='+str(with_index))
            header.append('directed='+str(self.directed))
            if io.v['kinetic']:
                header.append('kinetic='+str(self.kinetic))
            print >> fff, '@ ',', '.join(header)

            io.v['with_weight']=True
            io.v['with_coors']=False
            io.v['with_cluster']=with_cluster
            io.v['with_color']=False
            io.v['with_size']=False

            if pymol:
                io.v['with_coors']=True
            
            line=[]
            line.append('node')
            if io.v['with_weight']: line.append('weight')
            if io.v['with_cluster']: line.append('cluster')
            if io.v['with_size']: line.append('size')
            if io.v['with_color']: line.append('color')
            if io.v['with_coors']: line.append('coorx'); line.append('coory'); line.append('coorz')
            print >> fff,'# ',' '.join(line)

            for ii in range(self.num_nodes):
                line=[]
                node=self.node[ii]
                if io.v['with_index']:
                    line.append(str(ii))
                else:
                    line.append(node.label)
                if io.v['with_weight']: line.append(str(node.weight))
                if io.v['with_cluster']: line.append(str(node.cluster))
                if io.v['with_size']: line.append(str(node.size))
                if io.v['with_color']: line.append(str(node.color))
                if io.v['with_coors']:
                    line.append('  ')
                    aux=[0.0,0.0,0.0]
                    for jj in range(len(node.coors)):
                        aux[jj]=node.coors[jj]
                    for jj in aux: line.append(str(jj))
                print >> fff,' '.join(line)
            
            print >> fff,' '
            print >> fff,'# node node weight'
            
            for ii in range(self.num_nodes):
                for (jj,v) in self.node[ii].link.iteritems():
                    line=[]
                    if io.v['with_index']:
                        line.append(str(ii))
                        line.append(str(jj))
                    else:
                        line.append(self.node[ii].label)
                        line.append(self.node[jj].label)
                    line.append(str(v))
                    print >> fff,' '.join(line)
            
            del(io)
            fff.close()
            




############## FUNCTIONS FOR NETWORKS


    def k_distribution (self, option='out', bins=20, segment=None, delta_x=None,norm=None):

        """ option=['total','in','out']"""

        for ii in self.node:
            ii.k_out=len(ii.link)

        scratch=[]
        if option=='out':
            for ii in self.node:
                scratch.append(ii.k_out)

        xx,yy=pyn_math.histogram(scratch,bins,segment,delta_x,None,norm,plot=False)

        return xx, yy

    def weight_distribution (self, option='all', bins=20, segment=None, delta_x=None,norm=None):

        """option=['all','self_links','all-self_links','links']  (I have to do the same for the links of a node)"""

        scratch=[]
        if option=='all':
            for ii in self.node:
                scratch.append(ii.weight)
                
        if option=='self_links':
            for ii in range(self.num_nodes):
                scratch.append(self.node[ii].link[ii])

        if option=='all-self_links':
            for ii in range(self.num_nodes):
                scratch.append(self.node[ii].weight-self.node[ii].link[ii])

        xx,yy=pyn_math.histogram(scratch,bins,segment,delta_x,None,norm,plot=False)
        
        return xx, yy

    def fpt (self, node_origin, node_sink, num_runs=200, option='mean', bins=20, segment=None, delta_x=None,norm=None):

        """option=['mean','distribution','both','raw']"""

        if self.Ts==False :

            self.build_Ts()

        scratch=[]
        for ii in range(num_runs):
            iseed=numpy.random.random_integers(0,4094,4)
            iseed[3]=(iseed[3]/2)*2+1
            aa=f_net.brownian_run_fpt(self.T_start,self.T_ind,self.T_wl,iseed,node_origin,node_sink,self.num_nodes,self.k_total)
            scratch.append(aa)

        if option in ['distribution','both']:
            xx,yy=pyn_math.histogram(scratch,bins,segment,delta_x,None,norm,plot=False)

        if option in ['mean','both']:
            yy_av,yy_sigma=pyn_math.average(scratch)

        if option=='distribution': return xx,yy
        if option=='mean': return yy_av, yy_sigma
        if option=='both': return yy_av, yy_sigma, xx,yy
        if option=='raw':  return scratch

        pass

    def brownian_walker (self,origin=0,length=None,self_links=True):

        if length==None:
            print '# length needed.'
            pass

        if self.Ts==False :
            self.build_Ts()

        iseed=numpy.random.random_integers(0,4094,4)
        iseed[3]=(iseed[3]/2)*2+1

        if self_links:
            scratch=f_net.brownian_run(1,self.T_start,self.T_ind,self.T_wl,iseed,origin,length,self.num_nodes,self.k_total)
        else:
            scratch=f_net.brownian_run(0,self.T_start,self.T_ind,self.T_wl,iseed,origin,length,self.num_nodes,self.k_total)

        return scratch

    def min_distance (self):

        pass


    def detailed_balance_distance(self,p=1.000):

        if self.Ts==False :

            self.build_Ts()
            
        db_dist=f_net.detailed_balance_distance(p,self.T_start,self.T_ind,self.T_wl,self.num_nodes,self.k_total)

        return db_dist

    def evolution_step(self,vect_in):

        if self.Ts==False:
            self.build_Ts()

        vect_out=f_net.evolution_step(self.T_start,self.T_ind,self.T_wl,vect_in,self.num_nodes,self.k_total)

        return vect_out

    def symmetrize(self,new=True,verbose=True):

        if self.Ts==False :
            self.build_Ts()


        aux_k_total=copy.deepcopy(self.k_total)

        if new:
            temp             = network(verbose=False)
            temp.labels      = copy.deepcopy(self.labels)
            temp.file_net    = copy.deepcopy(self.file_net)
            temp.file_labels = copy.deepcopy(self.file_labels)
            temp.num_nodes   = copy.deepcopy(self.num_nodes)
            temp.directed    = copy.deepcopy(self.directed)
            temp.kinetic     = copy.deepcopy(self.kinetic)
        else:
            temp             = self


        aux={}
        for ii in range(self.num_nodes):
            for jj in self.node[ii].link.keys():
                aux[(ii,jj)]=0
                aux[(jj,ii)]=0
        temp.k_total=len(aux)
        del(aux)        


        pfff=f_net.symmetrize_net(temp.k_total,self.T_ind,self.T_wl,self.T_start,self.num_nodes,aux_k_total)

        temp.k_max=pfff[0]
        temp.T_wl=pfff[1]
        temp.T_ind=pfff[2]
        temp.T_start=pfff[3]
        temp.Ts=True
        temp.weight=0
        temp.node=[]
        for ii in range(temp.num_nodes):
            node=cl_node()
            node.weight=pfff[4][ii]
            for jj in range(temp.T_start[ii],temp.T_start[ii+1]):
                neigh=temp.T_ind[jj]
                node.link[neigh-1]=temp.T_wl[jj]
            node.k_out=temp.T_start[ii+1]-temp.T_start[ii]
            temp.weight+=node.weight
            temp.node.append(node)

        for kk,vv in temp.labels.iteritems():
            temp.node[vv].label=kk

        temp.__symmetric__=True

        if verbose==True :
            temp.info()

        del(pfff)

        if new:
            return temp
        else:
            pass


    def gradient_clusters(self,verbose=True):

        if self.Ts==False :

            self.build_Ts()


        self.num_clusters,pfff=f_net.grad(self.T_ind,self.T_wl,self.T_start,self.num_nodes,self.k_total)


        Clust={}

        for ii in range(self.num_nodes):
            try:
                Clust[pfff[ii]].append(ii)
            except:
                Clust[pfff[ii]]=[]
                Clust[pfff[ii]].append(ii)


        a=0
        for ii in Clust.keys():
            temp=cl_cluster()
            temp.label=self.node[int(ii)].label
            temp.nodes=Clust[ii]
            temp.num_nodes=len(temp.nodes)
            temp.weight=0
            for jj in temp.nodes:
                self.node[jj].cluster=a
                temp.weight+=self.node[jj].weight
            self.cluster.append(temp)
            a+=1


        # Output: self.clust_info, self.representants, self.node_belongs2, self.cluster_weight, self.num_clusters
        if verbose:
            print '# Number of clusters: ',self.num_clusters

    def clusters_links(self,verbose=True):

        if self.Ts==False:

            self.build_Ts()

        if self.num_clusters < 2:

            print '#Error: Number of clusters lower than 2'
            return

        for ii in self.node:
            c_a=ii.cluster
            for jj,ww in ii.link.items():
                c_b=self.node[jj].cluster
                try:
                    self.cluster[c_a].link[c_b]+=ww
                except:
                    self.cluster[c_a].link[c_b]=ww



    def cfep(self,mode='pfold',A=0,B=0,num_bins=100000,num_iter=200000):

        if self.Ts==False:

            self.build_Ts()

        if mode=='pfold':

            A=A+1
            B=B+1
            cfep_out,key_cfep1,key_cfep2=f_net.cfep_pfold(A,B,self.T_ind,self.T_wl,self.T_start,num_bins,num_iter,self.num_nodes,self.k_total)
            return cfep_out,key_cfep1,key_cfep2

        if mode=='mfpt':

            A=A+1
            cfep_out,key_cfep1,key_cfep2=f_net.cfep_mfpt(A,self.T_ind,self.T_wl,self.T_start,num_bins,num_iter,self.num_nodes,self.k_total)
            return cfep_out,key_cfep1,key_cfep2

    def dijkstra(self,node='all',alt_links=False):

        if self.Ts==False :
            self.build_Ts()

        if node=='all':
            node=-1
            dim_out=self.num_nodes
        else:
            node=node+1
            dim_out=1

        opt_directed=0
        if self.directed:
            opt_directed=1

        pfff=f_net.dijkstra(node,dim_out,opt_directed,self.T_start,self.T_ind,self.T_wl,self.num_nodes,self.k_total)

        return pfff

    def mds(self,dim=3,eigenvs='all',output=False,alt_links=False,distances=None,dijkstra=True,stress=False):

        #eigenvs=self.num_nodes
        if eigenvs in ['all','All']:
            eigenvs=self.num_nodes
        if eigenvs>self.num_nodes:
            print '# Error: eigenvs>num_nodes'
            return 
        if dim>eigenvs:
            print '# Error: dim>eigenvs'
            return

        if distances==None: #Just to fill the variable
            distances=numpy.zeros(shape=(1,1),dtype=float,order='Fortran')
            dim_distances=1
        else:
            dim_distances=len(distances)

        opt=0
        opt_stress=0
        opt_directed=0
        if self.directed:
            opt_directed=1

        if dijkstra:
            opt=1
        if stress:
            opt_stress=1

        if alt_links:
            alt_k_total, alt_T_start, alt_T_ind, alt_T_wl=self.build_Ts(alt_links=True)
            o_coors,o_eigenvals,o_eigenvects,o_stress=f_net.mds(opt_directed,opt,opt_stress,dim,eigenvs,alt_T_start,alt_T_ind,alt_T_wl,distances,self.num_nodes,alt_k_total,dim_distances)
            del(alt_k_total, alt_T_start, alt_T_ind, alt_T_wl)
        else:
            if self.Ts==False :
                self.build_Ts()
            o_coors,o_eigenvals,o_eigenvects,o_stress=f_net.mds(opt_directed,opt,opt_stress,dim,eigenvs,self.T_start,self.T_ind,self.T_wl,distances,self.num_nodes,self.k_total,dim_distances)

        for ii in range(self.num_nodes):
            self.node[ii].coors=o_coors[ii][:]

        if output:
            if stress:
                return o_eigenvals,o_eigenvects,o_stress
            else:
                return o_eigenvals,o_eigenvects
        else:
            pass


    def mcl(self,granularity=1.5,eps=0.005,iterations=0,pruning=True,alt_links=False,verbose=True,ss=True):

        ## I have to reset previous clusters

        if alt_links:

            alt_k_total, alt_T_start, alt_T_ind, alt_T_wl=self.build_Ts(alt_links=True)
            self.num_clusters,pfff=f_net.mcl(granularity,eps,iterations,alt_T_start,alt_T_ind,alt_T_wl,self.num_nodes,alt_k_total)
            del(alt_k_total, alt_T_start, alt_T_ind, alt_T_wl)

        else:
            if pruning:

                if ss:
                    fff=open('.input_mcl','w')
                    for ii in range(self.num_nodes):
                        for jj,kk in self.node[ii].link.iteritems():
                            if not jj<ii:
                                print >> fff, ii,jj,kk
                    fff.close()
                else:
                    fff=open('.input_mcl','w')
                    for ii in range(self.num_nodes):
                        for jj,kk in self.node[ii].link.iteritems():
                            print >> fff, ii,jj,kk
                    fff.close()

                comando='mcl .input_mcl --abc -I '+str(granularity)+' -o .output_mcl > /dev/null 2>&1'
                salida=system(comando)

                if salida!=0:
                    print '#Error'
                    exit()
                
                fff=open('.output_mcl','r')
                Clust={}

                cl=0
                for line in fff:    
                    line=line.split()
                    Clust[cl]=[]
                    for aa in line:
                        Clust[cl].append(int(aa))
                    cl+=1

                self.num_clusters=cl
                fff.close()
                salida=system('rm .output_mcl .input_mcl')


            else:
                if self.Ts==False :
                    self.build_Ts()
                    
                self.num_clusters,pfff=f_net.mcl(granularity,eps,iterations,self.T_start,self.T_ind,self.T_wl,self.num_nodes,self.k_total)

                Clust={}
                 
                for ii in range(self.num_nodes):
                    try:
                        Clust[pfff[ii]].append(ii)
                    except:
                        Clust[pfff[ii]]=[]
                        Clust[pfff[ii]].append(ii)

        aux=numpy.array(Clust.keys(),dtype=int)
        weight_clusts=numpy.zeros((self.num_clusters),dtype=float)
        for ii in range(aux.shape[0]):
            for jj in Clust[aux[ii]]:
                weight_clusts[ii]+=self.node[jj].weight

        tosort=weight_clusts.argsort(kind="mergesort")

        self.cluster=[]
        aa=0
        for ii in range(tosort.shape[0]-1,-1,-1):
            kk=tosort[ii]
            jj=aux[kk]
            temp=cl_cluster()
            temp.label=self.node[jj].label
            temp.nodes=Clust[jj]
            temp.num_nodes=len(temp.nodes)
            temp.weight=weight_clusts[kk]
            for ll in temp.nodes:
                self.node[ll].cluster=aa
            self.cluster.append(temp)
            aa+=1

        del(aux,weight_clusts,tosort,aa)
        #for ii in Clust.keys():
        #    temp=cl_cluster()
        #    temp.label=self.node[int(ii)].label
        #    temp.nodes=Clust[ii]
        #    temp.num_nodes=len(temp.nodes)
        #    temp.weight=0
        #    for jj in temp.nodes:
        #        self.node[jj].cluster=a
        #        temp.weight+=self.node[jj].weight
        #    self.cluster.append(temp)
        #    a+=1


        # Output: self.clust_info, self.representants, self.node_belongs2, self.cluster_weight, self.num_clusters
        if verbose:
            print '# Number of clusters: ',self.num_clusters

        pass

    def components(self,alt_links=False,verbose=True):

        if alt_links:
            alt_k_total, alt_T_start, alt_T_ind, alt_T_wl=self.build_Ts(alt_links=True)
            self.num_components,pfff=f_net.components(alt_T_start,alt_T_ind,alt_T_wl,self.num_nodes,alt_k_total)
            del(alt_k_total, alt_T_start, alt_T_ind, alt_T_wl)
        else:
            if self.Ts==False :
                self.build_Ts()
            self.num_components,pfff=f_net.components(self.T_start,self.T_ind,self.T_wl,self.num_nodes,self.k_total)

        Comp={}

        for ii in range(self.num_nodes):
            try:
                Comp[pfff[ii]].append(ii)
            except:
                Comp[pfff[ii]]=[]
                Comp[pfff[ii]].append(ii)


        a=0
        for ii in Comp.keys():
            temp=cl_cluster()
            temp.label=a
            temp.nodes=Comp[ii]
            temp.num_nodes=len(temp.nodes)
            temp.weight=0
            for jj in temp.nodes:
                self.node[jj].component=a
                temp.weight+=self.node[jj].weight
            self.component.append(temp)
            a+=1

        del(Comp)

        if verbose:
            print '# Number of components: ',self.num_components

        pass

#### External Functions

#def traj2net(filename=None,num_particles=0,num_frames=0,output=None):
# 
#    if output==None :
#        ii=filename.rfind('.')
#        if ii>0 :
#            output=filename[:ii]+'.pxn'
#        else:
#            output=filename+'.pxn'
# 
#    if filename.endswith('.bin'):
#        f_net.build_net_bin(filename,output,num_particles,num_frames)
#    else:
#        f_net.build_net(filename,output,num_particles,num_frames)
# 
# 
#    print ' # New network file:', output
#    return None

#class traj2net():
# 
#    def __init__(self,traj=None,num_frames=0,num_parts=0,dimension=0,optimized=False):
# 
#        self.optimized=False
#        self.init=False
#        self.num_frames=0
#        self.num_parts=0
#        self.dimension=0
#        self.keys=[]
#        self.coors=[]
#        self.nodes=[]
# 
#        if traj!=None:
#            self.init=True
#            self.coors=traj
#            self.num_frames=len(traj)
#            self.dimension=len(traj[0])
# 
#        if num_frames!=0 and num_parts!=0 and dimension!=0:
#            self.init=True
#            if optimized:
#                self.optimized=True
# 
#    def append_frame(self,frame=None):
# 
#        #if self.init:
#        #    for ii in range(num_parts):
#        #        for jj in range(dim):
#            
#        pass
# 
#    def append_frame(self):
# 
#        pass


def kinetic_network(traj=None,ranges=None,bins=None,traj_out=False,labels=True,verbose=True):

    prov_net=network(directed=True,kinetic=True,verbose=False)

    ranges=pyn_math.standard_ranges(ranges)
    dimensions=ranges.shape[0]
    traj=pyn_math.standard_traj(traj,dimensions)
    num_frames=traj.shape[0]
    num_parts=traj.shape[1]

    opt_labels=0
    if labels:
        opt_labels=1

    if bins!=None:
        if type(bins) in [int]:
            bins=[bins]
        if len(bins)!=dimensions:
            print '# The length of bins must be equal to the length of ranges'
            return
        bins=numpy.array(bins,dtype=int,order='F')
        traj_net=f_kin_anal.trajbinning2net(opt_labels,traj,ranges,bins,num_frames,num_parts,dimensions)
    else:
        traj_net=f_kin_anal.traj2net(opt_labels,traj,ranges,num_frames,num_parts,dimensions)

    traj_net=pyn_math.standard_traj(traj_net,particles=num_parts,dimensions=1)

    prov_net.Ts=True
    prov_net.T_ind=copy.deepcopy(f_kin_anal.t_ind)
    prov_net.T_wl=copy.deepcopy(f_kin_anal.t_tau)
    prov_net.T_start=copy.deepcopy(f_kin_anal.t_start)
    prov_net.build_from_Ts()
    if opt_labels:
        if bins==None:
            for ii in range(prov_net.num_nodes):
                label=str(f_kin_anal.labels[ii])
                prov_net.node[ii].label=label
                prov_net.labels[label]=ii
        else:
            for ii in range(prov_net.num_nodes):
                label=str(f_kin_anal.labels_daux[ii])
                prov_net.node[ii].label=label
                prov_net.labels[label]=ii


    f_kin_anal.free_memory_ts()
    if verbose:
        prov_net.info(update=False,verbose=True)
        pass
    else:
        pass

    if traj_out:
        return prov_net,traj_net
    else:
        del(traj_net)
        return prov_net


##class kinetic_network(network):
## 
##    def __init__(self,traj=None,ranges=None,verbose=True):
## 
##        self.__init_att__()
##        self.__init_Ts__()
## 
##        try:
##            rango_traj=traj.shape
##        except:
##            traj=array(traj,order='Fortran')
##            rango_traj=traj.shape
## 
##        try:
##            rango_ranges=ranges.shape
##        except:
##            ranges=array(ranges,order='Fortran')
##            rango_ranges=ranges.shape
## 
##        ###
## 
##        if len(rango_ranges)==1:
##            ranges=[ranges]
##            ranges=array(ranges,order='Fortran')
##            rango_ranges=ranges.shape
## 
##        dimensions=rango_ranges[0]
##        if (rango_ranges[1]!=2):
##            print '# Error with ranges'
##            pass
## 
##        if len(rango_traj)==1:
##            if dimensions==1:
##                traj=[expand_dims(traj,1)]
##                traj=array(traj,order='Fortran')
##                rango_traj=traj.shape
##            else:
##                print '# Error with dimension of traj and ranges'
##                pass
## 
##        if len(rango_traj)==2:
##            if dimensions>1:
##                traj=[traj]
##                traj=array(traj,order='Fortran')
##                rango_traj=traj.shape
##            else:
##                traj=expand_dims(traj,2)
##                traj=array(traj,order='Fortran')
##                rango_traj=traj.shape
## 
##        num_parts=rango_traj[0]
##        num_frames=rango_traj[1]
##        if (rango_traj[2]!=dimensions):
##            print '# Error with dimension of traj and ranges'
##            pass
## 
##        len_str=0
##        for aa in ranges:
##            bb=len(str(aa[0]))
##            if (bb>len_str): 
##                len_str=bb 
##            bb=len(str(aa[1]))
##            if (bb>len_str): 
##                len_str=bb 
##                
##        len_str=len_str+1
## 
##        #print num_parts,num_frames,dimensions
##        #print ranges
## 
##        f_kin_anal.traj2net(len_str,traj,ranges,num_parts,num_frames,dimensions)
##        
## 
##        self.Ts=True
##        self.T_ind=copy.deepcopy(f_kin_anal.t_ind)
##        self.T_wl=copy.deepcopy(f_kin_anal.t_tau)
##        self.T_start=copy.deepcopy(f_kin_anal.t_start)
##        self.build_from_Ts()
##        for ii in range(self.num_nodes):
##            label=str(f_kin_anal.labels[ii][:])
##            self.node[ii].label=label
##            self.labels[label]=ii
## 
## 
##        f_kin_anal.free_memory_ts()
##        if verbose:
##            self.info(update=False,verbose=True)
##            pass
##        else:
##            pass
        

    


