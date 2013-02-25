from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
from pymol import cmd
from numpy import *

def __init__(self):
   self.menuBar.addmenu('Network', 'Network Menu')

   self.menuBar.addmenuitem('Network', 'command',
                      'Add Network',
                      label='Add Network',
                     command = lambda : load_net())

   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add Links',
   #                   label='Add Links',
   #                  command = lambda : load_links())
   # 
   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add Weights',
   #                   label='Add Weights',
   #                  command = lambda : load_weights())
   # 
   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add sets',
   #                   label='Add sets',
   #                  command = lambda : load_sets())
   # 
   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add labels',
   #                   label='Add labels',
   #                  command = lambda : load_labels())

def open_it():
   f_name = tkFileDialog.askopenfilename()
   return f_name

def load_net():
   file_coors= open_it()
   net=network(file_coors)
   for ii in range(net.num_nodes):
      cmd.pseudoatom('net',pos=net.node[ii].coors,resi=ii,chain='A')
      selec='resi '+str(ii)
      cmd.alter(selec,'resn='+'"'+net.node[ii].label+'"')
      cmd.label(selec,'"'+net.node[ii].label+'"')
      cmd.alter(selec,'b='+str(net.node[ii].weight))
      cmd.alter(selec,'q='+str(net.node[ii].cluster))
   for ii in range(net.num_nodes):
      for jj in net.node[ii].link.keys():
         cmd.bond('resi '+str(ii),'resi '+str(jj))
   cmd.spectrum('b','blue_red')
   #cmd.spectrum('q','rainbow')
   #cmd.select('"'+dict1[ii]+'"','resn "'+dict1[ii]+'"')


#############################################################################
#############################################################################
########### From Pynoramix, pyn_cl_net.py

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
                    print '#Error in attribute of network:', opt
            elif att in ['num_nodes']:
                self.v[att]=int(opt)
            else:
                print '# Option in file not recognized:', att

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

    def __init__(self):

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

class cl_cluster():

    def __init__(self):

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
        self.k=0
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
 
        self.Ts=False
        self.T_ind=[]
        self.T_start=[]
        self.T_wl=[]
        self.T_wn=[]        

        pass

    def __init__(self,file_net=None,file_labels=None,directed=True,net_format='text',labels_format='text',verbose=True):

        self.__init_att__()
        self.directed=True

        if file_net!=None:
            self.read_net(file_net,net_format,verbose)
            if file_labels!=None:
                self.read_labels(file_labels,labels_format)
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

    def add_node(self, new_node, weight=0,iout=False):

        node=cl_node()
        node.label=str(new_node)
        node.weight=weight
        self.weight+=weight
        no=self.num_nodes
        self.node.append(node)
        self.labels[node.label]=no
        self.num_nodes+=1
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

    def read_net(self,name_file,format='text',verbose=True):

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
                        else:
                            if io.v['with_index']:
                                self.node[io.v['node1']].weight=io.v['weight']
                            else:
                                io.v['node1']=self.add_node(io.v['node1'],weight=io.v['weight'],iout=True)
                            if io.v['with_coors']: self.node[io.v['node1']].coors=[io.v['coorx'],io.v['coory'],io.v['coorz']]
                            if io.v['with_cluster']: self.node[io.v['node1']].cluster=io.v['cluster']
                            if io.v['with_color']: self.node[io.v['node1']].color=io.v['color']
                            if io.v['with_size']: self.node[io.v['node1']].size=io.v['size']

            self.directed=True
            self.kinetic=False

            if io.v['directed']:
                self.directed=True

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
            
            self.T_ind=zeros(shape=(k_total),dtype=int,order='Fortran')
            self.T_start=zeros(shape=(self.num_nodes+1),dtype=int,order='Fortran')
            self.T_wl=zeros(shape=(k_total),dtype=int,order='Fortran')
            
            data=fff.read()

            data2=[]
            for aa in data.split():
                data2.append(int(aa))
                

            jj=-1
            sumk=0
            for ii in range(self.num_nodes):
                jj+=1

                node=cl_node()
                node_ind=data2[jj]
                k_out=data2[jj+1]
                weight=data2[jj+2]
                self.T_start[ii]=sumk
                
                jj=jj+2
                for kk in range(k_out):
                    jj+=1
                    neigh=data2[jj]
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

    def read_labels(self,name_file,format='text'):

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
                self.node[int(line[nind])].label=line[nlab]
                self.labels[line[nlab]]=int(line[nind])

        if len(self.labels)>self.num_nodes:
            print '# Some labels have no node'
        if len(self.labels)<self.num_nodes:
            print '# Some nodes have no label'


        fff.close()



###def load_sets():
###   dict1={'1':'00-0',
###          '33':'K2-2'}
###   file_sets= open_it()
###   i=0
###   for line in open(file_sets,'r'):
###      i=i+1
###      selec='resi '+str(i)
###      ocup=line.split()
###      cmd.alter(selec,"resn="+'"'+dict1[ocup[0]]+'"')
###      cmd.alter(selec,"q="+ocup[0])
###      
###   for ii in sorted(dict1.iterkeys()):
###      cmd.select('"'+dict1[ii]+'"','resn "'+dict1[ii]+'"')
###   cmd.spectrum("q", 'rainbow')
###   


   










###import tkFileDialog
### 
#### examples:
### 
###def open_it():
### 
###filename = tkFileDialog.askopenfilename()
### 
###print filename # test
### 
### 
### 
###def save_it():
### 
###filename = tkFileDialog.askopenfilename()
### 
###print filename # test
### 
### 
### 
###def save_as():
### 
###filename = tkFileDialog.asksaveasfilename()
### 
###print filename # test
