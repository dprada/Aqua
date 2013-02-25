from pyn_fort_anal_trajs import glob as f_trajs
from pyn_cl_net import *
import pyn_math as pyn_math
import numpy
import copy 

## traj: frame,num_part,dim

class kinetic_1D_analysis():

    def __init__ (self,traject=None,columns=None,particles=None,dimensions=None,verbose=True):

        self.file_traj=''
        self.file_nodes=''
        self.file_clusters=''
        self.file_traj_columns=None
        self.traj=None
        self.dimensions=dimensions
        self.particles=particles
        self.frames=0
        self.traj_nodes=None
        self.traj_clusters=None
        self.network_nodes=None
        self.network_clusters=None
        self.__mode_in_ram__=True
        self.__mode_in_file__=False
        self.__mode_by_frame__=False
        self.__type_nodes__=None
        self.__type_clusters__=None
        self.__offset__=0

        if traject==None:
            print 'Trajectory needed (name of file or array)'
            return
        
        if type(columns) in [int]:
            self.particles=1
            self.dimensions=1
            columns=[columns]
        elif type(columns) in [tuple,list]:
            if (len(columns))==1:
                self.particles=1
                self.dimensions=1
            else:
                if particles==None and dimensions==None:
                    print '# Please, make use of the variables "particles" or/and "dimensions":'
                    print '#   traj:=[100 frames, 3 dimensions] --> "particles=1" or/and "dimensions=3"'
                    print '#   traj:=[100 frames, 8  particles] --> "particles=8" or/and "dimensions=1"'
                    return
                elif particles==1 or dimensions>=1:
                    self.particles=1
                    self.dimensions=len(columns)
                elif particles>1 or dimensions==1:
                    self.dimensions=1
                    self.particles=len(columns)

        elif columns in ['ALL','All','all']:
            fff=open(traject,'r')
            line=fff.readline().split()
            nn=len(line)
            fff.close()
            columns=[ii for ii in range(nn)]
            if nn>1:
                if particles==None and dimensions==None:
                    print '# Please, make use of the variables "particles" or/and "dimensions":'
                    print '#   traj:=[100 frames, 3 dimensions] --> "particles=1" or/and "dimensions=3"'
                    print '#   traj:=[100 frames, 8  particles] --> "particles=8" or/and "dimensions=1"'
                    return
                elif particles==1 or dimensions>=1:
                    self.particles=1
                    self.dimensions=nn
                elif particles>1 or dimensions==1:
                    self.dimensions=1
                    self.particles=nn
            else:
                self.dimensions=1
                self.particles=1

        if type(traject) in [str]:

            self.file_traj=traject
            self.file_column=columns

            self.frames = 0
            for line in open(traject): self.frames += 1

            self.traj=numpy.empty((self.frames,self.particles,self.dimensions),dtype=float,order='F')

            fff=open(traject,'r')

            gg=-1
            for line in fff:
                line=line.split()
                gg+=1
                for kk in range(len(columns)):
                    self.traj[gg,0,kk]=float(line[columns[kk]])
 
            fff.close()

            self.traj=pyn_math.standard_traj(self.traj,self.particles,self.dimensions)

        if type(traject) in [list,tuple,numpy.ndarray]:
            
            self.traj=pyn_math.standard_traj(traject,particles=self.particles,dimensions=self.dimensions)

        self.frames=self.traj.shape[0]
        self.particles=self.traj.shape[1]
        self.dimensions=self.traj.shape[2]

        if verbose:
            print '# Loaded:'
            self.info()

    def info(self):

        print '# ',self.frames,'frames,',self.particles,'particles,',self.dimensions,'dimensions.'

    def histogram(self,dimension=None,node=None,cluster=None,bins=20,segment=None,delta=None,select_dim=0,norm=False,cumul=False):

        if cluster==None and node==None:
            return pyn_math.histogram(self.traj,bins=bins,segment=segment,delta=delta,select_dim=select_dim,norm=norm,cumul=cumul)

        if cluster!=None:

            xx,yy=pyn_math.histogram_mask(self.traj,bins=bins,segment=segment,delta=delta,select_dim=select_dim,\
                                              traj_mask=self.traj_clusters,select_mask=cluster,offset_mask=self.__offset__,\
                                              norm=norm,cumul=cumul)
            return xx,yy

        if node!=None:

            xx,yy=pyn_math.histogram_mask(self.traj,bins=bins,segment=segment,delta=delta,select_dim=select_dim,\
                                              traj_mask=self.traj_nodes,select_mask=node,offset_mask=self.__offset__,\
                                              norm=norm,cumul=cumul)
            return xx,yy


    def life_time(self,traj=None,state=None,segment=None,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1

        if traj == None:
            traj_inp=self.traj
            aux_dims=self.dimensions
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=self.traj_clusters
            aux_dims=1
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=self.traj_nodes
            aux_dims=1
        else:
            print '# A readable traj is needed'
            return

        traj_inp=pyn_math.standard_traj(traj_inp,particles=self.particles,dimensions=aux_dims)

        if type(state) in [int,float]:
            num_states=1
            state=[state]
        elif type(state) in [list,tuple]:
            num_states=len(state)

        lt_mean=f_trajs.life_time_dist(opt_norm,traj_inp,state,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],num_states)

        lt_dist=copy.deepcopy(f_trajs.distrib)
        lt_x=copy.deepcopy(f_trajs.distrib_x)
        f_trajs.free_distrib()

        if verbose:
            print '# Mean life time:', lt_mean,'frames.'

        if mean:
            return lt_mean
        else:
            return lt_x, lt_dist

    def first_passage_time (self,traj=None,from_state=None,from_segment=None,to_state=None,to_segment=None,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_to_segment=0
        opt_from_segment=0
        opt_to_state=0
        opt_from_state=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1

        if from_state!=None:   
            opt_from_state=1
        else:
            from_state=0

        if from_segment!=None: 
            opt_from_segment=1
        else:
            from_segment=[0.0,0.0]

        if to_state!=None:     
            opt_to_state=1
        else:
            to_state=0

        if to_segment!=None:   
            opt_to_segment=1
        else:
            to_segment=[0.0,0.0]
        
        if opt_to_segment==0 and opt_to_state==0:
            print '# the input variable to_state or to_segment is needed'
            return

        if traj == None:
            traj_inp=self.traj
            aux_dims=self.dimensions
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=self.traj_clusters
            aux_dims=1
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=self.traj_nodes
            aux_dims=1
        else:
            print '# A readable traj is needed'
            return

        traj_inp=pyn_math.standard_traj(traj_inp,particles=self.particles,dimensions=aux_dims)

        if type(from_state) in [int,float]:
            from_num_states=1
            from_state=[from_state]
        elif type(state) in [list,tuple]:
            from_num_states=len(from_state)

        if type(to_state) in [int,float]:
            to_num_states=1
            to_state=[to_state]
        elif type(to_state) in [list,tuple]:
            to_num_states=len(to_state)


        fpt_mean=f_trajs.fpt_dist(opt_norm,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, \
                                                        from_state,from_segment,to_state,to_segment, \
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        from_num_states,to_num_states)

        fpt_dist=copy.deepcopy(f_trajs.distrib)
        fpt_x=copy.deepcopy(f_trajs.distrib_x)
        f_trajs.free_distrib()

        if verbose:
            print '# Mean first passage time:', fpt_mean,'frames.'

        if mean:
            return fpt_mean
        else:
            return fpt_x, fpt_dist


    def first_committed_passage_time (self,traj=None,states=None,segments=None,commitment=None,no_return=False,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_segments=0
        opt_states=0
        opt_noreturn=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1
        if (no_return):
            opt_noreturn=1

        if states!=None:   
            opt_states=1
            segments=[[0,0]]
            num_segments=1
        else:
            opt_segments=0
            states=[0]
            num_states=1

        if opt_segments==0 and opt_states==0:
            print '# the input variable states or segments is needed'
            return

        if traj == None:
            traj_inp=self.traj
            aux_dims=self.dimensions
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=self.traj_clusters
            aux_dims=1
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=self.traj_nodes
            aux_dims=1
        else:
            print '# A readable traj is needed'
            return

        traj_inp=pyn_math.standard_traj(traj_inp,particles=self.particles,dimensions=aux_dims)

        if type(states) in [int,float]:
            num_states=1
            states=[states]
        elif type(states) in [list,tuple]:
            num_states=len(states)

        if opt_segments:
            num_segments=len(segments)

        if commitment==None:
            if opt_segments:
                num_commits=num_segments
            else:
                num_commits=num_states
            commitment=[True for ii in range(num_commits)]
        else:
            num_commits=len(commitment)

        commitment_in=[int(ii) for ii in commitment]
        fcpt_mean=f_trajs.fcpt_dist(opt_norm,opt_noreturn,opt_states,opt_segments, states,segments,commitment_in,\
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        num_states,num_segments,num_commits)

        fcpt_dist=copy.deepcopy(f_trajs.distrib)
        fcpt_x=copy.deepcopy(f_trajs.distrib_x)
        f_trajs.free_distrib()

        if verbose:
            print '# Mean first passage time:', fcpt_mean,'frames.'

        if mean:
            return fcpt_mean
        else:
            return fcpt_x, fcpt_dist


    def trip_time (self,traj=None,from_state=None,from_segment=None,to_state=None,to_segment=None,no_return=False,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_no_return=0
        opt_to_segment=0
        opt_from_segment=0
        opt_to_state=0
        opt_from_state=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1
        if (no_return):
            opt_no_return=1

        if from_state!=None:   
            opt_from_state=1
        else:
            from_state=0

        if from_segment!=None: 
            opt_from_segment=1
        else:
            from_segment=[0.0,0.0]

        if to_state!=None:     
            opt_to_state=1
        else:
            to_state=0

        if to_segment!=None:   
            opt_to_segment=1
        else:
            to_segment=[0.0,0.0]
        
        if opt_to_segment==0 and opt_to_state==0:
            print '# the input variable to_state or to_segment is needed'
            return

        if traj == None:
            traj_inp=self.traj
            aux_dims=self.dimensions
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=self.traj_clusters
            aux_dims=1
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=self.traj_nodes
            aux_dims=1
        else:
            print '# A readable traj is needed'
            return

        traj_inp=pyn_math.standard_traj(traj_inp,particles=self.particles,dimensions=aux_dims)

        if type(from_state) in [int,float]:
            from_num_states=1
            from_state=[from_state]
        elif type(state) in [list,tuple]:
            from_num_states=len(from_state)

        if type(to_state) in [int,float]:
            to_num_states=1
            to_state=[to_state]
        elif type(to_state) in [list,tuple]:
            to_num_states=len(to_state)

        tt_mean=f_trajs.tt_dist(opt_norm,opt_no_return,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, \
                                                        from_state,from_segment,to_state,to_segment, \
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        from_num_states,to_num_states)

        tt_dist=copy.deepcopy(f_trajs.distrib)
        tt_x=copy.deepcopy(f_trajs.distrib_x)
        f_trajs.free_distrib()

        if verbose:
            print '# Mean first passage time:', tt_mean,'frames.'

        if mean:
            return tt_mean
        else:
            return tt_x, tt_dist

    def committed_trip_time (self,traj=None,states=None,segments=None,commitment=None,no_return=False,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_segments=0
        opt_states=0
        opt_noreturn=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1
        if (no_return):
            opt_noreturn=1

        if states!=None:   
            opt_states=1
            segments=[[0,0]]
            num_segments=1
        else:
            opt_segments=0
            states=[0]
            num_states=1

        if opt_segments==0 and opt_states==0:
            print '# the input variable states or segments is needed'
            return

        if traj == None:
            traj_inp=self.traj
            aux_dims=self.dimensions
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=self.traj_clusters
            aux_dims=1
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=self.traj_nodes
            aux_dims=1
        else:
            print '# A readable traj is needed'
            return

        traj_inp=pyn_math.standard_traj(traj_inp,particles=self.particles,dimensions=aux_dims)

        if type(states) in [int,float]:
            num_states=1
            states=[states]
        elif type(states) in [list,tuple]:
            num_states=len(states)

        if opt_segments:
            num_segments=len(segments)

        if commitment==None:
            if opt_segments:
                num_commits=num_segments
            else:
                num_commits=num_states
            commitment=[True for ii in range(num_commits)]
        else:
            num_commits=len(commitment)

        commitment_in=[int(ii) for ii in commitment]
        ctt_mean=f_trajs.ctt_dist(opt_norm,opt_noreturn,opt_states,opt_segments, states,segments,commitment_in,\
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        num_states,num_segments,num_commits)

        ctt_dist=copy.deepcopy(f_trajs.distrib)
        ctt_x=copy.deepcopy(f_trajs.distrib_x)
        f_trajs.free_distrib()

        if verbose:
            print '# Mean first passage time:', ctt_mean,'frames.'

        if mean:
            return ctt_mean
        else:
            return ctt_x, ctt_dist


    def kinetic_network(self,traj=None,ranges=None,verbose=False):

        if traj in ['CLUSTERS','Clusters','clusters']:
            if type(self.traj_clusters) not in [numpy.ndarray]:
                self.traj_clusters=numpy.array(self.traj_clusters,order="Fortran")

            self.network_clusters=kinetic_network(self.traj_clusters,ranges=[self.traj_clusters.min(),self.traj_clusters.max()],verbose=verbose)
            return

        elif traj in ['NODES','Nodes','nodes']:
            if type(self.traj_nodes) not in [numpy.ndarray]:
                self.traj_nodes=numpy.array(self.traj_nodes,order="Fortran")
            self.network_nodes=kinetic_network(self.traj_nodes,ranges=[self.traj_nodes.min(),self.traj_nodes.max()],verbose=verbose)
            return

        else:
            self.traj=pyn_math.standard_traj(self.traj,particles=self.particles,dimensions=self.dimensions)
            if ranges==None:
                ranges=pyn_math.build_ranges(self.traj)
            else:
                ranges=pyn_math.standard_ranges(ranges)
            self.network_nodes,self.traj_nodes=kinetic_network(self.traj,ranges=ranges,traj_out=True,verbose=verbose)                
            return

    def prada1(self,window=None,granularity=1.2,bins=20,ybins=10,segment=None,delta=None,clusters=True,verbose=False):

        if self.dimensions!=1:
            print '# Method not implemented yet for more than 1D.'
            return

        bins,mmx,mmn,delta=pyn_math.parameters_bins(self.traj,bins,segment,delta)

        rv_min=0
        rv_max=0

        if mmn>self.traj.min():
            rv_min=1
            bins+=1
        if mmx<self.traj.max():
            rv_max=1
            bins+=1

        if verbose:
            if rv_min:
                print '# Extra node for x <', mmn
            if rv_max:
                print '# Extra node for x >', mmx

        traj_aux=f_trajs.prada1(ybins,bins,mmn,mmx,delta,rv_min,rv_max,self.traj,window,self.particles,self.frames)

        ranges=pyn_math.build_ranges(traj_aux)

        self.network,self.traj_nodes=kinetic_network(traj_aux,ranges=ranges,traj_out=True,verbose=False)
         
        del(traj_aux)
         
        self.__offset__=window
         
        if clusters:
         
            self.network.symmetrize(new=False,verbose=verbose)
         
            self.network.mcl(granularity=granularity,pruning=True,verbose=verbose)
         
            num_nodes=self.network.num_nodes
            aux_list=numpy.empty(num_nodes,dtype=int,order='F')
            for ii in range(num_nodes):
                aux_list[ii]=self.network.node[ii].cluster
         
            new_num_frames=self.traj_nodes.shape[0]
            self.traj_clusters=f_trajs.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles)
            
            del(num_nodes,new_num_frames,aux_list)
         
            self.__type_clusters__='prada1'


    def prada2(self,window=None,granularity=1.2,bins=10,ybins=10,sbins=10,segment=None,delta=None,clusters=True,verbose=False):

        ref_max=self.traj.max()
        ref_min=self.traj.min()

        if segment==None:
            opt_range=0
            mmx=ref_max
            mmn=ref_min
        else:
            opt_range=1
            mmn=segment[0]
            mmx=segment[1]

        if delta!=None:
            opt=1
        else:
            delta=1.0 # Its given by gannas function
            opt=2

        bins,mmx,mmn,delta=pyn_math.parameters_bins(opt_range,opt,bins,mmn,mmx,delta)

        traj_aux=f_trajs.prada2(ybins,sbins,bins,mmn,mmx,delta,self.traj,window,self.particles,self.frames)

        ranges=pyn_math.build_ranges(traj_aux)
         
        self.network,self.traj_nodes=kinetic_network(traj_aux,ranges=ranges,traj_out=True,verbose=False)
         
        del(traj_aux)
         
        self.__offset__=window
         
        if clusters:
         
            self.network.symmetrize(new=False,verbose=verbose)
         
            self.network.mcl(granularity=granularity,pruning=True,verbose=verbose)
         
            num_nodes=self.network.num_nodes
            aux_list=numpy.empty(num_nodes,dtype=int,order='F')
            for ii in range(num_nodes):
                aux_list[ii]=self.network.node[ii].cluster
         
            new_num_frames=self.traj_nodes.shape[0]
            self.traj_clusters=f_trajs.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles,self.dimensions)
         
            del(num_nodes,new_num_frames,aux_list)
         
            self.__type_clusters__='prada2'



    def berezovska2012(self,window=None,ksi=0.5,granularity=1.2,bins=20,segment=None,delta=None,clusters=True,verbose=False):

        ref_max=self.traj.max()
        ref_min=self.traj.min()
        rv_min=0
        rv_max=0

        if segment==None:
            opt_range=0
            mmx=ref_max
            mmn=ref_min
        else:
            opt_range=1
            mmn=segment[0]
            mmx=segment[1]
            if mmn>ref_min:
                rv_min=1
            if mmx<ref_max:
                rv_max=1

        if delta!=None:
            opt=1
        else:
            delta=1.0 # Its given by gannas function
            opt=2

        if self.dimensions!=1:
            print '# Method not implemented yet for more than 1D.'
            return

        if verbose:
            if rv_min:
                print '# Extra node for x <', mmn
            if rv_max:
                print '# Extra node for x >', mmx

        self.traj_nodes=f_trajs.ganna(opt_range,opt,bins,mmn,mmx,delta,rv_min,rv_max,self.traj,ksi,window,self.particles,self.frames)
        self.__offset__=window

        self.network=kinetic_network(self.traj_nodes,ranges=[self.traj_nodes.min(),self.traj_nodes.max()],verbose=False)
        
        if clusters:

            self.network.symmetrize(new=False,verbose=verbose)

            self.network.mcl(granularity=granularity,pruning=True,verbose=verbose)

            num_nodes=self.network.num_nodes
            aux_list=numpy.empty(num_nodes,dtype=int,order='F')
            for ii in range(num_nodes):
                aux_list[ii]=self.network.node[ii].cluster

            new_num_frames=self.traj_nodes.shape[0]
            self.traj_clusters=f_trajs.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles)

            del(num_nodes,new_num_frames,aux_list)

            self.__type_clusters__='berezovska2012'


    def pca(self,num_eigenvs='ALL',verbose=True):


        if num_eigenvs in ['ALL','All','all']:
            num_eigenvs=self.particles*self.dimensions

        eigenvals,eigenvects=f_trajs.pca(num_eigenvs,self.traj,self.frames,self.particles,self.dimensions)

        return eigenvals,eigenvects



            


