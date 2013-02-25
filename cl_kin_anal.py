from pyn_fort_kin_anal import glob as f_kin_anal
from pyn_fort_kin_anal import bindata as libbin
from pyn_cl_net import *
import pyn_math as pyn_math
import numpy
import copy 

## traj: frame,num_part,dim

class pyn_file():

    def __init__ (self,name=None,frames=None,particles=None,dimensions=None,opt='r'):

        self.name=name
        self.opt=opt
        self.file=None
        self.binary=False
        self.unit=None
        self.io=False
        self.ioerr=False
        self.ioend=False
        self.iopos=None
        self.columns=None
        self.frames=frames
        self.particles=particles
        self.dimensions=dimensions
        self.coor_int=False

        if self.name.endswith('.bin'):
            self.binary=True


    def info(self):
        print '# File name:',self.name
        print '# Options:', self.opt
        print '# Binary:', self.binary
        print '# Opened:',self.io
        print '# Position:',self.iopos
        print '# End:',self.ioend
        if self.frames==None:
            print '# Frames: Needed!!!'
        else:
            print '# Frames:',self.frames
        print '# Particles:',self.particles
        print '# Dimensions:',self.dimensions

    def open(self):

        if self.binary:
            self.unit=len(pyn_math.pyn_f90units)+100
            pyn_math.pyn_f90units.append(self.unit)
            libbin.fopen_read(self.unit,self.name)
        else:
            self.file=open(self.name,self.opt)

        self.io=True
        self.ioend=False

    def close(self):

        if self.binary:
            libbin.fclose(self.unit)
            pyn_math.pyn_f90units.remove(self.unit)
            self.unit=None
        else:
            self.file.close()

        self.io=False
        self.ioend=False

    def extract(self):

        if self.binary:
            if self.coor_int:
                dd=numpy.empty((self.frames,self.particles,self.dimensions),dtype=int,order='F')
            else:
                dd=numpy.empty((self.frames,self.particles,self.dimensions),dtype=float,order='F')
            for ii in range(self.frames):
                dd[ii,:,:]=self.read_frame(frame=ii)
            return dd

    def read_frame(self,frame=0):

        if self.coor_int:
            return libbin.read_int_frame(self.unit,frame,self.particles,self.dimensions)
        else:
            return libbin.read_float_frame(self.unit,frame,self.particles,self.dimensions)

    def read_coor(self,frame=0,particle=0,dimension=0):

        if self.coor_int:
            return libbin.read_int_coor(self.unit,frame,particle,dimension,self.particles,self.dimensions)
        else:
            return libbin.read_float_coor(self.unit,frame,particle,dimension,self.particles,self.dimensions)

    def check_length (self):

        if self.coor_int:
            return libbin.check_int_length(self.unit,self.particles,self.dimensions)
        else:
            return libbin.check_float_length(self.unit,self.particles,self.dimensions)


class kinetic_analysis():

    def __init__ (self,traject=None,columns=None,frames=None,particles=None,dimensions=None,\
                      traject_type='Traj',in_ram=True,in_file=False,by_frame=False,verbose=True):

        self.dimensions=dimensions
        self.particles=particles
        self.frames=frames

        self.file_traj=None
        self.file_nodes=None
        self.file_clusters=None

        self.traj=None
        self.traj_nodes=None
        self.traj_clusters=None

        self.network=None
        self.network_clusters=None

        self.__tr_mode_in_ram__   = False
        self.__tr_mode_in_file__  = False
        self.__tr_mode_by_frame__ = False

        self.__no_mode_in_ram__   = False
        self.__no_mode_in_file__  = False
        self.__no_mode_by_frame__ = False

        self.__cl_mode_in_ram__   = False
        self.__cl_mode_in_file__  = False
        self.__cl_mode_by_frame__ = False

        self.__type_nodes__=None
        self.__type_clusters__=None
        self.__offset__=0

        if in_file: in_ram=False

        if traject_type in ['TRAJ','Traj','traj']:
            self.__tr_mode_in_ram__   = in_ram
            self.__tr_mode_in_file__  = in_file
            self.__tr_mode_by_frame__ = by_frame
        elif traject_type in ['NODES','Nodes','nodes']:
            self.__no_mode_in_ram__   = in_ram
            self.__no_mode_in_file__  = in_file
            self.__no_mode_by_frame__ = by_frame
        elif traject_type in ['CLUSTERS','Clusters','clusters']:
            self.__cl_mode_in_ram__   = in_ram
            self.__cl_mode_in_file__  = in_file
            self.__cl_mode_by_frame__ = by_frame

        if traject==None:
            if self.frames!=None:
                if self.dimensions==None:
                    self.dimensions=1
                if self.particles==None:
                    self.particles=1

                if traject_type in ['TRAJ','Traj','traj']:
                    self.traj=numpy.zeros((self.frames,self.particles,self.dimensions),dtype=float,order='F')
                elif traject_type in ['NODES','Nodes','nodes']:
                    self.traj_nodes=numpy.zeros((self.frames,self.particles,self.dimensions),dtype=float,order='F')
                elif traject_type in ['CLUSTERS','Clusters','clusters']:
                    self.traj_clusters=numpy.zeros((self.frames,self.particles,self.dimensions),dtype=float,order='F')

        if type(traject) in [str]:

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

            if self.__tr_mode_in_ram__ and not self.__tr_mode_by_frame__:

                if traject.endswith('.bin'):
                    self.file_traj=pyn_file(name=traject,frames=frames,particles=particles,dimensions=dimensions,opt='r')
                    self.file_traj.open()
                    self.traj=self.file_traj.extract()

                else:

                    self.file_traj=pyn_file(name=traject,opt='r')
                    self.frames = 0
                    
                    self.file_traj.open()
                    for line in self.file_traj.file: self.frames += 1
                    self.file_traj.close()
                    
                    self.traj=numpy.empty((self.frames,self.particles,self.dimensions),dtype=float,order='F')
                    
                    self.file_traj.open()
                    
                    gg=-1
                    for line in self.file_traj.file:
                        line=line.split()
                        gg+=1
                        for kk in range(len(columns)):
                            self.traj[gg,0,kk]=float(line[columns[kk]])
                            
                self.file_traj.close()

                self.traj=pyn_math.standard_traj(self.traj,self.particles,self.dimensions)

            if self.__tr_mode_in_file__:

                self.file_traj=pyn_file(name=traject,frames=frames,particles=particles,dimensions=dimensions)

                if self.dimensions==None or self.particles==None:
                    print '# Error: Input variables "dimensions" and "particles" needed.'
                    return

        if type(traject) in [list,tuple,numpy.ndarray]:
            
            self.traj=pyn_math.standard_traj(traject,particles=self.particles,dimensions=self.dimensions)

            self.frames=self.traj.shape[0]
            self.particles=self.traj.shape[1]
            self.dimensions=self.traj.shape[2]

        if verbose:
            self.info()

    def info(self):

        if self.__tr_mode_in_file__:
            print '# In file', self.file_traj.name,':'

        print '#',self.frames,'frames,',self.particles,'particles,',self.dimensions,'dimensions.'

    def histogram(self,dimension=None,node=None,cluster=None,bins=20,segment=None,delta=None,select_dim=0,norm=False,cumul=False):

        if cluster==None and node==None:
            if self.__tr_mode_in_file__:
                return pyn_math.histogram(self.traj,bins=bins,segment=segment,delta=delta,select_dim=select_dim,norm=norm,cumul=cumul,\
                                          in_file=self.file_traj,by_frame=self.__tr_mode_by_frame__)
            else:
                return pyn_math.histogram(self.traj,bins=bins,segment=segment,delta=delta,select_dim=select_dim,norm=norm,cumul=cumul,\
                                          in_file=False,by_frame=self.__tr_mode_by_frame__)

        if cluster!=None:

            return pyn_math.histogram_mask(self.traj,bins=bins,segment=segment,delta=delta,select_dim=select_dim,\
                                              traj_mask=self.traj_clusters,select_mask=cluster,offset_mask=self.__offset__,\
                                              norm=norm,cumul=cumul)

        if node!=None:

            return pyn_math.histogram_mask(self.traj,bins=bins,segment=segment,delta=delta,select_dim=select_dim,\
                                              traj_mask=self.traj_nodes,select_mask=node,offset_mask=self.__offset__,\
                                              norm=norm,cumul=cumul)


    def flux_cut(self,traj=None,cut=None,verbose=False):

        if cut==None:
            print "# Cut value needed."
            return

        if self.dimensions>1:
            print "# Not implemented yet"
            return

        if traj == None:

            traj_inp=pyn_math.standard_traj(self.traj,particles=self.particles,dimensions=self.dimensions)
            return f_kin_anal.flux_cut(traj_inp,cut,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2])

    def life_time(self,traj=None,state=None,segment=None,sel_dim=0,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_segm=0
        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1
        if (segment):
            opt_segm=1

        if type(sel_dim) in [int,float]:
            sel_dim=[sel_dim]
            num_sel_dim=1
        elif type(sel_dim) in [list,tuple]:
            num_sel_dim=len(sel_dim)

        segments=numpy.zeros((self.dimensions,2),dtype=float,order='F')
        if (opt_segm):
            segment=numpy.array(segment)
            if len(segment.shape)==1:
                segment.resize(1,segment.shape[0])
            for ii in range(num_sel_dim):
                segments[sel_dim[ii],0]=segment[ii,0]
                segments[sel_dim[ii],1]=segment[ii,1]
            num_states=1
            state=[0]
        else:
            if type(state) in [int,float]:
                num_states=1
                state=[state]
            elif type(state) in [list,tuple]:
                num_states=len(state)

        if traj == None:
            if self.__tr_mode_in_file__:
                infile=self.file_traj
                infile.unit=len(pyn_math.pyn_f90units)+1
                pyn_math.pyn_f90units.append(infile.unit)
                lt_mean=f_kin_anal.life_time_dist_infile(infile.name,infile.binary,infile.unit,opt_norm,opt_segm,state,segments,\
                                                             sel_dim,self.particles,self.dimensions,num_states,num_sel_dim)
                pyn_math.pyn_f90units.remove(infile.unit)
                infile.unit=None

            else:
                traj_inp=pyn_math.standard_traj(self.traj,particles=self.particles,dimensions=self.dimensions)
                lt_mean=f_kin_anal.life_time_dist(opt_norm,opt_segm,traj_inp,state,segments,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],num_states)

        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=pyn_math.standard_traj(self.traj_clusters,particles=self.particles,dimensions=1)
            lt_mean=f_kin_anal.life_time_dist(opt_norm,traj_inp,state,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],num_states)
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=pyn_math.standard_traj(self.traj_nodes,particles=self.particles,dimensions=1)
            lt_mean=f_kin_anal.life_time_dist(opt_norm,traj_inp,state,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],num_states)
        else:
            print '# A readable traj is needed'
            return


        lt_dist=copy.deepcopy(f_kin_anal.distrib)
        lt_x=copy.deepcopy(f_kin_anal.distrib_x)
        f_kin_anal.free_distrib()

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


        fpt_mean=f_kin_anal.fpt_dist(opt_norm,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, \
                                                        from_state,from_segment,to_state,to_segment, \
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        from_num_states,to_num_states)

        fpt_dist=copy.deepcopy(f_kin_anal.distrib)
        fpt_x=copy.deepcopy(f_kin_anal.distrib_x)
        f_kin_anal.free_distrib()

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
        fcpt_mean=f_kin_anal.fcpt_dist(opt_norm,opt_noreturn,opt_states,opt_segments, states,segments,commitment_in,\
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        num_states,num_segments,num_commits)

        fcpt_dist=copy.deepcopy(f_kin_anal.distrib)
        fcpt_x=copy.deepcopy(f_kin_anal.distrib_x)
        f_kin_anal.free_distrib()

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

        tt_mean=f_kin_anal.tt_dist(opt_norm,opt_no_return,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, \
                                                        from_state,from_segment,to_state,to_segment, \
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        from_num_states,to_num_states)

        tt_dist=copy.deepcopy(f_kin_anal.distrib)
        tt_x=copy.deepcopy(f_kin_anal.distrib_x)
        f_kin_anal.free_distrib()

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
        ctt_mean=f_kin_anal.ctt_dist(opt_norm,opt_noreturn,opt_states,opt_segments, states,segments,commitment_in,\
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
                                                        num_states,num_segments,num_commits)

        ctt_dist=copy.deepcopy(f_kin_anal.distrib)
        ctt_x=copy.deepcopy(f_kin_anal.distrib_x)
        f_kin_anal.free_distrib()

        if verbose:
            print '# Mean first passage time:', ctt_mean,'frames.'

        if mean:
            return ctt_mean
        else:
            return ctt_x, ctt_dist


    def kinetic_network(self,traj=None,ranges=None,bins=None,labels=True,verbose=False):

        if traj in ['CLUSTERS','Clusters','clusters']:
            if type(self.traj_clusters) not in [numpy.ndarray]:
                self.traj_clusters=numpy.array(self.traj_clusters,order="Fortran")

            self.network_clusters=kinetic_network(self.traj_clusters,ranges=[self.traj_clusters.min(),self.traj_clusters.max()],verbose=verbose)
            return

        elif traj in ['NODES','Nodes','nodes']:
            if type(self.traj_nodes) not in [numpy.ndarray]:
                self.traj_nodes=numpy.array(self.traj_nodes,order="Fortran")
            self.network=kinetic_network(self.traj_nodes,ranges=[self.traj_nodes.min(),self.traj_nodes.max()],verbose=verbose)
            return

        else:
            self.traj=pyn_math.standard_traj(self.traj,particles=self.particles,dimensions=self.dimensions)
            if ranges==None:
                ranges=pyn_math.build_ranges(self.traj)
            else:
                ranges=pyn_math.standard_ranges(ranges)
            self.network,self.traj_nodes=kinetic_network(self.traj,ranges=ranges,bins=bins,traj_out=True,labels=labels,verbose=verbose)                
            return

    def prada1_largo(self,window=None,granularity=1.2,bins=20,ybins=10,segment=None,delta=None,extra_min=False,extra_max=False,\
                         ram=4,increment=1,clusters=True,verbose=False):

        if self.dimensions!=1:
            print '# Method not implemented yet for more than 1D.'
            return

        bins,mmx,mmn,delta=pyn_math.parameters_bins(False,bins,segment,delta)

        rv_min=0
        rv_max=0

        if extra_min:
            rv_min=1
            bins+=1
        if extra_max:
            rv_max=1
            bins+=1

        if verbose:
            if rv_min:
                print '# Extra node for x <', mmn
            if rv_max:
                print '# Extra node for x >', mmx

        if not self.__tr_mode_in_file__:
            '# Error: the trajectory must be in a file'
            return

        self.network=network(kinetic=True,verbose=False)

        self.file_traj.open()

        mmram=ram*1024*1024*1024
        iterations=int(mmram/(self.particles*bins*8))
        
        b_frame=window*increment
        e_frame=self.file_traj.frames-1-window*increment

        first_period=1
        salida=1

        if (b_frame+iterations*increment)>e_frame:
            iterations=((e_frame-b_frame)/increment)

        while (iterations>0):

            print b_frame+iterations*increment
            traj_aux=f_kin_anal.prada1_infile(self.file_traj.unit,ybins,bins,segment[0],delta,rv_min,rv_max,\
                                                  b_frame,iterations,increment,window,self.particles,self.dimensions)
            ranges=pyn_math.build_ranges(traj_aux)
            network_per=kinetic_network(traj_aux,ranges=ranges,traj_out=False,verbose=False)
            self.network.merge_net(network_per,verbose=True)
            del(traj_aux)
            del(network_per)

            b_frame+=iterations*increment
            if (b_frame+iterations*increment)>e_frame:
                iterations=((e_frame-b_frame)/increment)
                
        # cerrar unidad
        self.file_traj.close()

        self.__offset__=window

        if clusters:

            print '# Entra a clusters'
         
            self.network.symmetrize(new=False,verbose=verbose)
         
            self.network.mcl(granularity=granularity,pruning=True,verbose=verbose)
         
            #num_nodes=self.network.num_nodes
            #aux_list=numpy.empty(num_nodes,dtype=int,order='F')
            #for ii in range(num_nodes):
            #    aux_list[ii]=self.network.node[ii].cluster
            # 
            #new_num_frames=self.traj_nodes.shape[0]
            #self.traj_clusters=f_kin_anal.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles)
            # 
            #del(num_nodes,new_num_frames,aux_list)
         
            self.__type_clusters__='prada1'



    def prada1_shell(self,window=None,granularity=1.2,bins=20,ybins=10,segment=None,delta=None,extra_min=False,extra_max=False,\
                         ram=4,increment=1,clusters=True,verbose=False):

        if self.dimensions!=1:
            print '# Method not implemented yet for more than 1D.'
            return

        bins,mmx,mmn,delta=pyn_math.parameters_bins(False,bins,segment,delta)

        rv_min=0
        rv_max=0

        if extra_min:
            rv_min=1
            bins+=1
        if extra_max:
            rv_max=1
            bins+=1

        if verbose:
            if rv_min:
                print '# Extra node for x <', mmn
            if rv_max:
                print '# Extra node for x >', mmx

        if not self.__tr_mode_in_file__:
            '# Error: the trajectory must be in a file'
            return

        self.network=network(kinetic=True,verbose=False)

        self.file_traj.open()

        mmram=ram*1024*1024*1024
        iterations=int(mmram/(self.particles*bins*8))
        
        b_frame=window*increment
        e_frame=self.file_traj.frames-1-window*increment

        first_period=1
        salida=1

        if (b_frame+iterations*increment)>e_frame:
            iterations=((e_frame-b_frame)/increment)

        while (iterations>0):

            print b_frame+iterations*increment
            traj_aux=f_kin_anal.prada1_infile(self.file_traj.unit,ybins,bins,segment[0],delta,rv_min,rv_max,\
                                                  b_frame,iterations,increment,window,self.particles,self.dimensions)
            ranges=pyn_math.build_ranges(traj_aux)
            network_per=kinetic_network(traj_aux,ranges=ranges,traj_out=False,verbose=False)
            self.network.merge_net(network_per,verbose=True)
            del(traj_aux)
            del(network_per)

            b_frame+=iterations*increment
            if (b_frame+iterations*increment)>e_frame:
                iterations=((e_frame-b_frame)/increment)
                
        # cerrar unidad
        self.file_traj.close()

        self.__offset__=window

        if clusters:

            print '# Entra a clusters'
         
            self.network.symmetrize(new=False,verbose=verbose)
         
            self.network.mcl(granularity=granularity,pruning=True,verbose=verbose)
         
            #num_nodes=self.network.num_nodes
            #aux_list=numpy.empty(num_nodes,dtype=int,order='F')
            #for ii in range(num_nodes):
            #    aux_list[ii]=self.network.node[ii].cluster
            # 
            #new_num_frames=self.traj_nodes.shape[0]
            #self.traj_clusters=f_kin_anal.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles)
            # 
            #del(num_nodes,new_num_frames,aux_list)
         
            self.__type_clusters__='prada1'


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

        traj_aux=f_kin_anal.prada1(ybins,bins,mmn,mmx,delta,rv_min,rv_max,self.traj,window,self.particles,self.frames)

        ranges=pyn_math.build_ranges(traj_aux)

        self.network,self.traj_nodes=kinetic_network(traj_aux,ranges=ranges,traj_out=True,verbose=verbose)
         
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
            self.traj_clusters=f_kin_anal.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles)
            
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

        traj_aux=f_kin_anal.prada2(ybins,sbins,bins,mmn,mmx,delta,self.traj,window,self.particles,self.frames)

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
            self.traj_clusters=f_kin_anal.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles,self.dimensions)
         
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

        self.traj_nodes=f_kin_anal.ganna(opt_range,opt,bins,mmn,mmx,delta,rv_min,rv_max,self.traj,ksi,window,self.particles,self.frames)
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
            self.traj_clusters=f_kin_anal.trajnodes2trajclusters(aux_list,self.traj_nodes,num_nodes,new_num_frames,self.particles)

            del(num_nodes,new_num_frames,aux_list)

            self.__type_clusters__='berezovska2012'


    def pca(self,num_eigenvs='ALL',verbose=True):


        if num_eigenvs in ['ALL','All','all']:
            num_eigenvs=self.particles*self.dimensions

        eigenvals,eigenvects=f_kin_anal.pca(num_eigenvs,self.traj,self.frames,self.particles,self.dimensions)

        return eigenvals,eigenvects



            


