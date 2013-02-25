import numpy
import copy
from pyn_fort_math import glob as fort_math
try:
    import pylab
    wpylab=1
except:
    wpylab=0


pyn_f90units=[100]    # always > 100.

def build_ranges(traj):

    dims=traj.shape[2]
    num_parts=traj.shape[1]

    ranges=numpy.zeros((dims,2),dtype=traj.dtype,order='Fortran')
    for jj in range(dims):
        ranges[jj,0]=traj[:,0,jj].min()
        ranges[jj,1]=traj[:,0,jj].max()
    for ii in range(1,num_parts):
        for jj in range(dims):
            posmin=traj[:,ii,jj].min()
            posmax=traj[:,ii,jj].max()
            if (posmin<ranges[jj,0]): ranges[jj,0]=posmin
            if (posmax>ranges[jj,1]): ranges[jj,1]=posmax

    return ranges

def standard_ranges(ranges):

    if type(ranges) not in [numpy.ndarray]:
        ranges=numpy.array(ranges,order='Fortran')

    if len(ranges.shape)==1:
        ranges.resize(1,ranges.shape[0])

    if ranges.shape[1]!=2:
        print '# Error with ranges'
        pass

    return ranges



def standard_traj(traj,particles=None,dimensions=None):

    if type(traj) not in [numpy.ndarray]:
        traj=numpy.array(traj,dtype=float,order='F')

    if len(traj.shape)==1:
        traj.resize(traj.shape[0],1,1)

    elif len(traj.shape)==2:
        if particles==None and dimensions==None:
            nn=str(traj.shape[-1])
            print '# Error! The trajectory format should be [frames,particles,dimensions].'
            print '# Your input has the shape ['+str(traj.shape[0])+','+str(traj.shape[1])+']:'
            print '# '+nn+' particles with a reaction coordinate? or just a '+nn+'-D reaction coordinate?'
            print '# '
            print '# Please, make use of the variables "particles" or/and "dimensions":'
            print '#   traj:=[100 frames, 3 dimensions] --> "particles=1" or/and "dimensions=3"'
            print '#   traj:=[100 frames, 8  particles] --> "particles=8" or/and "dimensions=1"'
            print '# '
            return 0
        elif particles==1 or dimensions>1:
            traj.resize(traj.shape[0],1,traj.shape[1])
        elif particles>1 or dimensions==1:
            traj.resize(traj.shape[0],traj.shape[1],1)
        elif particles==1 and dimensions==1:
            traj.resize(traj.shape[0],traj.shape[1],1)

    if not numpy.isfortran(traj):
        traj=numpy.array(traj,dtype=float,order='F')

    return traj

#def standard_traj_nodes(traj,particles=None,dimensions=None):
# 
#    if type(traj) not in [numpy.ndarray]:
#        traj=numpy.array(traj,dtype=int,order='F')
# 
#    if len(traj.shape)==1:
#        traj.resize(traj.shape[0],1)
# 
#    if not numpy.isfortran(traj):
#        traj=numpy.array(traj,dtype=int,order='F')
# 
#    return traj


def average(a):

    leng=len(a)

    if leng > 0 :

        v,sigm=fort_math.average(a,leng)

    else :
        
        v,sigm=0.0,0.0

    return v,sigm

def parameters_bins(traj=None,bins=None,segment=None,delta=None):

    if traj==None:
        if segment==None:
            print '# Not implemented yet.'
            return

    if segment==None:
        opt_range=0
        mmx=traj.max()
        mmn=traj.min()
    else:
        opt_range=1
        mmn=segment[0]
        mmx=segment[1]

    if delta!=None:
        opt_delta=1
    else:
        delta=1.0
        opt_delta=2

    bins,mmx,mmn,delta=fort_math.parameters_bins(opt_range,opt_delta,bins,mmn,mmx,delta)

    return bins,mmx,mmn,delta

def histogram(traj,bins=20,segment=None,delta=None,select_dim=0,norm=False,cumul=False,in_file=False,by_frame=False):
    
    infile=in_file
    opt_norm=0
    if norm:
        opt_norm=1
    
    opt_cumul=0
    if cumul:
        opt_cumul=1

    if by_frame:
        print '# Not implemented yet.'
        return

    else:

        if infile:
            infile.unit=len(pyn_f90units)+100
            pyn_f90units.append(infile.unit)
            select_dim+=1
            bins,mmx,mmn,delta=parameters_bins(False,bins,segment,delta)
            fort_math.histogram1d_infile(infile.name,infile.binary,infile.unit,opt_norm,opt_cumul,bins,mmn,mmx,delta,select_dim,\
                                      in_file.frames,in_file.particles,in_file.dimensions)

            h_x=copy.deepcopy(fort_math.histo_x)
            h_y=copy.deepcopy(fort_math.histo_y)
            fort_math.free_mem()

            pyn_f90units.remove(infile.unit)
            infile.unit=None
            return h_x,h_y

        else:
            
            traj=standard_traj(traj)
            
            bins,mmx,mmn,delta=parameters_bins(traj,bins,segment,delta)
            select_dim+=1
            fort_math.histogram1d(opt_norm,opt_cumul,traj,bins,mmn,mmx,delta,select_dim,\
                                      traj.shape[0],traj.shape[1],traj.shape[2])
            
            h_x=copy.deepcopy(fort_math.histo_x)
            h_y=copy.deepcopy(fort_math.histo_y)
            fort_math.free_mem()
            
            return h_x,h_y

    


def histogram_mask(traj,bins=20,segment=None,delta=None,select_dim=0,traj_mask=None,select_mask=None,offset_mask=None,norm=False,cumul=False):
    
    traj=standard_traj(traj)
    traj_mask=standard_traj(traj_mask)

    opt_norm=0
    if norm:
        opt_norm=1
    
    opt_cumul=0
    if cumul:
        opt_cumul=1

    bins,mmx,mmn,delta=parameters_bins(traj,bins,segment,delta)

    if type(select_mask) in [int]:
        select_mask=[select_mask]

    select_mask=numpy.array(select_mask,dtype=int,order='F')

    select_dim+=1
    fort_math.histogram1d_mask(opt_norm,opt_cumul,traj,bins,mmn,mmx,delta,select_dim,\
                                traj_mask,select_mask,offset_mask,traj_mask.shape[0],select_mask.shape[0],\
                                   traj.shape[0],traj.shape[1],traj.shape[2])

    h_x=copy.deepcopy(fort_math.histo_x)
    h_y=copy.deepcopy(fort_math.histo_y)
    fort_math.free_mem()

    return h_x,h_y


def histogram2D(traj,bins=[20,20],segment=None,delta_x=None,prec=None,norm=False,plot=False):

    leng=len(traj)

    if norm==False:
        opt_norm=0
    else:
        opt_norm=1

    if prec==None:
        opt_prec=0
        prec=1.0
    else:
        opt_prec=1

    if segment==None:
        opt_range=0
        mmx0=max(traj[:,0])
        mmx1=max(traj[:,1])
        mmn0=min(traj[:,0])
        mmn1=min(traj[:,1])
    else:
        opt_range=1
        mmn0=segment[0][0]
        mmn1=segment[1][0]
        mmx0=segment[0][1]
        mmx1=segment[1][1]

    if delta_x!=None:
        opt=1
    else:
        delta_x=[1.0,1.0]
        opt=2

    fort_math.histograma_2d(opt_norm,opt_prec,opt_range,opt,traj,bins,[mmn0,mmn1],[mmx0,mmx1],delta_x,prec,leng)

    h_x=copy.deepcopy(fort_math.histo_x)
    h_y=copy.deepcopy(fort_math.histo_y)
    h_z=copy.deepcopy(fort_math.histo_z)
    fort_math.free_mem()
    if plot and wpylab:
        pylab.plot(h_x,h_y,'ro-')

    return h_x,h_y,h_z

def binning(traj=None,bins=20,segment=None,delta_x=None,prec=None):


    if segment==None:
        opt_range=0
        mmx=0.0
        mmn=0.0
    else:
        opt_range=1
        mmn=segment[0]
        mmx=segment[1]

    if delta_x==None:
        opt_delta_x=0
        delta_x=1.0
    else:
        opt_delta_x=1

    if traj==None:

        o_delta_x=fort_math.binning_x(opt_range,opt,bins,mmn,mmx,delta_x)
        h_x=copy.deepcopy(fort_math.histo_x)
        fort_math.free_mem()
        
        return h_x

    else:

        tray_bins=fort_math.binning(opt_range,opt_delta_x,traj,bins,mmn,mmx,delta_x,len(traj))

        h_x=copy.deepcopy(fort_math.histo_x)
    
        fort_math.free_mem()

        return h_x,tray_bins


