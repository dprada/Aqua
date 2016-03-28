from frame import *
from numpy import *
from ctypes import *
from numpy.ctypeslib import ndpointer

import os
cwd=os.path.dirname(os.path.abspath(__file__))


# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame

exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, exdrINT, \
exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC, \
exdrNOMEM, exdrENDOFFILE, exdrNR = range(13)


try: 
    xdr=cdll.LoadLibrary(cwd+'/libxdrfile.so')
except:
    raise IOError("libxdrfile.so can't be loaded")


xdr=cdll.LoadLibrary(cwd+'/libxdrfile.so')

def open_traj_read(file_name):

    io_vars=[0 for ii in range(30)]
    io_pos=0
    io_err=0

    funit=xdr.xdrfile_open(file_name,"r")
    if not funit:
        io_err=1
    else:
        natoms=c_int()
        io_err=xdr.read_xtc_natoms(file_name,byref(natoms))
        if io_err:
            return funit,io_vars,io_pos,io_err
        else:
            io_vars[0]=natoms.value
            
        #for NumPy define argtypes - ndpointer is not automatically converted to POINTER(c_float)
        #alternative of ctypes.data_as(POINTER(c_float)) requires two version for numpy and c_float array
            xdr.read_xtc.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float), \
                ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float)]

    return funit,io_vars,io_pos,io_err

def read_all(file_unit,io_vars=None,io_pos=None,wrap=True):

    temp=[]
    while 1:
        temp_frame,io_pos,io_err,io_end=read_next(file_unit,io_vars,io_pos,wrap)
        if io_end or io_err:
            break
        temp.append(temp_frame)

    return temp,io_err,io_end


def read_next (file_unit,io_vars=None,io_pos=None,wrap=True):

    io_err=0
    io_end=0

    step = c_int()
    time = c_float()
    prec = c_float()
    lam = c_float()
    natoms=io_vars[0]

    temp_frame=cl_frame()
    temp_frame.coors=empty((natoms,3),dtype=float32)
    temp_frame.box=empty((3,3),float32)

    result = xdr.read_xtc(file_unit,natoms,byref(step),byref(time),temp_frame.box,
                               temp_frame.coors,byref(prec))

    if result==exdrENDOFFILE:
        io_end=1
    elif result!=exdrOK: 
        io_err=1
    else:
        temp_frame.precision=prec.value
        temp_frame.time=time.value
        temp_frame.step=step.value
        temp_frame.coors=array(10.0*temp_frame.coors,dtype=float,order='F')
        temp_frame.box=array(10.0*temp_frame.box,dtype=float,order='F')
        temp_frame.box2cell()
        temp_frame.box2invbox()
        if wrap: 
            temp_frame.wrap()

    return temp_frame,io_pos,io_err,io_end

#  /* Read one frame of an open xtc file */
#  extern int read_xtc(XDRFILE *xd,int natoms,int *step,float *time,
#                      matrix box,rvec *x,float *prec);
  

        
def read_frame(file_unit,frame,io_vars=None,io_pos=None,wrap=True):

    print 'Not implemented yet.'
    io_err=1
    io_end=0
    return None,io_pos,io_err,io_end

def close_traj(file_unit):

    io_err=0
    xdr.xdrfile_close(file_unit)
    return io_err
    

def write_all():

    pass


# esta funcion a ver...


def open_traj_write(file_name,io_vars):

    io_err=0
    funit=xdr.xdrfile_open(file_name,"w")
    xdr.write_xtc.argtypes=[c_int,POINTER(c_int),POINTER(c_int),POINTER(c_float), \
                           ndpointer(dtype=float32),ndpointer(dtype=float32),POINTER(c_float)]

    if not funit:
        io_err=1

    return funit,io_err

def write_frame(file_unit,temp_frame):

#extern int write_xtc(XDRFILE *xd,
#                       int natoms,int step,float time,
#                       matrix box,rvec *x,float prec);

#    result = xdr.write_xtc(file_unit,natoms,byref(step),byref(time),temp_frame.box,
#                               temp_frame.coors,byref(prec))

    step = c_int()
    time = c_float()
    prec = c_float()
    time.value = 0
    step.value = 1
    prec.value = 0.001

    coors=array(temp_frame.coors/10.0,dtype=float32,order='C')
    box=array(temp_frame.box/10.0,dtype=float32,order='C')

    result = xdr.write_xtc(file_unit,coors.shape[0],step,time,box,coors,prec)




#    pass
# 
# 
# 
#    io_err=0
#    io_end=0
# 
#    step = c_int()
#    time = c_float()
#    prec = c_float()
#    lam = c_float()
#    natoms=io_vars[0]
# 
#    temp_frame=cl_frame()
#    temp_frame.coors=empty((natoms,3),dtype=float32)
#    temp_frame.box=empty((3,3),float32)
# 
#    result = xdr.read_xtc(file_unit,natoms,byref(step),byref(time),temp_frame.box,
#                               temp_frame.coors,byref(prec))
# 
#    if result==exdrENDOFFILE:
#        io_end=1
#    elif result!=exdrOK: 
#        io_err=1
#    else:
#        temp_frame.precision=prec.value
#        temp_frame.time=time.value
#        temp_frame.step=step.value
#        temp_frame.coors=array(10.0*temp_frame.coors,dtype=float,order='F')
#        temp_frame.box=array(10.0*temp_frame.box,dtype=float,order='F')
#        temp_frame.box2cell()
#        temp_frame.box2invbox()
#        if wrap: 
#            temp_frame.wrap()
