from frame import *
from numpy import array
import libdcdfile as libdcd

# io_w_vars/io_vars[0]: Number of atoms
# io_w_vars/io_vars[1]: delta_t
# io_w_vars/io_vars[2]: pos_header
# io_w_vars/io_vars[3]: pos_frame
# io_w_vars/io_vars[4]: box initial frame 0
# -----
# io_w_vars/io_vars[10] : Number of frames in the file  (INT)
# io_w_vars/io_vars[11] : Number of previous integration steps  (INT)
# io_w_vars/io_vars[12] : Frequency (integration steps) to save this file  (INT)
# io_w_vars/io_vars[13] : Number of integration steps in the run to create this file  (INT)
# io_w_vars/io_vars[14] : Frequency of coordinate saving  (INT)
# io_w_vars/io_vars[17] : Number of degrees of freedom during the run  (INT)
# io_w_vars/io_vars[18] : Number of fixed atoms  (INT)
# io_w_vars/io_vars[19] : Timestep in AKMA-units. Bit-copy from the 32-bit real number  (INT)
# io_w_vars/io_vars[20] : if crystal lattice information is present in the frames  (INT)
# io_w_vars/io_vars[21] : if this is a 4D trajectory  (INT)
# io_w_vars/io_vars[22] : if fluctuating charges are present  (INT)
# io_w_vars/io_vars[23] : if trajectory is the result of merge without consistency checks  (INT)

# open
# > io_file,io_vars,io_pos,io_err

# close
# read_all
# > [frames],io_err

# read_next
# > temp_frame,io_pos,io_err

# read_frame
# write_all
# write_frame


##

# FFF: unit number for Fortran

def open_traj_read(file_name):

    io_vars=[0 for ii in range(30)]
    io_pos=0
    io_err=0

    funit,o_vars,o_natom,o_delta_t,io_pos=libdcd.open_read(len(file_name),str(file_name))

    if not funit:
        io_err=1
    else:
        io_vars[10:30]=o_vars[0:20]
        io_vars[0]=o_natom
        io_vars[1]=o_delta_t
        io_vars[2]=io_pos
        io_vars[3]=(o_vars[10]*7*8+3*(o_natom*4+8))

    return funit,io_vars,io_pos,io_err  # io_file,io_vars,io_pos,io_err


def read_all(file_unit,io_vars=None,io_pos=None):

    temp=[]
    while 1:
        temp_frame,io_pos,io_err,io_end=read_next(file_unit,io_vars,io_pos)
        if io_end or io_err:
            break
        temp_frame.cell2box()
        temp_frame.wrap()
        temp.append(temp_frame)

    return temp,io_err,io_end   # io_file,io_err,io_end

def read_next(file_unit,io_vars=None,io_pos=None):

    temp_frame=cl_frame()
    #temp_frame.coors=empty((natoms,3),dtype=float32)
    #temp_frame.box=empty((3,3),float32)

    io_pos,temp_frame.cell,temp_frame.coors,io_err,io_end=libdcd.read(file_unit,io_vars[0],io_vars[20],io_pos)
    temp_frame.cell2box()
    temp_frame.wrap()
    return temp_frame,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def read_frame(file_unit,frame,io_vars=None,io_pos=None):

    temp_frame=cl_frame()
    io_pos=io_vars[2]+frame*io_vars[3]
    io_pos,temp_frame.cell,temp_frame.coors,io_err,io_end=libdcd.read(file_unit,io_vars[0],io_vars[20],io_pos)
    temp_frame.cell2box()
    temp_frame.wrap()
    return temp_frame,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def close_traj(file_unit):

    io_err=libdcd.close(file_unit)
    return io_err   #io_err=0 good

def open_traj_write(file_name,io_w_vars,origin_name):

    io_err=0
    funit=libdcd.open_write(len(file_name),str(file_name),io_w_vars[10:30],io_w_vars[0],io_w_vars[1],str(origin_name))

    if not funit:
        io_err=1

    return funit,io_err

def write_frame(file_unit,temp_frame):

    io_err=libdcd.write(file_unit,temp_frame.cell,temp_frame.coors,len(temp_frame.coors))
    return io_err

def close_traj_write(file_unit,io_w_vars):

    io_err=0
    io_err=libdcd.close_write(file_unit,io_w_vars[10:30],io_w_vars[0],io_w_vars[1])
    return io_err   #io_err=0 good

