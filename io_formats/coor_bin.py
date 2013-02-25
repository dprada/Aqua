from frame import *
from numpy import array
import libbinfile as libbin

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame

### This format should be deprected.

#if self.type in ['bin','gro','xtc']:
#    self.coors=10.0*self.coors
#    self.box=10.0*self.box

# io_vars[0]: Number of atoms
# io_vars[1]: box
# io_vars[2]: pos_header
# io_vars[3]: pos_frame

def open_traj_read(file_name):

    io_vars=[0 for ii in range(30)]
    io_pos=0
    io_err=0

    funit,o_natom,o_cell,io_pos=libbin.open_read(len(file_name),str(file_name))
    
    if not funit:
        io_err=1
    else:
        io_vars[0]=o_natom
        io_vars[2]=io_pos
        io_vars[3]=(13+o_natom*3)*4
        io_vars[4]=o_cell

    return funit,io_vars,io_pos,io_err

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
    io_pos,temp_frame.step,temp_frame.time,temp_frame.precision,temp_frame.cell,temp_frame.coors,io_err,io_end=libbin.read(file_unit,io_vars[0],io_pos)
    temp_frame.cell2box()
    temp_frame.wrap()
    return temp_frame,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def read_frame(file_unit,frame,io_vars=None,io_pos=None):

    temp_frame=cl_frame()
    io_pos=io_vars[2]+frame*io_vars[3]
    io_pos,temp_frame.step,temp_frame.time,temp_frame.prec,temp_frame.cell,temp_frame.coors,io_err,io_end=libbin.read(file_unit,io_vars[0],io_pos)
    temp_frame.cell2box()
    temp_frame.wrap()
    return temp_frame,io_pos,io_err,io_end

def close_traj(file_unit):

    io_err=libbin.close(file_unit)
    return io_err

def write_all():
    print '# Not Supported'
    pass

def write_frame():
    print '# Not Supported'
    pass

