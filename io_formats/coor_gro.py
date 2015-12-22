from frame import *
from numpy import array
from numpy import dtype

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame

def open_traj_read(file_name):

    io_vars=[0 for ii in range(30)]
    io_pos=0
    io_err=0

    funit=open(file_name,'r')

    if not funit:
        io_err=1
    else:
        io_err=0
        io_pos=funit.tell()
        line=funit.readline()       # Header of the gro file
        natoms=int(funit.readline())    # Number of atoms
        io_vars[0]=natoms
        io_vars[2]=io_pos
        funit.seek(io_pos)

    return funit,io_vars,io_pos,io_err  # io_file,io_vars,io_pos,io_err

def read_all(file_unit,io_vars=None,io_pos=None,wrap=True):

    io_err=0
    io_end=0

    temp=[]
    box_pdb=[0,0,0]
    cell_pdb=[0,0,0]
    read_new=False
    frame=cl_frame()
    while True:
        line=file_unit.readline()                                          # Header of the gro file
        if len(line)>0: 
            line=file_unit.readline()                         
            frame.num_atoms=int(line)                                        # Number of atoms
            for i in range(frame.num_atoms): 
                line=file_unit.readline()
                frame.coors.append(map(float,[line[20:28],line[28:36],line[36:44]]))

            frame.coors=10.0*array(frame.coors,order='F')
            line=file_unit.readline().split()                         # Reading the size of the cell

            frame.box[0,0]=10.0*float(line[0])         
            frame.box[1,1]=10.0*float(line[1])         
            frame.box[2,2]=10.0*float(line[2])       
            if len(line)==9:
                frame.box[0,1]=10.0*float(line[3])
                frame.box[0,2]=10.0*float(line[4])
                frame.box[1,0]=10.0*float(line[5])
                frame.box[1,2]=10.0*float(line[6])
                frame.box[2,0]=10.0*float(line[7])
                frame.box[2,1]=10.0*float(line[8])
            frame.box2cell()
            frame.box2invbox()
            frame.coors=array(frame.coors,dtype=float,order='F')
            if wrap:
                frame.wrap()
            temp.append(frame)
            read_new=False
            frame=cl_frame()
            io_pos=file_unit.tell()
        else:
            del(frame)
            break

    io_end=True

    return temp,io_err,io_end   # io_file,io_err,io_end
        
def read_next(file_unit,io_vars=None,io_pos=None):

    temp_frame,io_pos,io_err,io_end=read_aux(file_unit,io_vars,io_pos)
    return temp_frame,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def read_frame(file_unit,frame,io_vars=None,io_pos=None):

    io_err=1
    io_end=0
    return None,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end


def read_aux(file_unit,io_vars=None,io_pos=None):

    io_end=0
    io_err=0
    temp_frame=cl_frame()

    line=file_unit.readline()                                          # Header of the gro file
    if len(line)==0: 
        io_end=1
    else:
        line=file_unit.readline()                         
        temp_frame.num_atoms=int(line)                                        # Number of atoms

        for i in range(temp_frame.num_atoms): 
            line=file_unit.readline().split()
            temp_frame.coors.append(map(float,line[3:6]))

        temp_frame.coors=10.0*array(temp_frame.coors,order='F')
        line=file_unit.readline().split()                         # Reading the size of the cell

        temp_frame.box[0,0]=10.0*float(line[0])         
        temp_frame.box[1,1]=10.0*float(line[1])         
        temp_frame.box[2,2]=10.0*float(line[2])       
        if len(line)==9:
            temp_frame.box[0,1]=10.0*float(line[3])
            temp_frame.box[0,2]=10.0*float(line[4])
            temp_frame.box[1,0]=10.0*float(line[5])
            temp_frame.box[1,2]=10.0*float(line[6])
            temp_frame.box[2,0]=10.0*float(line[7])
            temp_frame.box[2,1]=10.0*float(line[8])
        temp_frame.box2cell()
        temp_frame.box2invbox()
        temp_frame.coors=array(temp_frame.coors,dtype=float,order='F')
        temp_frame.wrap()
        io_pos=file_unit.tell()

    return temp_frame,io_pos,io_err,io_end   # io_file,io_err,io_end


def close_traj(file_unit):

    io_err=0
    file_unit.close()
    return io_err  #io_err

def write_all():

    pass

def write_frame():

    pass
