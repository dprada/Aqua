from frame import *
from numpy import array

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
    return funit,io_vars,io_pos,io_err # io_file,io_vars,io_pos,io_err

###def read_all(file_unit,io_vars=None,io_pos=None):
### 
###    io_err=0
###    io_end=0
### 
###    model_inds=[]
###    file_unit.seek(0)
###    while True:
###        line=file_unit.readline()
###        if line.startswith('MODEL'):
###            model_inds.append([int(line.split()[1]),file_unit.tell()])
###        if len(line)==0:
###            break
### 
###    #pos=file_unit.tell()
###    #file_unit.seek(pos)
### 
###    if len(model_inds)==0: model_inds.append([1,0])
### 
###    file_unit.seek(0)
###    temp=[]
###    for ref_mod in model_inds:
###        frame=cl_frame()
###        file_unit.seek(ref_mod[1])
###        while True:
###            line=file_unit.readline()
###            ii=line.split()
###            if ii[0]!= 'ATOM':
###                print ii[0]
###            if ii[0]=='CRYST1':
###                print 'dentro'
###                frame.box[0][0]=float(ii[1])
###                frame.box[1][1]=float(ii[2])
###                frame.box[2][2]=float(ii[3])
###                frame.cell=frame.box
###                frame.cell[0][1]=float(ii[4])
###                frame.cell[0][2]=float(ii[5])
###                frame.cell[1][2]=float(ii[6])
###            if (ii[0] in ['ATOM','HETATM']):
###                aux=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
###                frame.coors.append(aux)
###            if ii[0] in ['MODEL','END','ENDMDL']:
###                break
### 
###        frame.coors=array(frame.coors,order='Fortran')
###        temp.append(frame)
### 
###    io_end=True
### 
###    return temp,io_err,io_end # io_file,io_err,io_end

def read_all(file_unit,io_vars=None,io_pos=None):

    io_err=0
    io_end=0

    temp=[]
    box_pdb=[0,0,0]
    cell_pdb=[0,0,0]
    read_new=False
    frame=cl_frame()
    while True:
        line=file_unit.readline()
        if len(line)>0:
            ii=line.split()
            if ii[0]=='CRYST1':
                box_pdb=[float(ii[1]),float(ii[2]),float(ii[3])]
                cell_pdb=[float(ii[4]),float(ii[5]),float(ii[6])]
            if (ii[0] in ['ATOM','HETATM']):
                read_new=True
                aux=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
                frame.coors.append(aux)
            if (ii[0] in ['END','ENDMDL']):
                frame.cell[0,0] = box_pdb [0]
                frame.cell[1,1] = box_pdb [1]
                frame.cell[2,2] = box_pdb [2]
                frame.cell[0,1] = cell_pdb[0]
                frame.cell[0,2] = cell_pdb[1]
                frame.cell[1,2] = cell_pdb[2]
                frame.coors=array(frame.coors,order='Fortran')
                frame.cell2box()
                frame.wrap()
                temp.append(frame)
                read_new=False
                frame=cl_frame()
        else:
            if read_new:
                frame.cell[0,0] = box_pdb [0]
                frame.cell[1,1] = box_pdb [1]
                frame.cell[2,2] = box_pdb [2]
                frame.cell[0,1] = cell_pdb[0]
                frame.cell[0,2] = cell_pdb[1]
                frame.cell[1,2] = cell_pdb[2]
                frame.coors=array(frame.coors,order='Fortran')
                frame.cell2box()
                frame.wrap()
                temp.append(frame)
            else:
                del(frame)
            break

    io_end=True
    return temp,io_err,io_end # io_file,io_err,io_end

def read_next(file_unit,io_vars=None,io_pos=None):

    io_err=1
    io_end=0

    return None,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def read_frame(file_unit,frame,io_vars=None,io_pos=None):

    return None,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def close_traj(file_unit):

    io_err=0
    file_unit.close()
    return io_err  #io_err
    

def write_all():

    pass

def write_frame():

    pass
