import io_formats as io
from os import path

# io_vars[0]: Number of atoms
# io_vars[1]: delta_t
# io_vars[2]: pos_header
# io_vars[3]: pos_frame
# ---
# io_vars[10] : Number of frames in the file  (INT)
# io_vars[11] : Number of previous integration steps  (INT)
# io_vars[12] : Frequency (integration steps) to save this file  (INT)
# io_vars[13] : Number of integration steps in the run to create this file  (INT)
# io_vars[14] : Frequency of coordinate saving  (INT)
# io_vars[17] : Number of degrees of freedom during the run  (INT)
# io_vars[18] : Number of fixed atoms  (INT)
# io_vars[19] : Timestep in AKMA-units. Bit-copy from the 32-bit real number  (INT)
# io_vars[20] : if crystal lattice information is present in the frames  (INT)
# io_vars[21] : if this is a 4D trajectory  (INT)
# io_vars[22] : if fluctuating charges are present  (INT)
# io_vars[23] : if trajectory is the result of merge without consistency checks  (INT)

class traj():
    def __init__(self,file_input=None,frame=None,wrap=True,verbose=False):

        self.io_file=None
        self.io_opened=0
        self.io_end=0
        self.io_err=0
        self.io_pos=None
        self.io_vars=[0 for ii in range(30)]

        self.name=None
        self.type=None
        self.title=None
        self.num_frames=0
        self.total_frames=0
        self.precision=None
        self.frame=[]

        if file_input!=None:
            self.name=file_input
            self.type=file_input.split('.')[-1]
            self.open()

        if frame!=None:
            self.upload_frame(frame,wrap)

        if verbose:
            self.info()
            

    def info(self,index=None):

        if index==None:
            print '#',self.num_frames,'frames/models loaded.'
        else:
            print '#',self.num_frames,'frames/models in traj',index
            

    def open(self):
        
        self.io_file,self.io_vars,self.io_pos,self.io_err=getattr(io,'coor_'+self.type).open_traj_read(self.name)
        if self.io_err: print '# Error opening the file'; return
        self.io_opened=1
        self.total_frames=self.io_vars[10]


    def close(self):

        self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
        if self.io_err: print '# Error closing the file'; return
        self.io_opened=0

    def delete_frame(self,frame='ALL'):
     
        if frame in ['all','All','ALL']:
            del(self.frame)
            self.frame=[]
            self.num_frames=0
            return
     
        elif type(frame) in [list,tuple]:
            for ii in frame:
                self.frame.__delitem__(ii)
                self.num_frames-=1
            return
     
        pass

    def reload_frame(self,frame='next',old=0,wrap=True):

        if frame=='next':
            temp_frame,self.io_pos,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_next(self.io_file,self.io_vars,self.io_pos,wrap)
            if self.io_err: 
                print '# Error reading the file'
                return
            if self.io_end: 
                self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
                if self.io_err: return '# Error closing file'
                self.io_opened=0
                self.io_file=None
                return
            try:
                self.frame[old]=temp_frame
            except:
                self.frame.append(temp_frame)
                self.num_frames+=1
            del(temp_frame)
        else:
            if type(frame) not in [int]:
                print '# Not supported yet.'
                return
            else:
                temp_frame,self.io_pos,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_frame(self.io_file,frame,self.io_vars,self.io_pos,wrap)
                if self.io_err: 
                    return '# Error reading file'
                if self.io_end: 
                    self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
                    if self.io_err: return '# Error closing file'
                    self.io_opened=0
                    self.io_file=None
                    return
                try:
                    self.frame[old]=temp_frame
                except:
                    self.frame.append(temp_frame)
                    self.num_frames+=1
                del(temp_frame)

    def upload_frame(self,frame='next',wrap=True):

        if type(frame) in [int]: frame=[frame]

        ### Uploading next frame
        if frame in ['next','Next','NEXT']:
            temp_frame,self.io_pos,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_next(self.io_file,self.io_vars,self.io_pos,wrap)
            if self.io_err: 
                print '# Error reading the file'
                return
            if self.io_end: 
                self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
                if self.io_err: return '# Error closing file'
                self.io_opened=False
                self.io_file=None
                return
            self.frame.append(temp_frame)
            self.num_frames+=1
        ### Uploading all frames
        elif frame in ['all','All','ALL']:
            temp_frames,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_all(self.io_file,self.io_vars,self.io_pos,wrap)
            if self.io_err: return '# Error reading file'
            self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
            if self.io_err: return '# Error closing file'
            self.io_opened=False
            self.io_file=None
            self.num_frames+=len(temp_frames)
            for ii in temp_frames:
                self.frame.append(ii)
            self.precision=self.frame[0].precision
            del(temp_frames)
        ### Uploading a list of frames
        elif type(frame) in [tuple,list]:
            for ii in frame:
                temp_frame,self.io_pos,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_frame(self.io_file,ii,self.io_vars,self.io_pos,wrap)
                if self.io_err: return '# Error reading file'
                self.frame.append(temp_frame)
                self.num_frames+=1
        ### Error
        else:
            print '# Error in the frame selection'

        pass

    def convert(self,file_name=None,frame='ALL',wrap=True,verbose=False):

        temp_wtraj=wtraj(file_name)
        if temp_wtraj.io_err:
            return

        while 1:
            self.reload_frame(wrap=wrap)
            if self.io_end:
                temp_wtraj.close(verbose)
                break
            temp_wtraj.write(self.frame[0])


class wtraj():
    def __init__(self,file_name=None,verbose=False):

        self.io_file=None
        self.io_opened=0
        self.io_end=0
        self.io_err=0
        self.io_pos=None
        self.io_vars=[0 for ii in range(30)]

        self.name=None
        self.type=None
        self.title=None
        self.precision=None

        if file_name!=None:
            self.open(file_name)

        if verbose:
            self.info()

    def open(self,file_name):

        if path.exists(file_name):
            print '# File',file_name,'already exists.'
            self.io_err=1
            return

        if self.io_opened: print '# This file is already opened.'; self.io_err=1; return
        self.name=file_name
        self.type=file_name.split('.')[-1]
        self.io_vars     = [0 for ii in range(30)]
        self.io_vars[20] = 1     # 1 if crystal lattice information is present in the frames (INT)
        self.io_vars[0]  = self.io_vars[0]  # Num Atoms
        self.io_vars[1]  = 1.000000 # delta_t (by the moment 1.0d0) 
        self.io_vars[12] = 1 # delta steps to save the data
        self.io_vars[29] = 1007 # id number to identify AquaLab traj. 


    def info(self):

        if not self.io_opened:
            print '# No file opened to be written.'
            return
        else:
            print '# File',self.name,'opened to be written.'
            return

    def write(self,frame=None,sel='ALL'):
        if sel in ['ALL','All','all']:
            if not self.io_opened:
                self.io_vars[0]=frame.coors.shape[0]
                self.io_file,self.io_err=getattr(io,'coor_'+self.type).open_traj_write(self.name,self.io_vars)
                if self.io_err: print '# Error writting header of file.'; return
                self.io_opened=1
            self.io_vars[10]+=1
            self.io_err=getattr(io,'coor_'+self.type).write_frame(self.io_file,frame)
        else:
            if type(sel) in [list,tuple]:
                if not self.io_opened:
                    self.io_vars[0]=len(sel)
                    self.io_file,self.io_err=getattr(io,'coor_'+self.type).open_traj_write(self.name,self.io_vars)
                    if self.io_err: print '# Error writting header of file.'; return
                    self.io_opened=1
                self.io_vars[10]+=1
                self.io_err=getattr(io,'coor_'+self.type).write_frame_sel(self.io_file,frame,sel)
            else:
                print '# Input selection as list of atom indexes needed.'
                self.io_err=1

        if self.io_err: print '# Error writting frame in file.'

    def close(self,verbose=False):
        self.io_vars[13]=self.io_vars[10]  # Number of integration steps in the run to create this file  (INT)
        self.io_err=getattr(io,'coor_'+self.type).close_traj_write(self.io_file,self.io_vars)
        if self.io_err: print '# Error closing the file'; return
        self.io_opened=0
        if verbose: '#',self.io_vars[10],'frames written in',self.name

        
