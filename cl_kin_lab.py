import numpy


class traj():

    def __init__(self,coors=None,frames=None,particles=None,dimensions=None,dtype=float,verbose=False):

        self.index=None
        self.coors=None
        self.frames=0
        self.particles=0
        self.dimensions=0


        self.__process_coors__(coors,frames,particles,dimensions,dtype)


        if verbose:
            self.info()

    def __process_coors__(self,coors=None,frames=None,particles=None,dimensions=None,dtype=float):

        xx,ff,pp,dd=__process_coors__(coors,self.frames,self.particles,self.dimensions,dtype)
        self.coors=xx, self.frames=ff, self.particles=pp, self.dimensions=dd


    def info(self):
        
        print '#',self.frames,'frames',self.particles,'particles and',self.dimensions,'dimensions.'

    def append(self,coors=None,dtype=float):

        coors=numpy.array(coors,dtype=dtype,order='F')

        if self.frames==0:
            self.coors=coors
            self.frames=1

        else:
            self.coors=numpy.append(self.coors,coors)



class kinetic_lab():

    def __init__(self,num_trajs=None,frames=None,particles=None,dimensions=None,verbose=False):

        self.traj=[]
        self.num_trajs=0
        self.frames=0
        self.particles=0
        self.dimensions=0

        if frames:
            self.frames=frames

        if particles:
            self.particles=particles

        if dimensions:
            self.dimensions=dimensions

        if num_trajs:
            for ii in range(num_trajs):
                self.traj.append(frames=self.frames,particles=self.particles,dimensions=self.dimensions)
        
        if verbose:
            self.info()

    def new_traj(coors=None,frames=None,particles=None,dimensions=None,dtype=float,verbose=False):

        ii=self.num_trajs
        self.traj.append(traj(coors,frames,particles,dimensions,dtype,verbose))
        self.traj[ii].index=ii
        self.num_trajs=ii+1

    def append_coors(coors=None,traj=0,frames=None,particles=None,dimensions=None,dtype=float,verbose=False):

        if not traj<self.num_trajs:
            for ii in range(self.num_trajs,traj+1):
                self.traj.append()
                self.traj[ii].index=ii
            self.num_trajs=traj+1

        self.traj[traj].append(coors,frames,particles,dimensions,dtype,verbose)


    def info(self):

        self.num_trajs=len(self.traj)

        if self.num_trajs==0:
            print '# No trajectories stored.'
        else:
            for ii in self.traj:
                print '# traj',ii.index,'with:',ii.frames,'frames',ii.particles,'particles and',ii.dimensions,'dimensions.'






################################
################################

def __process_coors__(coors=None,frames=None,particles=None,dimensions=None,dtype=float):

    ff=frames
    pp=particles
    dd=dimensions

    if not ff:
        ff=0
    if not pp:
        pp=1
    if not dd:
        dd=1

    if coors==None:

        if ff:
            if pp==1:
                if dd==1:
                    xx=numpy.zeros((ff),dtype=dtype,order='F')
                else:
                    xx=numpy.zeros((ff,dd),dtype=dtype,order='F')
            else:
                if dd==1:
                    xx=numpy.zeros((ff,pp),dtype=dtype,order='F')
                else:
                    xx=numpy.zeros((ff,pp,dd),dtype=dtype,order='F')
            return xx,ff,pp,dd
        else:
            return None,0,0,0

    else:

        xx=numpy.array(coors,dtype=dtype,order='F')
        xxrange=len(xx.shape)

        if xxrange==1:
            ff=xx.shape[0]
            pp=1
            dd=1

        elif xxrange==2:
            if ff==0:
                ff=xx.shape[0]
                if pp==1 and dd==1:
                    dd=xx.shape[1]
                elif pp==1:
                    if dd!=xx.shape[1]:
                        print '# Error: Number of dimensions do not match'
                        return
                elif dd==1:
                    if pp!=xx.shape[1]:
                        print '# Error: Number of particles do not match'
                        return
                else:
                    print '# Error'
                    return

        elif xxrange==3:
            ff=xx.shape[0]
            pp=xx.shape[1]
            dd=xx.shape[2]
                

        return xx,ff,pp,dd

