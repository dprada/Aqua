from numpy import zeros
from numpy import array
import libcell2box as libcell

class cl_frame():

    def __init__(self):
        self.time=None
        self.step=None
        self.precision=None
        self.lam=None                    # Coming from trr files... what is this?
        self.model=None
        self.coors=[]
        self.box=zeros(shape=(3,3),dtype=float,order='F')
        self.cell=zeros(shape=(3,3),dtype=float,order='F')
        self.orthogonal=0
        self.volume=0.0

    def cell2box(self):
        self.box,self.volume,self.orthogonal=libcell.cell2box(self.cell)

    def box2cell(self):
        self.cell,self.volume,self.orthogonal=libcell.box2cell(self.box)

    def wrap(self):
        libcell.wrap(self.coors,self.box,self.orthogonal,self.coors.shape[0])

