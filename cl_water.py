from numpy import *
import pyn_water as f_water
from pyn_cl_set import *
from pyn_cl_net import *
import copy as ccopy
import pickle as pic

#####################################################################################
##### Water tools
#####################################################################################

# Turning on the pynoramix option in the fortran code:

f_water.main.python=1

#########################

def hbonds_type(option=None,verbose=True):

    hbs_type={}
    hbs_info={}
    hbs_type['Skinner']=1; hbs_info['Skinner']='R.Kumar, J.R. Schmidt and J.L. Skinner. J. Chem. Phys. 126, 204107 (2007)' 
    hbs_type['R(o,h)']=2;  hbs_info['R(o,h)']='V. J. Buch. J. Chem. Phys. 96, 3814-3823 (1992)'
    hbs_type['R(o,o)-Ang(o,o,h)']=3; hbs_info['R(o,o)-Ang(o,o,h)']='A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)'
    hbs_type['Donor-Acceptor-Number']=4; hbs_info['Donor-Acceptor-Number']='A. D. Hammerich, V. J. Buch. J. Chem. Phys. 128, 111101 (2008)'
    hbs_type['Topological']=5; hbs_info['Topological']='R. H. Henchman and S. J. Irudayam. J. Phys. Chem. B. 114, 16792-16810 (2010)'
    hbs_type['Donor-Number-Ang(o,o,h)']=6; hbs_info['Donor-Number-Ang(o,o,h)']='J. D. Smith, C. D. Cappa, et al. Proc. Natl. Acad. Sci. U.S.A. 102, 14171 (2005).'
    hbs_type['Nearest-Neighbour']=7; hbs_info['Nearest-Neighbour']='This is not a hydrogen bond definition but just a topological characterization.'



    if verbose:
        if option not in hbs_type.keys():
            for ii in hbs_type.keys():
                if len(ii)<=12: tab='\t\t\t'
                if 12<len(ii)<=18: tab='\t\t'
                if 18<len(ii): tab='\t'
                print '  ',ii,tab+'[',hbs_info[ii],']'
        return

    if option != None :
        if option not in hbs_type.keys():
            print option, ': Hbond type not defined.'
            print 'List of definitions:'
            for ii in hbs_type.keys():
                if len(ii)<=12: tab='\t\t\t'
                if 12<len(ii)<=18: tab='\t\t'
                if 18<len(ii): tab='\t'
                print '  ',ii,tab+'[',hbs_info[ii],']'
            return 0
        return hbs_type[option]

def hbonds_water(definition=None,system1=None,system2=None,frame=None,sk_param=0.00850,roh_param=2.3000,roo_param=3.5,angooh_param=30.0,optimize=False,verbose=False):

    # Setting up the hbond definition:

    f_water.hbonds.hb_def=hbonds_type(definition,verbose=False)
    if f_water.hbonds.hb_def == 0 : return
    if f_water.hbonds.hb_def == 1 : f_water.hbonds.sk_param=sk_param
    if f_water.hbonds.hb_def == 2 : f_water.hbonds.roh_param= roh_param
    if f_water.hbonds.hb_def == 3 : f_water.hbonds.roo_param, f_water.hbonds.cos_angooh_param= roo_param, cos(radians(angooh_param))
    if f_water.hbonds.hb_def == 4 : pass
    if f_water.hbonds.hb_def == 5 : pass
    if f_water.hbonds.hb_def == 6 : f_water.hbonds.cos_angooh_param= cos(radians(angooh_param))
    if f_water.hbonds.hb_def == 7 : pass

    # Reset of previous hbonds

    for ii in system1.atom :
        ii.hbonds=[]

    for ii in system1.water :
        ii.O.hbonds=[]
        ii.H1.hbonds=[]
        ii.H2.hbonds=[]

    # Frame to be analysed:

    if system2==None and optimize==False:
        if frame==None:
            frame=system1.last_frame

    # Setting up the fortran variables:

    if not optimize :
        f_water.main.list_neighbours=0     # Optimization for hbonds=False 

    f_water.main.nw=system1.num_waters
    f_water.main.initialize_coors_memory()

    __coors2fortran(system1,frame=0)

    f_water.hbonds.initialize_hbonds_memory()

    # Analysis

    f_water.hbonds.hbonds_box()

    # hbonds already in fortran variables... :
    ## f_water.hbonds.num_o2h[i]                 number of hbonds of the atom O of water i
    ## f_water.hbonds.o2h[i,dim=6]               list of water index bonded to water i because of the O
    ## f_water.hbonds.o2which[i,dim=6]           index of hydrogen -1 or 2- corresponding to f_water.hbonds.o2h
    ## f_water.hbonds.strength_o2h[i,dim=6]      skinner parameter corresponding to f_water.hbonds.o2h
    ## f_water.hbonds.num_h2i[i,j-1]             number of hbonds of the atom Hj of water i
    ## f_water.hbonds.h2o[i,j-1,dim=6]           list of water index bonded to water i because of the atom Hj
    ## f_water.hbonds.strength_h2o[i,j-1,dim=6]  skinner parameter corresponding to f_water.hbonds.h2o

    # Reformatting data:
    list_hbonds={}
    for ii in range(system1.num_waters):
        index_water_o=ii
        index_o=system1.water[ii].O.index
        for jj in range(f_water.hbonds.num_o2h[ii]):
            index_water_h=f_water.hbonds.o2h[ii,jj]-1
            if f_water.hbonds.o2which[ii,jj]==1:
                index_h=system1.water[index_water_h].H1.index
                aux_link_h=system1.water[index_water_h].H1.hbonds
            else:
                index_h=system1.water[index_water_h].H2.index
                aux_link_h=system1.water[index_water_h].H2.hbonds
            list_hbonds[str(index_o)+'-'+str(index_h)]=f_water.hbonds.strength_o2h[ii,jj]
            system1.water[index_water_o].O.hbonds.append([index_h,f_water.hbonds.strength_o2h[ii,jj]])
            system1.atom[index_o].hbonds.append([index_h,f_water.hbonds.strength_o2h[ii,jj]])
            aux_link_h.append([index_o,f_water.hbonds.strength_o2h[ii,jj]])
            system1.atom[index_h].hbonds.append([index_o,f_water.hbonds.strength_o2h[ii,jj]])

    # Free memory:
    f_water.main.free_coors_memory()
    f_water.hbonds.free_hbonds_memory()


    if verbose :
        return list_hbonds   # dict: 'index_O'-'indexH'=Skinner_parameter
    else:
        return

def skinner_parameter(system=None,index_wat_o=None,index_wat_h=None,index_h=None,frame=None):


    if frame==None:
        frame=system.last_frame

    f_water.main.nw=system.num_waters
    f_water.main.initialize_coors_memory()

    __coors2fortran(system,frame=0)


    sk=f_water.hbonds.skinner_parameter(index_wat_o,index_wat_h,index_h)

    f_water.main.free_coors_memory()

    return sk 



def mss_water (system=None,output_array=None,definition='Skinner',sk_param=0.00850,roh_param=2.3000,roo_param=3.5,angooh_param=30.0,verbose=True):

    """output_array=['None','microstates','indexes_waters'] """

    if system==None:
        print 'Error: input variables needed'
        print 'mss_water(system=None)'
        return None

    # Setting up the hbond definition:

    f_water.hbonds.hb_def=hbonds_type(definition,verbose=False)
    if f_water.hbonds.hb_def == 0 : return
    if f_water.hbonds.hb_def == 1 : f_water.hbonds.sk_param=sk_param
    if f_water.hbonds.hb_def == 2 : f_water.hbonds.roh_param= roh_param
    if f_water.hbonds.hb_def == 3 : f_water.hbonds.roo_param, f_water.hbonds.cos_angooh_param= roo_param, cos(radians(angooh_param))
    if f_water.hbonds.hb_def == 4 : pass
    if f_water.hbonds.hb_def == 5 : pass
    if f_water.hbonds.hb_def == 6 : f_water.hbonds.cos_angooh_param= cos(radians(angooh_param))
    if f_water.hbonds.hb_def == 7 : pass

    # Initialize Fortran objects:

    f_water.main.nw=system.num_waters
    f_water.main.initialize_coors_memory()

    # Data in Fortran for the frame:

    f_water.main.list_neighbours=0          ## Optimization for hbonds=False in first frame
    __coors2fortran(system,frame=0)

    # Analysis:

    if output_array=='indexes_waters' :
        mss=f_water.microstates.microstates_box_ind_wat(system.num_waters)
    else :
        mss=f_water.microstates.microstates_box(system.num_waters)
        for ii in range(system.num_waters):
            label=str(mss[ii][0])+' |'
            for jj in range(1,5):
                label+=' '+str(mss[ii][jj])
            label+=' |'
            for jj in range(5,8):
                label+=' '+str(mss[ii][jj])
            label+=' |'
            for jj in range(8,11):
                label+=' '+str(mss[ii][jj])
            label+=' |'
            for jj in range(11,14):
                label+=' '+str(mss[ii][jj])
            label+=' |'
            for jj in range(14,17):
                label+=' '+str(mss[ii][jj])
            system.water[ii].microstate=label
        if verbose: print "# Water microstates updated"
        
    # Deallocating Fortran Memory:
    f_water.main.free_coors_memory()

    if output_array in ['indexes_waters','microstates']:
        return mss
    else:
        pass
        

### Private methods for the code:

def coors2fortran(system,frame=None):

    f_water.main.lbox[:,:]=system.frame[frame].box[:,:]

    for jj in range(system.num_waters):
        f_water.main.xarr[jj,0,:]=system.frame[frame].coors[system.water[jj].O.index,:]
        f_water.main.xarr[jj,1,:]=system.frame[frame].coors[system.water[jj].H1.index,:]
        f_water.main.xarr[jj,2,:]=system.frame[frame].coors[system.water[jj].H2.index,:]
        
    pass
    

###class kinetic_network(network):
###    
###    def __init__(self,system=None,file_traj=None,begin=None,end=None,definition='Skinner',sk_param=0.00850,roh_param=2.3000,roo_param=3.5,angooh_param=30.0,verbose=True,optimize=None,memory=1):
###        """ optimize=[None,Disk,RAM], memory=x*Gb"""
### 
###        if optimize in ['RAM','Disk']:
###            l_frames=int((1073741824*float(memory))/(17*4*self.num_waters*1.0)) # 1Gb= 1073741824, 4 each integer
###            
### 
###        if system==None or file_traj==None or begin==None or end==None:
###            print 'Error: input variables needed'
###            print 'kinetic_network(system=None,file_traj=None,begin=None,end=None)'
###            return None
### 
###        # Setting up the hbond definition:
### 
###        f_water.hbonds.hb_def=hbonds_type(definition,verbose=False)
###        if f_water.hbonds.hb_def == 0 : return
###        if f_water.hbonds.hb_def == 1 : f_water.hbonds.sk_param=sk_param
###        if f_water.hbonds.hb_def == 2 : f_water.hbonds.roh_param= roh_param
###        if f_water.hbonds.hb_def == 3 : f_water.hbonds.roo_param, f_water.hbonds.cos_angooh_param= roo_param, cos(radians(angooh_param))
###        if f_water.hbonds.hb_def == 4 : pass
###        if f_water.hbonds.hb_def == 5 : pass
###        if f_water.hbonds.hb_def == 6 : f_water.hbonds.cos_angooh_param= cos(radians(angooh_param))
###        if f_water.hbonds.hb_def == 7 : pass
### 
###        # Frame to be analysed:
### 
###        system.last_frame=begin
### 
###        # Setting up the fortran variables:
### 
###        f_water.main.list_neighbours=0     # Optimization for hbonds=False for the first frame
### 
###        f_water.main.nw=system.num_waters
###        f_water.main.initialize_coors_memory()
### 
### 
###        ####### INITIALIZE NET#####
###        nodes_ant=[0 for ii in range(system.num_waters)]
### 
###        num_nodes=-1
### 
###        self=network(verbose=False)
###        self.file_traj=file_traj
### 
###        ####### INITIALIZE NET#####
### 
###        ###################################### first frame
###        system.delete_coors()
###        system.load_coors(file_traj)
### 
###        f_water.main.lbox[:,:]=system.frame[0].box[:,:]
###        
###        for jj in range(system.num_waters):
###            f_water.main.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
###            f_water.main.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
###            f_water.main.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]
### 
### 
###        mss=f_water.microstates.microstates_box(system.num_waters)
###  
###        ###### NET: 1ST FRAME NODES ########
###        for jj in range(system.num_waters):
###            aa=str(mss[jj])
###            try:
###                nodes_ant[jj]=self.labels[aa]
###            except:
###                num_nodes+=1
###                self.labels[aa]=num_nodes
###                nodes_ant[jj]=num_nodes
###                temp=cl_node()
###                temp.label=aa
###                self.node.append(temp)
### 
###        ###### NET: 1ST FRAME NODES ########
### 
###        system.delete_coors()
###        
###        ###################################### Remaining frames
### 
###        for ii in range(begin+1,end+1):
### 
###            system.load_coors(file_traj)
### 
###            f_water.main.lbox[:,:]=system.frame[0].box[:,:]
###            
###            for jj in range(system.num_waters):
###                f_water.main.xarr[jj,0,:]=system.frame[0].coors[system.water[jj].O.index,:]
###                f_water.main.xarr[jj,1,:]=system.frame[0].coors[system.water[jj].H1.index,:]
###                f_water.main.xarr[jj,2,:]=system.frame[0].coors[system.water[jj].H2.index,:]
### 
###            mss=f_water.microstates.microstates_box(system.num_waters)
### 
###            ###### NET: FRAME NODES ########
###            for jj in range(system.num_waters):
###                bb=str(mss[jj])
###                dd=nodes_ant[jj]
###                try:
###                    cc=self.labels[bb]
###                    try:
###                        self.node[dd].link[cc]+=1
###                    except:
###                        self.node[dd].link[cc]=1
###                    nodes_ant[jj]=cc
###                except:
###                    num_nodes+=1
###                    self.labels[bb]=num_nodes
###                    temp=cl_node()
###                    temp.label=bb
###                    self.node.append(temp)
###                    self.node[dd].link[num_nodes]=1
###                    nodes_ant[jj]=num_nodes
### 
###            ###### NET: FRAME NODES ########
### 
###            system.delete_coors()
###            
### 
###            
###        ################################################# END
### 
###        f_water.main.free_coors_memory()
### 
###        self.num_nodes=len(self.node)
###        self.num_links=0
###        self.weight=0
###        self.k_max=0
###        for ii in range(self.num_nodes):
###            self.node[ii].k_out=len(self.node[ii].link)
###            self.node[ii].weight=sum(self.node[ii].link.values())
###            self.num_links+=self.node[ii].k_out
###            self.weight+=self.node[ii].weight
###            if (self.k_max<self.node[ii].k_out): 
###                self.k_max=self.node[ii].k_out
###        self.k_total=self.num_links
###     
###        self.build_Ts()
###     
###        if verbose:
###            self.info()
### 
###        return 
        

class kinetic_network(network):
    
    def __init__(self,system=None,file_traj=None,begin=None,end=None,definition='Skinner',sk_param=0.00850,roh_param=2.3000,roo_param=3.5,angooh_param=30.0,verbose=True,optimize=None,memory=1):
        """ optimize=[None,Disk,RAM], memory=x*Gb"""

        if optimize in ['RAM','Disk']:
            l_frames=int((1073741824*float(memory))/(17*4*system.num_waters*1.0)) # 1Gb= 1073741824, 4 each integer
            print l_frames
            f_net.funcs.num_frames=l_frames
            f_net.funcs.num_parts=1024
            f_net.funcs.dim_mss=17
            f_net.funcs.init_traj_mss_2_net()
            return


        if system==None or file_traj==None or begin==None or end==None:
            print 'Error: input variables needed'
            print 'kinetic_network(system=None,file_traj=None,begin=None,end=None)'
            return None

        # Setting up the hbond definition:

        f_water.hbonds.hb_def=hbonds_type(definition,verbose=False)
        if f_water.hbonds.hb_def == 0 : return
        if f_water.hbonds.hb_def == 1 : f_water.hbonds.sk_param=sk_param
        if f_water.hbonds.hb_def == 2 : f_water.hbonds.roh_param= roh_param
        if f_water.hbonds.hb_def == 3 : f_water.hbonds.roo_param, f_water.hbonds.cos_angooh_param= roo_param, cos(radians(angooh_param))
        if f_water.hbonds.hb_def == 4 : pass
        if f_water.hbonds.hb_def == 5 : pass
        if f_water.hbonds.hb_def == 6 : f_water.hbonds.cos_angooh_param= cos(radians(angooh_param))
        if f_water.hbonds.hb_def == 7 : pass

        # Frame to be analysed:

        system.last_frame=begin

        # Setting up the fortran variables:

        f_water.main.list_neighbours=0     # Optimization for hbonds=False for the first frame

        f_water.main.nw=system.num_waters
        f_water.main.initialize_coors_memory()


        ####### INITIALIZE NET#####
        nodes_ant=[0 for ii in range(system.num_waters)]

        num_nodes=-1

        self=network(verbose=False)
        self.file_traj=file_traj

        ####### INITIALIZE NET#####

        ###################################### first frame
        system.delete_coors()
        system.load_coors(file_traj)

        coors2fortran(system,frame=0)

        mss=f_water.microstates.microstates_box(system.num_waters)
  
        ###### NET: 1ST FRAME NODES ########
        for jj in range(system.num_waters):
            aa=str(mss[jj])
            try:
                nodes_ant[jj]=self.labels[aa]
            except:
                nodes_ant[jj]=self.add_node(aa,iout=True)

        system.delete_coors()
        
        ###################################### Remaining frames

        for ii in range(begin+1,end+1):

            system.load_coors(file_traj)

            coors2fortran(system,frame=0)

            mss=f_water.microstates.microstates_box(system.num_waters)

            for jj in range(system.num_waters):
                aa,nodes_ant[jj]=self.add_link(nodes_ant[jj],str(mss[jj]),weight=1,index_origin=True,iout=True)

            system.delete_coors()
            

            
        ################################################# END

        f_water.main.free_coors_memory()

        # Computing weight nodes:
        for ii in self.node:
            ii.weight=sum(ii.link.values())


        self.info(update=True,verbose=verbose)

        return 
        




            


