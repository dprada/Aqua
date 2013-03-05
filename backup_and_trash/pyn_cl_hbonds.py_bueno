from numpy import *
#import pyn_water as f_water
from pyn_cl_set import *
#from pyn_cl_net import *
import copy
import pickle as pic
import pyn_hbonds as f_hb

#####################################################################################
##### General H-bonds functions (not only water) 
#####################################################################################

def hbonds_pack (system1=None,select1=None,system2=None,select2=None,verbose=True):
    
    pack_done=True

    if select2!=None and system2==None : 
        system2=system1
        
    f_hb.hbonds.free_memory_in()
    
    if select1==None: select1=system1.list_atoms

    ind_don1=[]
    ind_Hdon1=[]
    num_Hdon1=[]
    ind_acc1=[]

    ind_don2=[]
    ind_Hdon2=[]
    num_Hdon2=[]
    ind_acc2=[]


    if system2!=None:

        if select2==None: select2=system2.list_atoms

        for ii in select1:
            if system1.atom[ii].donor: 
                ind_don1.append(ii+1)
                ind_Hdon1.append(system1.donors_hydrogen[ii])
                num_Hdon1.append(len(system1.donors_hydrogen[ii]))
            if system1.atom[ii].acceptor: ind_acc1.append(ii+1)

        for ii in select2:
            if system2.atom[ii].donor: 
                ind_don2.append(ii+1)
                ind_Hdon2.append(system2.donors_hydrogen[ii])
                num_Hdon2.append(len(system2.donors_hydrogen[ii]))
            if system2.atom[ii].acceptor: ind_acc2.append(ii+1)

        same_set=False

    else:

        for ii in select1:
            if system1.atom[ii].donor and system1.atom[ii].acceptor:
                ind_don1.append(ii+1)
                ind_Hdon1.append(system1.donors_hydrogen[ii])
                num_Hdon1.append(len(system1.donors_hydrogen[ii]))
                ind_acc1.append(ii+1)
            elif system1.atom[ii].donor:
                ind_don2.append(ii+1)
                ind_Hdon2.append(system1.donors_hydrogen[ii])
                num_Hdon2.append(len(system1.donors_hydrogen[ii]))
            elif system1.atom[ii].acceptor:
                ind_acc2.append(ii+1)

        same_set=True


    f_hb.hbonds.ldon1=len(ind_don1)
    f_hb.hbonds.lacc1=len(ind_acc1)
    f_hb.hbonds.ldon2=len(ind_don2)
    f_hb.hbonds.lacc2=len(ind_acc2)
    f_hb.hbonds.maxh1=max(num_Hdon1)
    f_hb.hbonds.maxh2=max(num_Hdon2)

    f_hb.hbonds.initialize_memory_in()

    f_hb.hbonds.ind_don1=array(ind_don1,dtype=int,order='Fortran')
    f_hb.hbonds.num_hdon1=array(num_Hdon1,dtype=int,order='Fortran')
    f_hb.hbonds.ind_acc1=array(ind_acc1,dtype=int,order='Fortran')
    f_hb.hbonds.ind_don2=array(ind_don2,dtype=int,order='Fortran')
    f_hb.hbonds.ind_acc2=array(ind_acc2,dtype=int,order='Fortran')
    f_hb.hbonds.num_hdon2=array(num_Hdon2,dtype=int,order='Fortran')

    for i in range(len(ind_Hdon1)):
        for j in range(len(ind_Hdon1[i])):
            f_hb.hbonds.ind_hdon1[i,j]=ind_Hdon1[i][j]+1

    for i in range(len(ind_Hdon2)):
        for j in range(len(ind_Hdon2[i])):
            f_hb.hbonds.ind_hdon2[i,j]=ind_Hdon2[i][j]+1

    del(ind_don1,ind_Hdon1,num_Hdon1,ind_acc1)
    del(ind_don2,ind_Hdon2,num_Hdon2,ind_acc2)

    if verbose:
        return pack_done,same_set
                    


def hbonds(system1=None,select1=None,system2=None,select2=None,box=None,r_param=3.5,ang_param=30.0,pack=None,verbose=False):
    

    if pack :
        same_set=pack[1]
    else:
        same_set=hbonds_pack (system1=system1,select1=select1,system2=system2,select2=select2)[1]

    if select2!=None and system2==None : 
        system2=system1

    if same_set:

        f_hb.hbonds.natoms1=system1.num_atoms
        f_hb.hbonds.initialize_coors_memory1()
        f_hb.hbonds.lbox1[:,:]=system1.frame[0].box[:,:]
        f_hb.hbonds.coor1[:,:]=system1.frame[0].coors[:,:]  # El paso de las coordenadas del sistema  (esto deberia ser global, tengo que cambiarlo)

        f_hb.hbonds.same_set(r_param,ang_param)
        
        output1=copy.deepcopy(f_hb.hbonds.salida1)
        output2=copy.deepcopy(f_hb.hbonds.salida2)
        output3=copy.deepcopy(f_hb.hbonds.salida3)
        output4=copy.deepcopy(f_hb.hbonds.salida4)
         
        if pack==None:
            f_hb.hbonds.free_memory_in()

        f_hb.hbonds.free_memory_coor()
        f_hb.hbonds.free_memory_out()

        #print 'si'
        #return output1,output2,output3,output4
        return f_hb.hbonds.num_hbs1, f_hb.hbonds.num_hbs2, f_hb.hbonds.num_hbs3, f_hb.hbonds.num_hbs4

    else:

        f_hb.hbonds.natoms1=system1.num_atoms
        f_hb.hbonds.natoms2=system2.num_atoms
        f_hb.hbonds.initialize_coors_memory1()
        f_hb.hbonds.initialize_coors_memory2()


        f_hb.hbonds.lbox1[:,:]=system1.frame[0].box[:,:]
        f_hb.hbonds.coor1[:,:]=system1.frame[0].coors[:,:]  # El paso de las coordenadas del sistema  (esto deberia ser global, tengo que cambiarlo)
        f_hb.hbonds.lbox2[:,:]=system2.frame[0].box[:,:]
        f_hb.hbonds.coor2[:,:]=system2.frame[0].coors[:,:]


        f_hb.hbonds.diff_set(r_param,ang_param)

        output1=copy.deepcopy(f_hb.hbonds.salida1)
        output2=copy.deepcopy(f_hb.hbonds.salida2)

        if pack==None:
            f_hb.hbonds.free_memory_in()

        f_hb.hbonds.free_memory_coor()
        f_hb.hbonds.free_memory_out()

        return f_hb.hbonds.num_hbs1, output1, f_hb.hbonds.num_hbs2, output2

