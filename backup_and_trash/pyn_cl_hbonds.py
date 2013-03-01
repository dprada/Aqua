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

def hbonds_pack (system1=None,select1=None,system2=None,select2=None,verbose=False):
    
    if select2!=None and system2==None : 
        system2=system1
        
    ind_don1=[]
    ind_Hdon1=[]
    dind_Hdon1={}
    num_Hdon1=[]
    ind_acc1=[]

    ind_don2=[]
    ind_Hdon2=[]
    dind_Hdon2={}
    num_Hdon2=[]
    ind_acc2=[]

    if select1==None: select1=system1.list_atoms

    if system2!=None:

        if select2==None: select2=system2.list_atoms

        for ii in select1:
            if system1.atom[ii].donor: 
                ind_don1.append(ii)
                dind_Hdon1[ii]=system1.donors_hydrogen[ii]
                ind_Hdon1.append(system1.donors_hydrogen[ii])
                num_Hdon1.append(len(system1.donors_hydrogen[ii]))
            if system1.atom[ii].acceptor: ind_acc1.append(ii)

        for ii in select2:
            if system2.atom[ii].donor: 
                ind_don2.append(ii)
                dind_Hdon2[ii]=system2.donors_hydrogen[ii]
                ind_Hdon2.append(system2.donors_hydrogen[ii])
                num_Hdon2.append(len(system2.donors_hydrogen[ii]))
            if system2.atom[ii].acceptor: ind_acc2.append(ii)

        max_Hdon1=max(num_Hdon1)
        max_Hdon2=max(num_Hdon2)
        same_set=False

    else:

        for ii in select1:
            if system1.atom[ii].donor and system1.atom[ii].acceptor:
                ind_don1.append(ii)
                dind_Hdon1[ii]=system1.donors_hydrogen[ii]
                ind_Hdon1.append(system1.donors_hydrogen[ii])
                num_Hdon1.append(len(system1.donors_hydrogen[ii]))
                ind_acc1.append(ii)
            elif system1.atom[ii].donor:
                ind_don2.append(ii)
                dind_Hdon2[ii]=system1.donors_hydrogen[ii]
                ind_Hdon2.append(system1.donors_hydrogen[ii])
                num_Hdon2.append(len(system1.donors_hydrogen[ii]))
            elif system1.atom[ii].acceptor:
                ind_acc2.append(ii)

        max_Hdon1=max(num_Hdon1)
        max_Hdon2=max(num_Hdon2)
        same_set=True

    Ldon1=len(ind_don1)
    Lacc1=len(ind_acc1)
    Ldon2=len(ind_don2)
    Lacc2=len(ind_acc2)

    return Ldon1,Lacc1,Ldon2,Lacc2,ind_don1,ind_Hdon1,num_Hdon1,ind_acc1,ind_don2,ind_Hdon2,num_Hdon2,ind_acc2,dind_Hdon1,dind_Hdon2,max_Hdon1,max_Hdon2,same_set
                    


def hbonds(system1=None,select1=None,system2=None,select2=None,box=None,r_param=3.5,ang_param=30.0,pack=None,verbose=False):
    
    if box==None:
        box=system1.frame[0].box

    if select2!=None and system2==None : 
        system2=system1

    if pack :

        Ldon1,Lacc1,Ldon2,Lacc2=pack[0],pack[1],pack[2],pack[3]
        ind_don1,ind_Hdon1,num_Hdon1,ind_acc1=pack[4],pack[5],pack[6],pack[7]
        ind_don2,ind_Hdon2,num_Hdon2,ind_acc2=pack[8],pack[9],pack[10],pack[11]
        dind_Hdon1,dind_Hdon2=pack[12],pack[13]
        max_Hdon1,max_Hdon2=pack[14],pack[15]
        same_set=pack[16]

    else:
        
        ind_don1=[]
        dind_Hdon1={}
        ind_Hdon1=[]
        num_Hdon1=[]
        ind_acc1=[]

        ind_don2=[]
        dind_Hdon2={}
        ind_Hdon2=[]
        num_Hdon2=[]
        ind_acc2=[]

        if select1==None: select1=system1.list_atoms

        if system2!=None:

            if select2==None: select2=system2.list_atoms

            for ii in select1:
                if system1.atom[ii].donor: 
                    ind_don1.append(ii)
                    dind_Hdon1[ii]=system1.donors_hydrogen[ii]
                    ind_Hdon1.append(system1.donors_hydrogen[ii])
                    num_Hdon1.append(len(system1.donors_hydrogen[ii]))
                if system1.atom[ii].acceptor: ind_acc1.append(ii)

            for ii in select2:
                if system2.atom[ii].donor: 
                    ind_don2.append(ii)
                    dind_Hdon2[ii]=system2.donors_hydrogen[ii]
                    ind_Hdon2.append(system2.donors_hydrogen[ii])
                    num_Hdon2.append(len(system2.donors_hydrogen[ii]))
                if system2.atom[ii].acceptor: ind_acc2.append(ii)

            max_Hdon1=max(num_Hdon1)
            max_Hdon2=max(num_Hdon2)
            same_set=False

        else:

            for ii in select1:
                if system1.atom[ii].donor and system1.atom[ii].acceptor:
                    ind_don1.append(ii)
                    dind_Hdon1[ii]=system1.donors_hydrogen[ii]
                    ind_Hdon1.append(system1.donors_hydrogen[ii])
                    num_Hdon1.append(len(system1.donors_hydrogen[ii]))
                    ind_acc1.append(ii)
                elif system1.atom[ii].donor:
                    ind_don2.append(ii)
                    dind_Hdon2[ii]=system1.donors_hydrogen[ii]
                    ind_Hdon2.append(system1.donors_hydrogen[ii])
                    num_Hdon2.append(len(system1.donors_hydrogen[ii]))
                elif system1.atom[ii].acceptor:
                    ind_acc2.append(ii)

            max_Hdon1=max(num_Hdon1)
            max_Hdon2=max(num_Hdon2)
            same_set=True
        
        Ldon1=len(ind_don1)
        Lacc1=len(ind_acc1)
        Ldon2=len(ind_don2)
        Lacc2=len(ind_acc2)


    if same_set:

        f_hb.hbonds.ldon1=Ldon1
        f_hb.hbonds.lacc1=Lacc1
        f_hb.hbonds.ldon2=Ldon2
        f_hb.hbonds.lacc2=Lacc2
        f_hb.hbonds.maxh1=max_Hdon1
        f_hb.hbonds.maxh2=max_Hdon2

        f_hb.hbonds.initialize_memory()
        f_hb.hbonds.initialize_pbc(box)

        f_hb.hbonds.ind_don1=array(ind_don1,dtype=int,order='Fortran')
        f_hb.hbonds.num_hdon1=array(num_Hdon1,dtype=int,order='Fortran')
        f_hb.hbonds.ind_acc1=array(ind_acc1,dtype=int,order='Fortran')
        f_hb.hbonds.ind_don2=array(ind_don2,dtype=int,order='Fortran')
        f_hb.hbonds.num_hdon2=array(num_Hdon2,dtype=int,order='Fortran')
        f_hb.hbonds.ind_acc2=array(ind_acc2,dtype=int,order='Fortran')
       
        # I should get rid of this
        for jj in range(Ldon1):
            f_hb.hbonds.coor_don1[jj,:]=system1.frame[0].coors[ind_don1[jj],:]
            for kk in range(num_Hdon1[jj]):
                bb=ind_Hdon1[jj][kk]
                f_hb.hbonds.coor_hdon1[jj,kk,:]=system1.frame[0].coors[bb,:]
        for jj in range(Lacc1):
            f_hb.hbonds.coor_acc1[jj,:]=system1.frame[0].coors[ind_acc1[jj],:]
         
        for jj in range(Ldon2):
            f_hb.hbonds.coor_don2[jj,:]=system1.frame[0].coors[ind_don2[jj],:]
            for kk in range(num_Hdon2[jj]):
                bb=ind_Hdon2[jj][kk]
                f_hb.hbonds.coor_hdon2[jj,kk,:]=system1.frame[0].coors[bb,:]
        for jj in range(Lacc2):
            f_hb.hbonds.coor_acc2[jj,:]=system1.frame[0].coors[ind_acc2[jj],:]
         
        # I should correct the use of pointers in the function
        # apparently is slow.
        f_hb.hbonds.same_set(r_param,ang_param)
         
        ## I should get rid of this 
        for jj in range(f_hb.hbonds.num_hbs1):
            aa=f_hb.hbonds.salida1[jj,0]
            bb=f_hb.hbonds.salida1[jj,1]
            f_hb.hbonds.salida1[jj,1]=dind_Hdon1[aa][bb]
         
        for jj in range(f_hb.hbonds.num_hbs2):
            aa=f_hb.hbonds.salida2[jj,0]
            bb=f_hb.hbonds.salida2[jj,1]
            f_hb.hbonds.salida2[jj,1]=dind_Hdon1[aa][bb]
         
        for jj in range(f_hb.hbonds.num_hbs3):
            aa=f_hb.hbonds.salida3[jj,0]
            bb=f_hb.hbonds.salida3[jj,1]
            f_hb.hbonds.salida3[jj,1]=dind_Hdon2[aa][bb]
         
        for jj in range(f_hb.hbonds.num_hbs4):
            aa=f_hb.hbonds.salida4[jj,0]
            bb=f_hb.hbonds.salida4[jj,1]
            f_hb.hbonds.salida4[jj,1]=dind_Hdon2[aa][bb]
         
         
         
        output1=copy.deepcopy(f_hb.hbonds.salida1)
        output2=copy.deepcopy(f_hb.hbonds.salida2)
        output3=copy.deepcopy(f_hb.hbonds.salida3)
        output4=copy.deepcopy(f_hb.hbonds.salida4)
         
        f_hb.hbonds.free_memory()
        #print 'si'
#        return output1,output2,output3,output4
        return f_hb.hbonds.num_hbs1, f_hb.hbonds.num_hbs2, f_hb.hbonds.num_hbs3, f_hb.hbonds.num_hbs4

    else:

        f_hb.hbonds.ldon1=Ldon1
        f_hb.hbonds.lacc1=Lacc1
        f_hb.hbonds.ldon2=Ldon2
        f_hb.hbonds.lacc2=Lacc2
        f_hb.hbonds.maxh1=max_Hdon1
        f_hb.hbonds.maxh2=max_Hdon2

        f_hb.hbonds.initialize_memory()
        f_hb.hbonds.initialize_pbc(box)
        
        f_hb.hbonds.ind_don1=array(ind_don1,dtype=int,order='Fortran')
        f_hb.hbonds.num_hdon1=array(num_Hdon1,dtype=int,order='Fortran')
        f_hb.hbonds.ind_acc1=array(ind_acc1,dtype=int,order='Fortran')
        f_hb.hbonds.ind_don2=array(ind_don2,dtype=int,order='Fortran')
        f_hb.hbonds.num_hdon2=array(num_Hdon2,dtype=int,order='Fortran')
        f_hb.hbonds.ind_acc2=array(ind_acc2,dtype=int,order='Fortran')

        for jj in range(Ldon1):
            f_hb.hbonds.coor_don1[jj,:]=system1.frame[0].coors[ind_don1[jj],:]
            for kk in range(num_Hdon1[jj]):
                bb=ind_Hdon1[jj][kk]
                f_hb.hbonds.coor_hdon1[jj,kk,:]=system1.frame[0].coors[bb,:]
        for jj in range(Lacc1):
            f_hb.hbonds.coor_acc1[jj,:]=system1.frame[0].coors[ind_acc1[jj],:]

        for jj in range(Ldon2):
            f_hb.hbonds.coor_don2[jj,:]=system2.frame[0].coors[ind_don2[jj],:]
            for kk in range(num_Hdon2[jj]):
                bb=ind_Hdon2[jj][kk]
                f_hb.hbonds.coor_hdon2[jj,kk,:]=system2.frame[0].coors[bb,:]
        for jj in range(Lacc2):
            f_hb.hbonds.coor_acc2[jj,:]=system2.frame[0].coors[ind_acc2[jj],:]


        f_hb.hbonds.diff_set(r_param,ang_param)

        for jj in range(f_hb.hbonds.num_hbs1):
            aa=f_hb.hbonds.salida1[jj,0]
            bb=f_hb.hbonds.salida1[jj,1]
            f_hb.hbonds.salida1[jj,1]=dind_Hdon1[aa][bb]

        for jj in range(f_hb.hbonds.num_hbs2):
            aa=f_hb.hbonds.salida2[jj,0]
            bb=f_hb.hbonds.salida2[jj,1]
            f_hb.hbonds.salida2[jj,1]=dind_Hdon2[aa][bb]

        output1=copy.deepcopy(f_hb.hbonds.salida1)
        output2=copy.deepcopy(f_hb.hbonds.salida2)

        f_hb.hbonds.free_memory()

    return f_hb.hbonds.num_hbs1, output1, f_hb.hbonds.num_hbs2, output2

