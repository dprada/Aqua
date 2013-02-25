from numpy import *
import pyn_fort_enm as f_enm
import pyn_fort_general as f
from pyn_cl_set import *
import pylab


#####################################################################################
##### Gaussian Network Model
#####################################################################################


class gnm():
    
    def __init__(self,system=None,cutoff=10.0):

        self.contact_map=None

        if system!=None:
            if type(system) in [list,ndarray]:
                self.contact_map=array(contact_map,order='Fortran')
            else:
                self.contact_map=self.mk_contact_map(system,cutoff)
            self.eigenvals,self.eigenvects,self.freqs,self.bfacts,self.inverse,self.correl=self.rebuild()
            self.bfacts_pdb,self.factor,self.sqr_dev=self.fitt_bfacts()
            
    def rebuild(self):

        return f_enm.gnm(self.contact_map,len(self.contact_map[0]))
        

    def mk_contact_map(self,system,cutoff):

        self.system=system
        comap=f_enm.contact_map(cutoff,system.frame[0].coors,system.num_atoms)
        return comap

    def fitt_bfacts(self):
        
        bfacts_pdb=[]
        for ii in range(len(self.contact_map)):
            bfacts_pdb.append(self.system.atom[ii].bfactor)

        aa=0.0
        bb=0.0

        for ii in range(len(self.contact_map)):
            aa+=bfacts_pdb[ii]*self.bfacts[ii]
            bb+=self.bfacts[ii]*self.bfacts[ii]

        aa=aa/bb

        bb=0.0
        for ii in range(len(self.contact_map)):
            bb+=(bfacts_pdb[ii]-aa*self.bfacts[ii])**2
            
        return bfacts_pdb,aa,bb

    def plot_best_cutoff(self):

        ctoff=[]
        r_2=[]
        l=1.0*len(self.system.atom)
        for ii in arange(6.5,12.6,0.1):
            ctoff.append(ii)
            aa=gnm_classic(self.system,cutoff=ii)
            r_2.append((aa.sqr_dev)/l)
            del(aa)

        pylab.plot(ctoff,r_2,'yo')
        pylab.ylabel('<R^2>|atom')
        pylab.xlabel('Cut Off (A)')
        self.best_cutoff=[]
        self.best_cutoff.append(ctoff)
        self.best_cutoff.append(r_2)
        return pylab.show()

    def plot_bfacts(self):

        pylab.plot(self.bfacts_pdb,color="blue")
        pylab.plot(self.factor*self.bfacts,color="red")
        return pylab.show()

    def plot_dispersion_bfacts(self):

        pylab.plot(self.bfacts,self.bfacts_pdb,'yo')
        pylab.plot(self.bfacts,self.factor*self.bfacts,'r--')
        return pylab.show()

    def plot_contact_map(self):

        pylab.gray()
        pylab.matshow(self.contact_map,cmap='binary')
        return pylab.show()

    def plot_inverse(self):

        pylab.gray()
        pylab.matshow(self.inverse,cmap='binary')
        return pylab.show()
    
    def plot_correl_norm_2(self):
        
        #pylab.imshow(self.correl,origin='lower',interpolation=None) 
        vmin=ma.minimum(self.correl)
        vmax=ma.maximum(self.correl)
        vmax=max([abs(vmin),vmax])
        #ref_white=(-vmin)/(vmax-vmin)
        cdict = {                                                                            
            'red'  :  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,1.0,1.0)),
            'green':  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,0.0,0.0)),
            'blue' :  ((0.0,1.0,1.0), (0.5,1.0,1.0), (1.0,0.0,0.0))
            }
        my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        pylab.matshow(self.correl,cmap=my_cmap,vmin=-vmax,vmax=vmax)
        pylab.colorbar()
        return pylab.show()

    def plot_correl_norm(self):
        
        cdict = {                                                                            
            'red'  :  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,1.0,1.0)),
            'green':  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,0.0,0.0)),
            'blue' :  ((0.0,1.0,1.0), (0.5,1.0,1.0), (1.0,0.0,0.0))
            }
        my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        #pylab.matshow(self.correl,cmap='RdBu',vmin=-1.0,vmax=1.0)
        #pylab.pcolor(self.correl,cmap=my_cmap,vmin=-1.0,vmax=1.0)
        #nx,ny=self.correl.shape
        #pylab.matshow(self.correl,cmap=my_cmap,vmin=-1.0,vmax=1.0,extent=[0,nx,0,ny])
        pylab.matshow(self.correl,cmap=my_cmap,vmin=-1.0,vmax=1.0)
        #if hasattr(self,'system'):
        #    pylab.xticks(arange(0,len(self.system.atom[:]),5),[x.resid_pdb_index for x in self.system.atom[:]],rotation=90)
        #    pylab.yticks(arange(0,len(self.system.atom[:]),5),[x.resid_pdb_index for x in self.system.atom[:]])
        pylab.colorbar()
        return pylab.show()

    def plot_correl(self):
        
        vmin=ma.minimum(self.inverse)
        vmax=ma.maximum(self.inverse)     
        vmax=max([abs(vmin),vmax])
        cdict = {                                                                            
            'red'  :  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,1.0,1.0)),
            'green':  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,0.0,0.0)),
            'blue' :  ((0.0,1.0,1.0), (0.5,1.0,1.0), (1.0,0.0,0.0))
            }
        my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        pylab.matshow(self.inverse,cmap=my_cmap,vmin=-vmax,vmax=vmax)
        pylab.colorbar()
        return pylab.show()


    def write(self):

        f_map = open('contact_map.oup','w')
        
        for ii in range(len(self.contact_map)):
            for jj in range(ii+1,len(self.contact_map)):
                if self.contact_map[ii][jj] == True :
                    f_map.write("%s %s \n" %(self.system.atom[ii].pdb_index,self.system.atom[jj].pdb_index))
        f_map.close

        f_vects = open('gnm_vects.oup','w')

        f_vects.write("%s Modes, %s Nodes \n" %(len(self.contact_map),len(self.contact_map)))
        f_vects.write(" \n")
        for ii in range(len(self.contact_map)):
            for jj in range(len(self.contact_map)):
                f_vects.write("%s %f \n" %(self.system.atom[jj].pdb_index,self.eigenvects[ii][jj]))
            f_vects.write(" \n")




#####################################################################################
##### Anisotropic Network Model
#####################################################################################

class anm():
    
    def __init__(self,system=None,cutoff=10.0,contact_map=None):

        if system==None and contact_map==None:                 # Condition in case of error
            print 'Provide a system or a contact map:'
            print 'anm(system=foo,cutoff=10.0)'
            print 'anm(contact_map=matrix(real))'
            return

        self.num_nodes=None       # Number of nodes
        self.num_modes=None       # Number of modes
        self.contact_map=None     # Contact map (real matrix= kij; kij=force constant ij)
        self.parent=system        # System analysed
        self.eigenvals=None       # eigenvalues
        self.eigenvects=None      # eigenvectors as they come from the analysis vector: [x1,y1,z1,x2,y2,z2...]
        self.eigenvects_3d=None   # eigenvectors 3d: Nx3 [[x1,y1,z1],[x2,y2,z2]...]
        self.freqs=None           # frequencies (eigenvectors dimensionalized -fitting bfactors-)
        self.inverse=None         # inverse matrix from single value decomposition
        self.correl=None          # correlation matrix
        self.node=[]              # List of nodes

        for atom in system.atom:
            prov_node=labels_unit()
            prov_node.name=atom.name
            prov_node.index=atom.index
            prov_node.pdb_index=atom.pdb_index
            self.node.append(prov_node)

        self.system=system
        if contact_map!=None:                                  # Building the contact map
            self.contact_map=array(contact_map,order='Fortran')
        else:
            self.contact_map=self.make_contact_map(system,cutoff)

        self.num_nodes=len(self.contact_map)                   # Number of nodes
        self.eigenvals,self.eigenvects,self.freqs,self.bfacts,self.inverse,self.correl=self.build(system) # analysis

        self.eigenvects_3d=zeros(shape=(len(self.contact_map)*3,len(self.contact_map),3))  # eigenvectors_3d
        for aa in range(3*len(self.contact_map)):
            for ii in range(len(self.contact_map)):
                iii=(ii)*3
                for jj in range(3):
                    jjj=iii+jj
                    self.eigenvects_3d[aa,ii,jj]=self.eigenvects[aa,jjj]

        self.bfacts_pdb,self.factor,self.sqr_dev=self.fitt_bfacts()

    def build(self,system):

        return f_enm.anm(self.contact_map,system.frame[0].coors,len(self.contact_map[0]))

    def rebuild(self):

        self.eigenvals,self.eigenvects,self.freqs,self.bfacts,self.inverse,self.correl=self.build(self.parent)
        self.eigenvects_3d=zeros(shape=(len(self.contact_map)*3,len(self.contact_map),3))
        for aa in range(3*len(self.contact_map)):
            for ii in range(len(self.contact_map)):
                iii=(ii)*3
                for jj in range(3):
                    jjj=iii+jj
                    self.eigenvects_3d[aa,ii,jj]=self.eigenvects[aa,jjj]

        self.bfacts_pdb,self.factor,self.sqr_dev=self.fitt_bfacts()
        

    def make_contact_map(self,system,cutoff):

        comap=f_enm.contact_map(cutoff,system.frame[0].coors,system.num_atoms)
        return comap

    def fitt_bfacts(self):
        
        bfacts_pdb=[]
        for ii in range(len(self.contact_map)):
            bfacts_pdb.append(self.system.atom[ii].bfactor)

        aa=0.0
        bb=0.0

        for ii in range(len(self.contact_map)):
            aa+=bfacts_pdb[ii]*self.bfacts[ii]
            bb+=self.bfacts[ii]*self.bfacts[ii]

        aa=aa/bb

        bb=0.0
        for ii in range(len(self.contact_map)):
            bb+=(bfacts_pdb[ii]-aa*self.bfacts[ii])**2
            
        return bfacts_pdb,aa,bb

    def best_cutoff(self):

        ctoff=[]
        r_2=[]
        l=1.0*len(self.system.atom)
        for ii in arange(8.0,16.0,0.2):
            ctoff.append(ii)
            aa=anm_classic(self.system,cutoff=ii)
            r_2.append((aa.sqr_dev)/l)
            del(aa)

        pylab.plot(ctoff,r_2,'yo')
        pylab.ylabel('<R^2>|atom')
        pylab.xlabel('Cut Off (A)')
        self.best_cutoff=[]
        self.best_cutoff.append(ctoff)
        self.best_cutoff.append(r_2)
        return pylab.show()

    def plot_bfacts(self):

        pylab.plot(self.bfacts_pdb,color="blue")
        pylab.plot(self.factor*self.bfacts,color="red")
        return pylab.show()

    def plot_dispersion_bfacts(self):

        pylab.plot(self.bfacts,self.bfacts_pdb,'yo')
        pylab.plot(self.bfacts,self.factor*self.bfacts,'r--')
        return pylab.show()

    def plot_contact_map(self):

        pylab.gray()
        pylab.matshow(self.contact_map,cmap='binary')
        return pylab.show()

    def plot_inverse(self):

        pylab.gray()
        pylab.matshow(self.inverse,cmap='binary')
        return pylab.show()
    
    def plot_correl_norm_2(self):
        
        #pylab.imshow(self.correl,origin='lower',interpolation=None) 
        vmin=ma.minimum(self.correl)
        vmax=ma.maximum(self.correl)
        vmax=max([abs(vmin),vmax])
        #ref_white=(-vmin)/(vmax-vmin)
        cdict = {                                                                            
            'red'  :  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,1.0,1.0)),
            'green':  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,0.0,0.0)),
            'blue' :  ((0.0,1.0,1.0), (0.5,1.0,1.0), (1.0,0.0,0.0))
            }
        my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        pylab.matshow(self.correl,cmap=my_cmap,vmin=-vmax,vmax=vmax)
        pylab.colorbar()
        return pylab.show()

    def plot_correl_norm(self):
        
        cdict = {                                                                            
            'red'  :  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,1.0,1.0)),
            'green':  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,0.0,0.0)),
            'blue' :  ((0.0,1.0,1.0), (0.5,1.0,1.0), (1.0,0.0,0.0))
            }
        my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        #pylab.matshow(self.correl,cmap='RdBu',vmin=-1.0,vmax=1.0)
        #pylab.pcolor(self.correl,cmap=my_cmap,vmin=-1.0,vmax=1.0)
        #nx,ny=self.correl.shape
        #pylab.matshow(self.correl,cmap=my_cmap,vmin=-1.0,vmax=1.0,extent=[0,nx,0,ny])
        pylab.matshow(self.correl,cmap=my_cmap,vmin=-1.0,vmax=1.0)
        #if hasattr(self,'system'):
        #    pylab.xticks(arange(0,len(self.system.atom[:]),5),[x.resid_pdb_index for x in self.system.atom[:]],rotation=90)
        #    pylab.yticks(arange(0,len(self.system.atom[:]),5),[x.resid_pdb_index for x in self.system.atom[:]])
        pylab.colorbar()
        return pylab.show()

    def plot_correl(self):
        
        vmin=ma.minimum(self.inverse)
        vmax=ma.maximum(self.inverse)     
        vmax=max([abs(vmin),vmax])
        cdict = {                                                                            
            'red'  :  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,1.0,1.0)),
            'green':  ((0.0,0.0,0.0), (0.5,1.0,1.0), (1.0,0.0,0.0)),
            'blue' :  ((0.0,1.0,1.0), (0.5,1.0,1.0), (1.0,0.0,0.0))
            }
        my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        pylab.matshow(self.inverse,cmap=my_cmap,vmin=-vmax,vmax=vmax)
        pylab.colorbar()
        return pylab.show()


    def build_correl(self,modes='all'):

        list_modes=[]
        if modes=='all':
            num_modes=3*len(self.contact_map)
            for ii in range(num_modes):
                list_modes.append(ii)
        elif type(modes)==int:
            num_modes=1
            list_modes.append(modes)
        elif type(modes) in [list,tuple]:
            num_modes=len(modes)
            list_modes=list(modes)

        self.correl=f_enm.correlation(self.eigenvects,self.eigenvals,list_modes,self.num_nodes,num_modes)

    def involv_coefficient(self,modes='all',vect=None):

        if vect==None:
            print 'Error: vect=None'
            return

        self.ic={}
        list_modes=[]
        
        if modes=='all':
            for ii in range(self.num_modes):
                list_modes.append(ii+1)
        elif type(modes)==int:
            list_modes.append(modes)
        elif type(modes) in [list,tuple]:
            list_modes=modes

        for ii in list_modes:
            self.ic[ii]=f.aux_funcs_general.proj3d(self.eigenvects_3d[ii],vect,len(self.eigenvects_3d[ii]))

        return

        

    def write(self):

        f_map = open('contact_map.oup','w')
        
        for ii in range(len(self.contact_map)):
            for jj in range(ii+1,len(self.contact_map)):
                if self.contact_map[ii][jj] == True :
                    f_map.write("%s %s \n" %(self.system.atom[ii].pdb_index,self.system.atom[jj].pdb_index))
        f_map.close

        f_vects = open('anm_vects.oup','w')

        f_vects.write("%s Modes, %s Nodes \n" %(len(self.contact_map),len(self.contact_map)))
        f_vects.write(" \n")
        for aa in range(3*len(self.contact_map)):
            for ii in range(len(self.contact_map)):
                    f_vects.write("%s %f %f %f\n" %(self.system.atom[jj].pdb_index,self.eigenvects_3d[aa,ii,0],self.eigenvects_3d[aa,ii,1],self.eigenvects_3d[aa,ii,2]))
                    f_vects.write(" \n")




#####################################################################################
##### Building pdb movie from Anisotropic Network Model
#####################################################################################

def build_fluct_anm(system,anm,mode='all',output=None,amplitude=8.0,steps=60):

    prov_list=[]
    for ii in system.atom:
        if ii.name in ['N','CA','C','O'] and ii.type_pdb =='ATOM':
            prov_list.append(ii.index)

    prov_system=make_selection(system,prov_list)

    for aa in prov_system.atom[:]:
        if aa.type_pdb in ['ATOM']:
            if aa.chain.name not in prov_system.chains:
                prov_system.chains.append(aa.chain.name)

    ## building a reduced eigenvects if this is necessary:

    num_modes=len(anm.eigenvects_3d[:])    
    jj=0
    for ii in range(anm.num_nodes):
        if anm.node[ii].index in prov_list:
            jj+=1

    aux_vect_3d=npy.zeros(shape=(num_modes,jj,3),order='Fortran')

    jj=-1
    for ii in range(anm.num_nodes):
        if anm.node[ii].index in prov_list:
            jj+=1
            aux_vect_3d[:,jj,:]=anm.eigenvects_3d[:,ii,:]

    ## The 'anm.eigenvects_3d' has been replaced by aux_vect_3d[:,jj,:]
    ## since this point.

    list_modes=[]

    if mode=='all':
        for ii in range(num_modes):
            list_modes.append(ii)
    elif type(mode)==int:
        num_modes=1
        list_modes.append(mode)
    elif type(mode) in [list,tuple]:
        num_modes=len(mode)
        list_modes=list(mode)
    
    osc=zeros(shape=(prov_system.num_atoms,3))

    prefix=system.name
    if prefix[-1]=='.':
        prefix=prefix[:-1]

    in_net=[]
    for ii in anm.system.atom:
        if ii.type_pdb in ['ATOM']: ## Added
            in_net.append(ii.index)
    in_syst=[]
    for ii in prov_system.atom:
        in_syst.append(ii.index)

    tt=zeros(shape=(prov_system.num_atoms))


    jj=-1

    for chch in prov_system.chains :
        interr=-1
        for ii in anm.system.atom:
            if ii.type_pdb in ['ATOM']: ## Added
                if ii.chain.name == chch:
                    extreme=ii.index
                    extreme=in_net.index(extreme)

        for ii in prov_system.atom:
            if ii.chain.name == chch:
                jj+=1
        
                if ii.index in in_net:
                    interr+=1
                    net_initial=in_net.index(ii.index)
                    net_end=net_initial+1
                    if net_end>extreme:
                        interr=-1
                    else:
                        initial=jj
                        end=in_syst.index(in_net[net_end]) #

                if interr==-1:
                    tt[jj]=1.0
                else:
                    tt[jj]=dot((prov_system.frame[0].coors[jj]-prov_system.frame[0].coors[initial]),(prov_system.frame[0].coors[end]-prov_system.frame[0].coors[initial]))
                    tt[jj]=tt[jj]/(dot((prov_system.frame[0].coors[end]-prov_system.frame[0].coors[initial]),(prov_system.frame[0].coors[end]-prov_system.frame[0].coors[initial])))


    for ind_mode in list_modes :
        kk=0
        jj=-1
        for chch in prov_system.chains :
            interr=-1
            for ii in anm.system.atom:
                if ii.type_pdb in ['ATOM']: ## Added
                    if ii.chain.name == chch:
                        extreme=ii.index
                        extreme=in_net.index(extreme)
            ant=aux_vect_3d[ind_mode][kk]
            for ii in prov_system.atom:
                if ii.chain.name == chch:
                    jj+=1
        
                    if ii.index in in_net:
                        interr+=1
                        net_initial=in_net.index(ii.index)
                        net_end=net_initial+1
                        if net_end>extreme:
                            interr=-1
                            ant[:]=aux_vect_3d[ind_mode][kk]
                        else:
                            aaa=aux_vect_3d[ind_mode][net_initial]
                            bbb=aux_vect_3d[ind_mode][net_end]
                        kk+=1
                    if interr==-1:
                        osc[jj][:]=ant[:]
                    else:
                        osc[jj][:]=aaa[:]+(bbb[:]-aaa[:])*tt[jj]
            

        if output==None:
            file_name=prefix+'_anm_'+str(ind_mode)+'.pdb'
        else:
            file_name=output
        file=open(file_name,'w')

        a='HEADER    '+prefix+'     ANM: Mode '+str(ind_mode)+'\n'
        file.write(str(a))

        for ii in system.pdb_ss:
            file.write(str(ii))
 

        delta_f=2.0*pi/(steps*1.0)

        for frame in range(0,steps):

            a='MODEL '+str(frame)+'\n'
            file.write(str(a))


            for ii in range(prov_system.num_atoms):
                a='ATOM  '                                 # 1-6
                a+="%5d" % (ii+1)                          # 7-11
                #a+="%5d" % prov_system.atom[ii].pdb_index  # 7-11
                a+=' '                                     # 12
                a+=' '+"%-3s" % prov_system.atom[ii].name  # 13-16
                a+=' '                                     # 17
                a+="%3s" % prov_system.atom[ii].resid.name # 18-20
                a+=' '                                     # 21
                a+="%1s" % prov_system.atom[ii].chain.name # 22
                a+="%4d" % prov_system.atom[ii].resid.pdb_index # 23-26
                a+=' '                                     # 27
                a+='   '                                   # 28-30
                a+="%8.3f" % float(prov_system.frame[0].coors[ii][0]+amplitude*sin(delta_f*frame)*osc[ii][0]) # 31-38
                a+="%8.3f" % float(prov_system.frame[0].coors[ii][1]+amplitude*sin(delta_f*frame)*osc[ii][1]) # 39-46
                a+="%8.3f" % float(prov_system.frame[0].coors[ii][2]+amplitude*sin(delta_f*frame)*osc[ii][2]) # 47-54
                a+="%6.2f" % prov_system.atom[ii].occup    # 55-60
                a+="%6.2f" % prov_system.atom[ii].bfactor  # 61-66
                a+='          '                            # 67-76
                a+="%2s" % prov_system.atom[ii].elem_symb  # 77-78
                a+="%2s" % prov_system.atom[ii].charge     # 79-80

                if a[-1]!='\n' :
                    a+='\n'
                file.write(str(a))

            a='ENDMDL \n'
            file.write(str(a))

        file.close() 


    return
