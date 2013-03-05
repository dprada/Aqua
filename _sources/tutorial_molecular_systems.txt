
Tutorial Molecular Systems
**************************

First of all, lets load Pynoramix in our script or in a ipython session:

.. sourcecode:: ipython

     In [1]: from pynoramix_beta import *


Some basic notions on python will be assumed along this tutorial. If you just landed here without any idea on python, have a look to the section *First steps on python*.

.. todo:: Make a short tutorial on python, enough to run pynoramix.

----------------------
 

Loading/Writting the topology
=============================

Loading
+++++++

A system can be loaded from a file (pdb,gro) or downloaded from the Protein Data Bank.

.. sourcecode:: ipython

   In [2]: mol_test=molecule(download='2WC2')
   # File saved as 2WC2.pdb
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions
   # 20  frames/models in traj 0

   In [3]: mol_test=molecule('2WC2.pdb',verbose=False)

   In [4]: mol_test.info()
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions
   # 20  frames/models in traj 0


.. todo:: Complete the topology of residues and terminals. Apparently
   it works because the option "with_bonds=False" is by default. It needs to be watched.

This way the pdb file has been loaded together with the coordinates
present in the file, 20 different models in the case of 2WC2.  If the
coordinates are going to be loaded from a different file, the topology
can be created adding the option coors=False:

.. sourcecode:: ipython

   In [5]: mol_test2=molecule('2WC2.pdb',coors=False)
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions

   In [6]: mol_test2.info_trajs()
   # No coordinates

   In [7]: mol_test.info_trajs()
   # 20 frames/models in traj 0


Navigating
++++++++++



Writting
++++++++


sdfsdfdsf

----------------------

Loading/Writting a MD trajectory
================================

Loading
+++++++

Once a topology has been created a trajectory can be loaded from
different formats: pdb, gro, xtc, trr, dcd, bin (to be deprecated).

It is recommended the use of dcd files, the file is unformatted and
thereby it is small and easy to handle.

Along this section the different ways to do it will be illustrated
using the files :download:`GSGS.pdb <../tutorials/systems_tut1/GSGS.pdb>`
and :download:`GSGS.dcd <../tutorials/systems_tut1/GSGS.dcd>`.

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb')
   # System created from the file GSGS.pdb :
   # 4723  atoms
   # 1568  residues
   # 3  chains
   # 1560  waters
   # 4  ions
   # 1  frames/models in traj 0

   In [3]: GSGS.delete_traj()
    
   In [4]: GSGS.info_trajs()
   # No coordinates

   In [5]: GSGS.load_traj('GSGS.dcd','ALL')
   # 10 frames/models loaded.

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',coors=False,verbose=False)
    
   In [3]: GSGS.load_traj('GSGS.dcd',frame='ALL',verbose=False)
    
   In [4]: GSGS.info(); GSGS.info_trajs()
   # System created from the file GSGS.pdb :
   # 4723  atoms
   # 1568  residues
   # 3  chains
   # 1560  waters
   # 4  ions
   # 10  frames/models in traj 0

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',coors=False,verbose=False)

   In [3]: GSGS.load_traj('GSGS.dcd')
   # 0 frames/models in traj 0

   In [4]: print GSGS.traj[0].name, GSGS.traj[0].type, GSGS.traj[0].io_opened, GSGS.traj[0].io_end
   GSGS.dcd dcd True False

   In [5]: while 1:
     ....:     GSGS.traj[0].upload_frame()
     ....:     if GSGS.traj[0].io_end: break
     ....: 
   # End of file

   In [6]: GSGS.info_trajs()
   # 10 frames/models in traj 0

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',coors=False,verbose=False)

   In [3]: GSGS.load_traj('GSGS.dcd',frame=0)  # Or frame='Next'
   # 1 frames/models in traj 0

   In [4]: while GSGS.traj[0].io_opened:
      ...:     print GSGS.traj[0].frame[0].coors[0]
      ...:     GSGS.traj[0].reload_frame()
      ...: 
   [ -7.26851273  -8.12112617  10.57811832]
   [ -5.16595078  -9.8920269   12.24640751]
   [ -6.12880325  -9.20014763  15.28322697]
   [ -4.90646744  -8.31535339  12.97708988]
   [ -5.04781723  -9.68705559  14.15655327]
   [ -5.95707321  -8.45479965  17.51550102]
   [ -4.45994186 -10.63479614  16.19140053]
   [ -6.01659775 -13.60509872  16.98220253]
   [ -4.40946579 -13.10482597  17.12298393]
   [ -5.01924515 -13.77911949  15.64630699]
   # End of file

   In [5]: GSGS.info_trajs()
   # 1 frames/models in traj 0





Converting a trajectory into other format
+++++++++++++++++++++++++++++++++++++++++

Right now the output formats are only dcd files.

This way the original trajectory is stored in memory:

.. sourcecode:: ipython

   In [2]: ion=molecule('run_ion.gro',coors=False,verbose=False)
    
   In [3]: ion.load_traj('traj.xtc',frame='ALL',verbose=False)
    
   In [4]: ion.traj[0].write('new_traj.dcd',action='Open')
    
   In [5]: ion.traj[0].write(frame='ALL')
    
   In [6]: ion.traj[0].write(action='Close')

This way the original trajectory is not stored in memory:

.. sourcecode:: ipython

   In [2]: ion=molecule('run_ion.gro',coors=False,verbose=False)
    
   In [3]: ion.load_traj('traj.xtc',frame='Next',verbose=False)
    
   In [4]: ion.traj[0].write('new_traj.dcd',action='Open')
    
   In [5]: while ion.traj[0].io_opened:
      ...:     ion.traj[0].write()
      ...:     ion.traj[0].reload_frame()
      ...: 
   # End of file
    
   In [6]: ion.traj[0].write(action='Close')


How to make atoms selections
============================

The syntax is close to the pynoramix syntax.
There are few special key words.

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',verbose=False)
    
   In [3]: list1=GSGS.selection('backbone')

   In [4]: list2=GSGS.selection('atom.name N CA C O')
    
   In [5]: print list1; print list2
   [0, 4, 7, 8, 9, 11, 18, 19, 20, 22, 25, 26, 27, 30, 32]
   [0, 4, 7, 8, 9, 11, 18, 19, 20, 22, 25, 26, 27, 30, 32]

   In [5]: list1=GSGS.selection('sidechain')
   
   In [6]: list2=GSGS.selection('(atom.resid.type Protein and not atom.name N CA C O H1 H2)')
    
   In [7]: print list; print list2
   [1, 2, 3, 5, 6, 10, 12, 13, 14, 15, 16, 17, 21, 23, 24, 28, 29, 31, 33, 34, 35, 36, 37, 38]
   [1, 2, 3, 5, 6, 10, 12, 13, 14, 15, 16, 17, 21, 23, 24, 28, 29, 31, 33, 34, 35, 36, 37, 38]

We can also make use of the expression 'within X of', X is a float number indicating a distance threshold.


.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',verbose=False)
    
   In [3]: list1=GSGS.selection('atom.name OH2 within 3.0 of atom.resid.type Protein')

.. sourcecode:: ipython

   In [2]: ion=molecule('run_ion.gro',coors=False,verbose=False) 
    
   In [3]: ion.load_traj('traj.dcd',frame='ALL',verbose=False)
    
   In [4]: list1=ion.selection('atom.name OW within 3.0 of atom.resid.type Ion')

   In [5]: for ii in range(3):
      ...:     print len(list1[ii]),'waters below 3.0 in frame', ii
      ...: 
   6 waters below 3.0 in frame 0
   6 waters below 3.0 in frame 1
   6 waters below 3.0 in frame 2


Computing distances
===================

The distance between a set of n1 atoms, list1, and a set of n2 atoms,
list2.  If only one frame is analysed the output is a 2D matrix
[n1,n2].  This way the distance between the i-th atom in list1 and
j-th in list2 correspond to output[i,j].

If more than one frame is analysed the output is a 3D matrix
[n1,n2,number_frames] (following the previous notation).

It is faster if len(list1)<len(list2).


.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',coors=False,verbose=False)
    
   In [3]: GSGS.load_traj('GSGS.dcd',frame='ALL',verbose=False)
    
   In [4]: list1=GSGS.selection('atom.resid.type Protein')
    
   In [5]: list2=GSGS.selection('atom.resid.type Water and atom.type O')

   In [6]: result=GSGS.distance(list1,list2)
 
   In [7]: for ii in range(GSGS.traj[0].num_frames):
      ....:    print 'The distance between atoms index',list1[1],'and',list2[3],'is',result[1,3][ii],'in frame',ii
      ....: 
   The distance between atoms index 1 and 48 is 24.8435076017 in frame 0
   The distance between atoms index 1 and 48 is 23.6529328175 in frame 1
   The distance between atoms index 1 and 48 is 24.3209230117 in frame 2
   The distance between atoms index 1 and 48 is 21.5236312048 in frame 3
   The distance between atoms index 1 and 48 is 25.2685193116 in frame 4
   The distance between atoms index 1 and 48 is 28.2550958504 in frame 5
   The distance between atoms index 1 and 48 is 26.1290587977 in frame 6
   The distance between atoms index 1 and 48 is 20.9157208891 in frame 7
   The distance between atoms index 1 and 48 is 21.6473615840 in frame 8
   The distance between atoms index 1 and 48 is 18.5862638499 in frame 9

Neighbors and Ranked contacts
+++++++++++++++++++++++++++++

The function neighbs can help with its different options to approach this problems.
Notice that the cut-off here is the limit in the ranked list of closest neighbors or a given distance.
In the former case the output can be sorted or not by distance. 
If only contact map is need, maybe the following section is suitable because of its efficience.

Neighbors with a distance lower or equal than a certain threshold:



Contact Maps
============



Radial Distribution Funcions
============================

The radial distribution function g_{a,b}(r) represents the radial
concentration of atoms "b" around atoms "a", normalized by the
concentration of "b".


It is more efficient (fast and no memory consumming) when the trajectorie is read frame by frame, and
not loaded at a time.
As it happens with the distance function, when len(list1)<len(list2) it is faster.

.. warning:: Right now the function does not work properly if
   setA=setB. In addition, this should be efficient including the
   condition "same set" for function dists.

.. sourcecode:: ipython

   In [2]: ion=molecule('run_ion.gro',coors=False,verbose=False) 
    
   In [3]: ion.load_traj('traj.dcd',frame='ALL',verbose=False)
    
   In [4]: list1=ion.selection('atom.resid.type Ion')
    
   In [5]: list2=ion.selection('atom.name OW')
    
   In [6]: rdf_xx,rdf_yy=ion.rdf(setA=list1,setB=list2,bins=1500,segment=[0.0,30.0])

.. sourcecode:: ipython

   In [2]: ion=molecule('run_ion.gro',coors=False,verbose=False) 
    
   In [3]: ion.load_traj('traj.dcd',frame='Next',verbose=False)
    
   In [4]: list1=ion.selection('atom.name OW')
    
   In [5]: list2=ion.selection('atom.resid.type Ion')
    
   In [6]: rdf_xx=pyn_math.binning(bins=1500,segment=[0.0,30.0])
    
   In [7]: rdf_yy=zeros(shape=(1500),dtype=float,order='Fortran')
    
   In [8]: num_frames=0
    
   In [9]: while ion.traj[0].io_opened:
      ...:     rdf_yy+=ion.rdf(setA=list1,setB=list2,traj=0,frame=0,bins=1500,segment=[0.0,30.0])
      ...:     num_frames+=1
      ...:     ion.traj[0].reload_frame()
   # End of file
    
   In [10]: rdf_yy=rdf_yy/(1.0*num_frames)







