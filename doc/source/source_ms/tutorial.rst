
Tutorials
*********

First of all, lets load Aqua in our script or in a ipython session:

.. sourcecode:: ipython

     In [1]: from aqua import *

Some basic notions on python will be assumed along this tutorial. If
you just landed here without any idea on python, have a look to a
Python manual or tutorial before.

Some files are required to follow the tutorials:

- Met-Enkephalin (`PDB: 1PLW
  <http://www.rcsb.org/pdb/explore.do?structureId=1PLW>`_ ):
  :download:`metenk.pdb <tutorial/metenk.pdb>` and
  :download:`traj_metenk.xtc <tutorial/traj_metenk.xtc>`

- Catabolite Activator Protein (`PDB: 2WC2
  <http://www.rcsb.org/pdb/explore.do?structureId=2WC2>`_ ):
  :download:`2WC2.pdb <tutorial/metenk.pdb>`.

----------------------
 

Loading/Writting the topology
=============================

Loading
+++++++

A system can be downloaded from the Protein Data Bank.

.. sourcecode:: ipython

   In [2]: mol_test=msystem(download='2WC2',verbose=True)
   # File saved as 2WC2.pdb
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions
   # 20  frames/models in traj 0


Or 

.. sourcecode:: ipython

   In [3]: mol2_test=msystem('2WC2.pdb')

   In [4]: mol2_test.info()
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
without coordinates can be created adding the option coors=False:

.. sourcecode:: ipython

   In [5]: mol3_test=msystem('2WC2.pdb',coors=False)
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions

   In [6]: mol3_test.info_trajs()
   # No coordinates

   In [7]: mol2_test.info_trajs()
   # 20 frames/models in traj 0


.. Note:: For further information on the functions and attributes used
   visit: :class:`msystem`, :meth:`msystem.info`, :meth:`msystem.info_trajs`.



Navigating
++++++++++



Writting
++++++++
XXX

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
using the files :download:`GSGS.pdb <../../tutorials/systems_tut1/GSGS.pdb>`
and :download:`GSGS.dcd <../../tutorials/systems_tut1/GSGS.dcd>`.

.. sourcecode:: ipython

   In [2]: GSGS=msystem('GSGS.pdb')
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

   In [5]: GSGS.msystem('GSGS.dcd','ALL')
   # 10 frames/models loaded.

.. sourcecode:: ipython

   In [2]: GSGS=msystem('GSGS.pdb',coors=False,verbose=False)
    
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

   In [2]: GSGS=msystem('GSGS.pdb',coors=False,verbose=False)

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

   In [2]: GSGS=msystem('GSGS.pdb',coors=False,verbose=False)

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

.. _ms-tut-convert-traj:

Converting a trajectory into other format
+++++++++++++++++++++++++++++++++++++++++

Right now the output format is only dcd.

.. sourcecode:: ipython

   In [2]: metenk=msystem('metenk.pdb')
    
   In [3]: metenk.load_traj('traj_metenk.xtc')
   
   In [4]: metenk.info_trajs()
   # 0 frames/models in traj 0

   In [5]: metenk.traj[0].write('traj_metenk.dcd',action='Open')
    
   In [6]: while 1:
      ...:     metenk.traj[0].reload_frame()
      ...:     if metenk.traj[0].io_end:
      ...:         metenk.traj[0].write(action='Close')
      ...:         break
      ...:     metenk.traj[0].write(frame=0)

.. seealso:: :class:`msystem`, :class:`traj`, :meth:`msystem.load_traj`, :meth:`msystem.info_trajs`, :meth:`traj.reload_frame`, :meth:`traj.write`

.. _ms-tut-selections:

How to make atoms selections
============================

The syntax is close to the aqua syntax.
There are few special key words.

.. sourcecode:: ipython

   In [2]: metenk=msystem('metenk.pdb')
    
   In [3]: list1=metenk.selection('backbone')

   In [4]: list2=metenk.selection('atom.name N CA C O')
    
   In [5]: print list1; print list2
   [0, 4, 21, 22, 23, 25, 28, 29, 30, 32, 35, 36, 37, 39, 55, 56, 57, 59, 72]
   [0, 4, 21, 22, 23, 25, 28, 29, 30, 32, 35, 36, 37, 39, 55, 56, 57, 59, 72]

   In [5]: metenk.selection('resid.name PHE and backbone')
   Out[5]: [37, 39, 55, 56]

   In [6]: metenk.selection('chain.name A and atom.donor')
   Out[6]: [0, 23, 30, 37, 57]

   In [7]: metenk.selection('(atom.resid.name GLY and not atom.name N CA C O H) or (atom.name O1)')  
   Out[7]: [26, 27, 33, 34, 73]

We can also make use of the expression 'within X of', X is a float number indicating a distance threshold.

.. sourcecode:: ipython

   In [8]: metenk.load_traj('traj_metenk.xtc',frame='ALL',verbose=False)

   In [9]: list0=metenk.selection('atom.name OW within 3.0 of resid.type Protein')

   In [10]: for ii in range(metenk.traj[0].num_frames):
      ....:     print len(list0[ii]), 'waters below 3.0 in frame', ii
      ....: 
   25 waters below 3.0 in frame 0
   30 waters below 3.0 in frame 1
   28 waters below 3.0 in frame 2
   ...

   In [11]: list0=metenk.selection('(atom.name CA and not resid.name TYR) \
      ....:    within 6.0 of (atom.name CA and resid.name TYR)')

   In [12]: for ii in range(metenk.traj[0].num_frames):
      ....:    print 'At time:', metenk.traj[0].frame[ii].time
      ....: for jj in list0[ii]:
      ....:     print '   ',metenk.atom[jj].info(),'within 6.0 of',metenk.atom[4].info()
   At time: 0.0
      CA-26/GLY-2 within 6.0 of CA-5/TYR-1
   At time: 10.0
      CA-26/GLY-2 within 6.0 of CA-5/TYR-1
   At time: 20.0
      CA-26/GLY-2 within 6.0 of CA-5/TYR-1
      CA-33/GLY-3 within 6.0 of CA-5/TYR-1
   ...


.. seealso:: :meth:`msystem.selection`



Distances
=========

.. _ms-tut-dists:

Computing distances
+++++++++++++++++++

The distance between a set of n1 atoms, sel1, and a set of n2 atoms,
sel2.  If only one frame is analysed the output is a 2D matrix
[n1,n2].  This way the distance between the i-th atom in sel1 and
j-th in sel2 correspond to the output element [i,j].

If more than one frame is analysed the output is a 3D matrix
[n1,n2,number_frames] (following the previous notation).

The method by construction is faster if n1<n2.

.. sourcecode:: ipython

   In [2]: metenk=msystem('metenk.pdb')

   In [3]: metenk.load_traj('traj_metenk.xtc',frame='ALL')

   In [4]: dists,keys1,keys2=metenk.distance('resid.pdb_index 1','resid.pdb_index 4',legend=True)
   
   In [5]: for ii in range(metenk.traj[0].num_frames):
      ...:         print 'The distance between atoms index',keys1[10],\
      ...:               'and',keys2[5],'is',dists[ii,10,5],'in frame',ii
      ...: 
   The distance between atoms index 10 and 42 is 7.30732658606 in frame 0
   The distance between atoms index 10 and 42 is 7.54015357744 in frame 1
   The distance between atoms index 10 and 42 is 7.22121193005 in frame 2
   ...

If only a sets of atoms is provided, the distance among them is computed:

.. sourcecode:: ipython

   In [6]: CAs=metenk.selection('atom.name CA')

   In [7]: dists=metenk.distance(sel1=CAs,frame=10)

   In [8]: print metenk.atom[CAs[2]].info(),'--',metenk.atom[CAs[4]].info(), dists[2,4]
   CA-33/GLY-3 -- CA-60/MET-5 5.81840212291

.. seealso:: :class:`msystem`, :class:`traj`, :class:`atom`, :meth:`msystem.load_traj`, :meth:`msystem.info_trajs`, :meth:`msystem.selection`, :meth:`msystem.distance`


Neighbors and Ranked contacts
+++++++++++++++++++++++++++++

The function neighbs can help with its different options to approach this problems.
Notice that the cut-off here is the limit in the ranked list of closest neighbors or a given distance.
In the former case the output can be sorted or not by distance. 
If only contact map is need, maybe the following section is suitable because of its efficience.

Neighbors with a distance lower or equal than a certain threshold:



Contact Maps
============





Dihedral Angles
===============

.. _ms-tut-rama-map:

Ramachandran Map
++++++++++++++++

.. sourcecode:: ipython

   In [2]: metenk=msystem('metenk.pdb',coors=True)

   In [3]: angs,keys=metenk.ramachandran_map(legend=True)

   In [4]: print keys[3][0],'=',angs[3][0],'  ', keys[3][1],'=',angs[3][1]
   Phi 3 = -1.92163955875    Psi 3 = -0.642196409677

   In [5]: print metenk.ramachandran_map(resid=3)
   [-1.92163956 -0.64219641]

   In [6]: dangs,dkeys=metenk.ramachandran_map(resid=4,pdb_index=True,legend=True)

   In [7]: print keys[0],'=',angs[0],'  ', keys[1],'=',angs[1]
   Phi 4 = -1.92163955875    Psi 4 = -0.642196409677

In case of having different frames to analyse coming from a trajectory:

.. sourcecode:: ipython

   In [8]: metenk.info_trajs()
   # 1 frames/models in traj 0

   In [9]: print metenk.traj[0].name
   metenk.pdb

   In [10]: metenk.delete_traj()

   In [11]: metenk.info_trajs()
   # No coordinates

   In [12]: metenk.load_traj('traj_metenk.xtc',frame='ALL',verbose=True)
   # 21 frames/models in traj 0

   In [13]: angs,keys=metenk.ramachandran_map(resid=[3,5],frame=[0,10,20],pdb_index=True,legend=True)

   In [14]: print keys[0][1],'(step=[0,10,20]) =',angs[:,0,1]
   Psi 3 (step=[0,10,20]) = [ 0.87594459  0.64741795  1.22817325]

   In [15]: angs=metenk.ramachandran_map(resid=2,frame='ALL')

   In [16]: for ii in range(0,metenk.traj[0].num_frames,10):
      ....:    print metenk.resid[2].info(),'| t=', metenk.traj[0].frame[ii].time, \
      ....:          '| [phi,psi]=', angs[ii,:]
   GLY-3 | t= 0.0 | [phi,psi]= [ 2.01215258  0.87594459]
   GLY-3 | t= 100.0 | [phi,psi]= [ 1.90493036  0.64741795]
   GLY-3 | t= 200.0 | [phi,psi]= [ 1.9432302   1.22817325]


.. seealso:: :class:`msystem`, :class:`traj`, :class:`resid`, :meth:`msystem.load_traj`, :meth:`msystem.info_trajs`, :meth:`msystem.delete_traj`, :meth:`msystem.ramachandran_map`

.. note:: Depending on how this method is used, it can result with a
   low performance. Check the next section :ref:`ms-tut-any-dihang` for a better performance.


.. _ms-tut-any-dihang:

Computing any dihedral angle
++++++++++++++++++++++++++++

A dihedral angle can be computed with a set of 4 atoms. This way, the
computation of dihedral angles goes beyond the analysis done with
:meth:`msystem.ramachandran_map`.

There is a method to extract sets of 4 atoms corresponding to
consecutive covalently bound chains of atoms:

.. sourcecode:: ipython

   In [2]: metenk=msystem('metenk.pdb')

   In [3]: phi_3=metenk.selection_covalent_chains(chain=['C','N','CA','C'],select='resid.index 2 3')

   In [4]: psi_3=metenk.selection_covalent_chains(chain=['N','CA','C','N'],select='resid.index 3 4')

   In [5]: xi1_3=metenk.selection_covalent_chains(chain=['N','CA','CB','CG'],select='resid.index 3')

   In [5]: omegas=metenk.selection_covalent_chains(chain=[['CA','CH3'],'C','N',['CA','CH3']])

   In [6]: print omegas
   [[4, 21, 23, 25], [25, 28, 30, 32], [32, 35, 37, 39], [39, 55, 57, 59]]

With this sets of 4-tuples the dihedral angles can be computed as:

.. sourcecode:: ipython

   In [7]: metenk.load_traj('traj_metenk.xtc',frame='ALL')

   In [8]: angs_phi_3=metenk.dihedral_angle(phi_3)
   
   In [9]: for ii in range(metenk.traj[0].num_frames):
     ....:     print 'Phi 3 (step='+str(ii)+')', angs_phi_3[ii]
     ....: 
   Phi 3 (step=0) -1.92022469777
   Phi 3 (step=1) -1.59101737012
   Phi 3 (step=2) -1.61015325568
   Phi 3 (step=3) -1.99394603011
   ...

   In [10]: print metenk.dihedral_angle(covalent_chain=psi_3,frame=10)
   -0.23084216

   In [11]: ang_omegas=metenk.dihedral_angle(covalent_chain=omegas,frame=[0,10,20])

   In [12]: for ii in range(len(omegas)):
      ....:     print 'omega_'+str(ii+1)+' (step=[0,10,20])= ',str(ang_omegas[:,ii])
      ....: 
   omega_1 (step=[0,10,20])=  [ 3.13688483 -3.08321194 -3.07475646]
   omega_2 (step=[0,10,20])=  [ 3.1379404   3.08377778 -3.11113559]
   omega_3 (step=[0,10,20])=  [ 3.13383903 -3.11193193 -3.0400114 ]
   omega_4 (step=[0,10,20])=  [ 3.13693669 -3.0892444  -3.03312968]

.. seealso:: :class:`msystem`, :class:`traj`, :meth:`msystem.load_traj`, :meth:`msystem.selection_covalent_chains`, :meth:`msystem.dihedral_angle`

.. note:: If the angles phi and psi are the goal of the analysis, the method :meth:`msystem.ramachandran_map` can be another tool to be considered.




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


