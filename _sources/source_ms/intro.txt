Introduction
************

A quick view
============
General description with a general list of what can be done and a first example on whats a molecule.


Example
+++++++

Download the files :download:`metenk.pdb <tutorial/metenk.pdb>` and
:download:`traj_metenk.xtc <tutorial/traj_metenk.xtc>` to run a short
example with the Methionine-Enkephaline coming from the 1PLW pdb.

.. sourcecode:: ipython

     In [1]: from aqua import *

     In [2]: metenk=msystem('metenk.pdb',verbose=True)
     # System created from the file metenk.pdb :
     # 6015  atoms
     # 1490  residues
     # 2  chains
     # 1485  waters
     # 0  ions

     In [3]: idS=metenk.selection('atom.type S')

     In [4]: print 'The atom',metenk.atom[idS[0]].name+'-'+str(metenk.atom[idS[0]].pdb_index)+'/'+\
        ...: metenk.atom[idS[0]].resid.name+'-'+str(metenk.atom[idS[0]].resid.pdb_index),\
 	...: 'has the internal index:', idS[0]
     The atom SD-68/MET-5 has the internal index: 67

     In [5]: metenk.info_trajs()
     # No coordinates

     In [6]: metenk.load_traj('traj_metenk.xtc',frame='ALL',verbose=True)
     # 21 frames/models in traj 0

     In [7]: for ii in metenk.traj[0].frame:
        ...:     print '# Step',ii.step,'(t='+str(ii.time)+')','SD-68/MET-5 at', ii.coors[idS[0]]
	...: 
     Step 0 (t=0.0) SD-68/MET-5 at [ 19.8200016   12.90999985	 19.6400013 ]
     Step 5000 (t=10.0) SD-68/MET-5 at [ 19.87000084	13.27000046  19.44000053]
     Step 10000 (t=20.0) SD-68/MET-5 at [ 20.09000015	 12.96000004  18.93000031]
     ...

.. seealso:: :class:`msystem`, :class:`atom` , :meth:`msystem.selection`, :meth:`msystem.load_traj`, :meth:`msystem.info_trajs`, :class:`traj`, :class:`frame`



Units
=====

**Aqua** works with the following units no matter the input data.

- Length: Angstroms
- Angles: Degrees
- Time:   Picoseconds

Periodic Box and Cell
=====================

**Aqua** sets the coordinates of a molecular system always in the
  positive cuadrant. This way any atom is placed is found in
  {[0,Lx),[0,Ly),[0,Lz)}, where Lx, Ly and Lz are the dimensions of
  the cubic box. It implies that **molecules can be split if they are
  over the box edges**.

**Aqua** works with cubic boxes at the moment. Orthorombic cells with
  angles different from 90 degrees are not supported yet.



