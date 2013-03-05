Molecular Systems
*****************

Info about Molecular Systems
============================


Units
+++++

Length: Angstroms.
Angles: Degrees.

Periodic Box and Cell
+++++++++++++++++++++

In frame, class trajectory, box is the orthorombic set of vectors and
cell (unit cell) has the length of three vectors vi in the elements cell[i,i] and
the angles alpha, beta and gamma in cell[1,2], cell[1,3] and cell[2,3]
respectively.

See: http://www.mail-archive.com/gmx-users@gromacs.org/msg28019.html
See also: SUBROUTINE TRICLINIC (cell,box) in libdcdfile.f90 

Pynoramix sets the box always in the positive cuadrant. This way any
atom is always in {[0,Lx),[0,Ly),[0,Lz)}. It implies that molecules
can be split if they are over the cell edges.

----------------------


Molecule class
==============

Subclasses
++++++++++

Functions
+++++++++

Trajectory class
================

Functions
+++++++++

Functions
=========

Selection
+++++++++



 


