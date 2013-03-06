Introduction
************

A quick view
============
General description with a general list of what can be done and a first example on whats a molecule.


Example
+++++++


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



