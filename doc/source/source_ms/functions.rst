Functions
*********

Topology
========

Loading
+++++++

This function should be avoided by the user if there is not a clear
reason. When a new object of class msystem is initialized this
function is called implicitly by construction.

.. method:: msystem.load_topol(name_file)

   :arg str name_file: Name of input file with topology
   :return: None

See topology formats for additional info.


Selection
+++++++++

.. method:: msystem.selection(condition=None [,traj=0,frame='ALL',pbc=True])

   :arg str condition: Logical sentence with selection.
   :arg int traj: Index of trajectory to analyse if a condition over coordinates needs to be evaluated.
   :arg frame: Index of frame to analyse if a condition over coordinates needs to be evaluated.
   :type frame: list, int or str
   :arg bool pbc: Periodic Boundary Conditions if a condition over coordinates needs to be evaluated.
   :return: list of atom indexes.
   :rtype: list of integers.

Syntaxis
--------

The syntaxis used by the function msystem.selection() is given by the
classes atom, resid and chain.  In addition, boolean operations can be
done using the words: "AND", "OR", "NOT", "IN", "WITHIN", "OF"; no
matter if are written in capital or lower case, and the characters:
"[", "]", "(", ")" and ",".

Some special words can be used although they have an easy translation:

     - "backbone":  '(atom.name N CA C O)'
     - "sidechain": '(resid.type Protein and not atom.name N CA C O H1 H2)'
     - "protein":   '(resid.type Protein)'
     - "water":     '(resid.type Water)'

Find some examples in section Tutorial-Syntaxis.


Info
++++

.. method:: msystem.info()

   :return: None
   :prints out: General information of the molecular system.

Writting pdb files
++++++++++++++++++


Trajectories
============

Load Trajectory
+++++++++++++++

.. method:: msystem.load_traj(file_input=None,[frame=None,begin=None,end=None,increment=1,units='frames',]verbose=False)

   :arg str file_input: Name of trajectory file (pdb, gro, xtc, trr, dcd) 
   :arg frame: Index of frames to be loaded.
   :type frame: int or list.
   :arg begin: Frame to start the loading process.
   :type begin: int.
   :arg end:  Frame to finish the loading process.
   :type end: int.
   :arg increment: Interval to load frames.
   :type increment: int.
   :arg units: Units for "begin", "end" and "increment". [only units='frames' available]
   :type units: str.
   :arg bool verbose: Prints out general information.
   :return: None
   :prints out: general information if verbose=True.


Info Trajectories
+++++++++++++++++

.. method:: msystem.info_trajs()

   :return: None
   :prints out: General information of the trajectories loaded.





Writting dcd files
++++++++++++++++++

Writting pdb files
++++++++++++++++++

