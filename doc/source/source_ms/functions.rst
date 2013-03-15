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

.. seealso:: :ref:`ms-tut-selections`


Info
++++

.. method:: msystem.info()

   :return: None
   :prints out: General information of the molecular system.

Writting pdb files
++++++++++++++++++

xxx

------------------------

Trajectories
============

Load Trajectory
+++++++++++++++

.. method:: msystem.load_traj(file_input=None [,frame=None,begin=None,end=None,increment=1,units='frames',verbose=False])

   :arg str file_input: Name of trajectory file in the following formats: pdb, gro, xtc, trr, dcd.
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
   :prints out: General information of the whole set of trajectories loaded.

Delete Trajectory
+++++++++++++++++

.. method:: msystem.delete_traj(index='ALL')

   :arg index: index of trajectory to be deleted in the list msystem.traj
   :arg index: int or 'ALL'
   :return: None


Load Frames
+++++++++++

To append a new frame to the list traj.frame:

.. method:: traj.upload_frame(frame='next'[,begin=None,end=None,increment=1,units=None])

   :arg frame: Index of frames to be loaded.
   :type frame: int or list or 'next' or 'all'.
   :arg begin: Frame to start the loading process.
   :type begin: int.
   :arg end:  Frame to finish the loading process.
   :type end: int.
   :arg increment: Interval to load frames.
   :type increment: int.
   :arg units: Units for "begin", "end" and "increment". [only units='frames' available]
   :type units: str.
   :return: None

To rewrite an old frame of the traj.frame list with a new frame:

.. method:: traj.reload_frame(frame='next',old=0)

   :arg frame: Index of frames to be loaded.
   :type frame: int or 'next'.
   :arg int old: Index of frame in the list traj.frame to be rewritten.

Info Frames
+++++++++++

.. method:: traj.info()

   :return: None
   :prints out: General information of the whole set of the trajectory.


Delete Frames
+++++++++++++

.. method:: traj.delete_frame(frame='ALL'[,begin=None,end=None,increment=1,units=None])

   :arg frame: Index of frame in the list traj.frame to be deleted.
   :type frame: int or list or 'ALL'.
   :arg begin: Frame to start the deleting process. [not available yet]
   :type begin: int.
   :arg end:  Frame to finish the deleting process. [not available yet]
   :type end: int.
   :arg increment: Interval to delete frames. [not available yet]
   :type increment: int.
   :arg units: Units for "begin", "end" and "increment". [only units='frames' available]
   :type units: str.
   :return: None

.. seealso:: .


Writting dcd files
++++++++++++++++++

.. method:: traj.write(file_name=None,frame='ALL',begin=None,end=None,increment=1,units=None,action=None)

   :arg str file_name: Name of new trajectory file.
   :arg frame: Indexes of frames to be written in the new trajectory file.
   :type frame: int or list or 'ALL'
   :arg begin: Frame to start the deleting process. [not available yet]
   :type begin: int.
   :arg end:  Frame to finish the deleting process. [not available yet]
   :type end: int.
   :arg increment: Interval to delete frames. [not available yet]
   :type increment: int.
   :arg units: Units for "begin", "end" and "increment". [only units='frames' available]
   :type units: str.
   :arg action: "OPEN" or "CLOSE" to open or close a file without writting, None to writte frames in the opened new trajectory file.
   :type action: str.
   :return: None

.. seealso:: :ref:`ms-tut-convert-traj`


------------------------

Analysis
========

Ramachandran Map 
++++++++++++++++ 

The function computes the pairs of angles phi-psi for any list of
residues and frames.  Since each dihedral angle is computed with atoms
from 2 different residues, the relationship between angle and
residue is given by the atom 'CA'.

.. method:: msystem.ramachandran_map(resid='ALL',traj=0,frame='ALL',pdb_index=False,legend=False)

   :arg resid: List of residue indexes or pdb indexes.
   :type resid: int, list[int] or 'ALL'
   :arg int traj: Index of trajectory to be analysed.
   :arg frame: List of frame indexes.
   :type frame: int, list[int] or 'ALL'
   :arg bool pdb_index: Residues in input and output are identified by the pdb indexes if True.
   :arg bool legend: Key of angles in output if True.
   :returns: **angles**, **keys** [if legend]; List of pairs phi-psi for the corresponding residues and frames, and key legend if choosen. 
   :rtype: 
   	 * angles: 
	   	   * numpy.array[num_frames,num_resids,2].
   	           * numpy.array[num_resids,2] (if num_frames=1)
		   * numpy.array[num_frames,2] (if num_resids=1)
	           * numpy.array[2]	       (if num_frames=num_resids=1)
	 * keys: 
	   	   * list[num_resids][2]
		   * list[num_resids][2] (if num_resids=1)


.. seealso:: :ref:`ms-tut-rama-map`

.. note:: Depending on how this method is used, it can result with a
   low performance. Check :meth:`msystem.dihedral_angle` for a better performance.


Covalent chains
+++++++++++++++

.. method:: msystem.selection_covalent_chains(chain=None,select='protein')

   :arg chain: List of atom names to find as covalently bonded chains in "select".
   :type chain: list[str]
   :arg select: Selection or set of atoms where the method looks for covalent chains. (see: :meth:`msystem.selection`)
   :type select: str
   :returns: **chains**; List of chains found. Each chain is a list of the corresponding atom indexes.
   :rtype: 
   	   * chains:
			* list[num_chains][atom_indexes]

.. seealso:: :ref:`ms-tut-any-dihang`


Dihedral Angles
+++++++++++++++

The method computes the value of any list of 4 atoms defining a
dihedral angle. The output

.. method:: msystem.dihedral_angle(covalent_chain=None,traj=0,frame='ALL')

   :arg covalent_chain: List of 4 atoms defining the dihedral angle, or list of lists. (see: :meth:`msystem.selection_covalent_chains`)
   :type covalent_chain: list[int] or list[list[int]]
   :arg int traj: Index of traj to be analysed.
   :arg frame: List of frame indexes.
   :type frame: int, list[int] or 'ALL'
   :returns: **angles**; List of angle values (radians) in the same order found in input parameter covalent_chain.
   :rtype: 
   	   * angles:
			* numpy.array[num_frames,num_angs].
   	     		* numpy.array[num_angs] (if num_frames=1)
   	     		* numpy.array[num_frames] (if num_angs=1)

.. seealso:: :ref:`ms-tut-any-dihang`
