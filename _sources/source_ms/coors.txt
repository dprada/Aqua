Trajectories
************

Trajectory
++++++++++


.. class:: traj(file_input=None[,frame=None,begin=None,end=None,increment=1,units=None,verbose=True])

   All data related with a trajectory and the coordinates of a system
   is sort in this class. Although this class can be independently
   initalized, it is recommended the use of the method
   msystem.load_traj() for the same sake.

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

.. attribute:: traj.name

   Name of the input file with the trajectory. (String)

.. attribute:: traj.type

   Extension of the topology file: pdb, gro, xtc, trr, dcd. (String)

.. attribute:: traj.title

   Title or header of the file if there exist. (String)

.. attribute:: traj.num_frames

   Number of frames loaded from the input file and thereby saved in the list traj.frame. (int)

.. attribute:: traj.total_frames

   Total number of frames in the input file. (int)

.. attribute:: traj.precision

   Precision of the coordinates found in the input file. (float)

.. attribute:: traj.frame

   List of frames. Each fream is an object of frame class. (List)

.. attribute:: traj.io_opened

   Input/Output control variable: 1 -if the file is opened- or 0 -if it is closed-. (int)

.. attribute:: traj.io_end

   Input/Output control variable: 1 -if the file is at the end- or 0 -if it is not-. (int)

.. attribute:: traj.io_err

   Input/Output contro variable: 1 -if the i/o found and error- or 0 -if not-. (int)


Frame
+++++

Frames are represented by objects of class frame with common attributes and methods:

.. class:: frame()

.. attribute:: frame.time

   Time value of frame for this frame. (float)

.. attribute:: frame.step

   Integration step found in the trajectory for this frame.

.. attribute:: frame.precision

   Precision of coordinates as it was found in the trajectory.

.. attribute:: frame.lam

   --

.. attribute:: frame.model

   Index of model if the input file is a pdb with different models.

.. attribute:: frame.coors

   Array with the coordinates in the frame. (numpy.ndarray(shape=(num_atoms,3),dtype=float))

.. attribute:: frame.box

   Array with dimensions of the box in the frame. (numpy.ndarray(shape=(3,3),dtype=float))

.. attribute:: frame.cell

   Array with dimensions and angles of the unit chell for the frame. (numpy.ndarray(shape=(3,3),dtype=float))

.. attribute:: frame.orthogonal

   Control variable: 1 -if the box is orthogonal- or 0 -if not-. (int)

.. attribute:: frame.volume

   Volume of the box. (float)


Box and Cell
++++++++++++

In frame, class trajectory, box is the orthorombic set of vectors and
cell (unit cell) has the length of three vectors vi in the elements cell[i,i] and
the angles alpha, beta and gamma in cell[1,2], cell[1,3] and cell[2,3]
respectively.

See: http://www.mail-archive.com/gmx-users@gromacs.org/msg28019.html
See also: SUBROUTINE TRICLINIC (cell,box) in libdcdfile.f90 

Pynoramix sets the box always in the positive cuadrant. This way any
atom is always in {[0,Lx),[0,Ly),[0,Lz)}. It implies that molecules
can be split if they are over the cell edges.



 


