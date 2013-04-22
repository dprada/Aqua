Molecular Systems
+++++++++++++++++

The topology of my molecular system is not recognized by AquaLab.
=================================================================

When a :class:`msystem` is initialized with a pdb or gro file, the topology
-covalent bonds- and the properties of atoms are taken from the data
stored in the directory Aqua/top_par. These pre-defined topologies have
been implemented according to the residue and atom names found in the
most popular water and protein models included in GROMACS, CHARMM or NAMD.

If your system is not described by AquaLab there are three ways to
proceed:

1. Create your own psf file to initialize the :class:`msystem`. If the
   psf file does not include a list of donors and acceptors, they have to
   be set up (see...).

2. Find a template to create a new topology in
   :download:`Aqua/top_par/top_par_template.py<static_ms/top_par_template.py>`. Include
   an edited copy with the description of your molecular system in the
   directory Aqua/top_par. And add to the file Aqua/top_par/top_par.py
   the line: "user_topol.append('your_top_par_file_name')". It is time
   now to run your script calling :class:`msystem`.

3. Find a template to create a new topology in the file
   :download:`Aqua/top_par/top_par_template.py<static_ms/top_par_template.py>`. Include
   in your working directory an edited copy with the description of
   your molecular system. Run your script loading the new topology
   with xxx before initializing the :class:`msystem`.




