Molecular system
****************

**Aqua** defines a moleculer system with the class *msystem*. This python class is the parent object where all info about the molecular system (topology, attributes, coordinates) as well as functions for its analysis is found.

.. class:: msystem(input_file=None,download=None,coors=False,with_bonds=True,missing_atoms=True,verbose=False)

   Supra class containing the whole data related with the system
   (topology, trajectory,...)  A file containing the description of
   the system must be provided to initialize this object, either a
   pdb, gro or psf file as input_file, either a pdb code to be
   downloaded from the `Protein Data Bank <http://www.rcsb.org>`_.

   :arg input_file: Input file in the following formats: pdb, gro, psf.
   :type input_file: String
   :arg download: PDB code to be downloaded.
   :type download: String
   :arg coors: Loading coordinates from the topology file (if there exist).
   :type coors: Bool
   :arg with_bonds: Automatic covalent bonds detection.
   :type with_bonds: Bool
   :arg missing_atoms: Printing out missing atoms.
   :type missing_atoms: Bool
   :arg verbose: Printing out information about the loaded system.
   :type verbose: Bool

.. attribute:: msystem.file_topol
   
   Name of the input file with topology.  (String)

.. attribute:: msystem.file_topol_type 

   Extension of the topology file: pdb, gro or psf. (String)

.. attribute:: msystem.atom

   List of atoms. Each atom is an object of unit class. (List)

.. attribute:: msystem.resid

   List of residues. Each residue is an object of residue class. (List)

.. attribute:: msystem.chain

   List of chains. Each chain is an object of xx class. (List)

.. attribute:: msystem.chains

   List of chain names. (List)

.. attribute:: msystem.ion 

   List of ions. Each ion is an object of xx class. (List)

.. attribute:: msystem.water

   List of waters. Each water is an object of xx class.

.. attribute:: msystem.water_model

   Water model found in the system. (String)

.. attribute:: msystem.acceptors

   List of indexes of acceptor atoms. (List)

.. attribute:: msystem.donors

   List of pairs: [index of donor atom, index of H atom bond to donor]. (list of lists)

.. attribute:: msystem.num_atoms

   Total number of atoms. (int)

.. attribute:: msystem.dimensionality

   Dimensionality of the molecular system: 3*num_atoms (int)

.. attribute:: msystem.name

   Name of the system. Taken from file_topol. (String)

.. attribute:: msystem.num_residues

   Total number of residues. (int)

.. attribute:: msystem.num_waters

   Total number of water molecules. (int)

.. attribute:: msystem.num_chains

   Total number of chains in the system. (int)

.. attribute:: msystem.num_ions

   Total number of ions in the system. (int)

.. attribute:: msystem.list_atoms

   List of atom indexes. If msystem was created from a file: range(msystem.num_atoms). (List)

.. attribute:: msystem.traj

   List of trajectories. Each trajectory is an object of xx class. (List)

.. attribute:: msystem.pdb_header

   Header of the a pdb file. (String)

.. attribute:: msystem.pdb_ss

   Secondary structure description from the pdb file. (String)


Atom
++++

Atoms are represented by objects of class atom with the common attributes:

.. attribute:: atom.name
   
   Name of atom as in msystem.input_file. (String)

.. attribute:: atom.type

   Type of atom. (String)

.. attribute:: atom.index

   Index of atom in msystem.atom[]. (int)

.. attribute:: atom.pdb_index

   Index of atom in msystem.input_file. (int)

.. attribute:: atom.resid

   Residue to which the atom belongs. (class residue)

.. attribute:: atom.chain

   Chain to which the atom belongs. (class chain)

.. attribute:: atom.covalent_bonds

   List of atom indexes with a covalent bond. (List)

.. attribute:: atom.acceptor

   The atom is acceptor for hydrogen bonds. (Bool)

.. attribute:: atom.donor

   The atom is donor for hydrogen bonds. (Bool)

.. attribute:: atom.mass

   Mass of atom. (Float)

.. attribute:: atom.charge

   Charge of atom. (Float)

.. attribute:: atom.vdw

   Van der Waals radius of atom. (Float)

.. attribute:: atom.occup

   Occupancy of atom in a PDB. (Float)

.. attribute:: atom.bfactor

   B-factor of atom in a PDB. (Float)

.. attribute:: atom.alt_log

   Alternate location in a PDB. 

.. attribute:: atom.code_ins_res

   Code of insertion of residues in a PDB.

.. attribute:: atom.seg_ident

   Index of segment in a PDB. (String)

.. attribute:: atom.elem_symb

   Element symbol in a PDB. (String)

.. attribute:: atom.type_pdb 

   Type of atom for the PDB: ATOM, HETATM,... (String)


Residue
+++++++

Residues are represented by objects of class atom with the common attributes:


.. attribute:: resid.name
   
   Name of residue as it comes from msystem.input_file. (String)

.. attribute:: resid.type

   Type of residue. (String)

.. attribute:: resid.index

   Index of residue in msystem.resid[]. (int)

.. attribute:: resid.pdb_index

   Index of residue in msystem.input_file. (int)

.. attribute:: resid.num_atoms

   Number of atoms in residue. (int)

.. attribute:: resid.list_atoms

   List of atom indexes in the residue. (List)

.. attribute:: resid.chain

   Chain to which the residue belongs. (class chain)

Chain
+++++

Chains are represented by objects of class chain with the common attributes:

.. attribute:: chain.name
   
   Name of chain as it comes from msystem.input_file. (String)

.. attribute:: chain.index

   Index of chain in msystem.chain[]. (int)

.. attribute:: chain.list_atoms

   List of atom indexes in the chain. (List)

Water
+++++

Water molecules are represented as residues in msystem.resid. In
addition there is an special class to include other common attributes
of these molecules.

.. attribute:: water.index

   Index of water molecule in msystem.resid list. (int)

.. attribute:: water.name

   Name of water molecule. Equal than resid.name. (String)

.. attribute:: water.model

   Name of water model. (String)

.. attribute:: water.list_atoms

   List of atom indexes in the water molecule. (List)

.. attribute:: water.O.index

   Index of atom Oxygen in msystem.atom list. (int)

.. attribute:: water.O.name

   Name of atom. Equal than atom.name (String)

.. attribute:: water.H1.index

   Index of atom Hydrogen 1 in msystem.atom list. (int)

.. attribute:: water.H1.name

   Name of atom. Equal than atom.name (String)

.. attribute:: water.H2.index

   Index of atom Hydrogen 2 in msystem.atom list. (int)

.. attribute:: water.H2.name

   Name of atom. Equal than atom.name (String)

.. attribute:: water.uvect_norm

   Perpendicular normal vector OH1xOH2. (Numpy array)







 


