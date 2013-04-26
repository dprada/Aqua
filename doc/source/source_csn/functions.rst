Functions
*********

Topology
========

Loading network
+++++++++++++++

.. method:: network.load_net(name_file,format='text',verbose=True)

   :arg str name_file: Name of input file with network
   :arg str format: Format of input file: 'text' or 'native'
   :verbose: Prints out general information of the network loaded.
   :return: None
   :prints out: general information if verbose=True.

Loading labels
++++++++++++++

.. method:: network.load_labels(name_file,format='text')

   :arg str name_file: Name of input file with network
   :arg str format: Format of input file: 'text' or 'water'
   :return: None
   :prints out: None

Adding node
+++++++++++

.. method:: network.add_node(label=None, weight=0,iout=False)

   :arg label: Label of new node.
   :type label: str or int.
   :arg float weight: Weight of new node.
   :arg bool iout: Returning index of new node in network.
   :return: Index of new node if iout=True.
   :prints out: None

Adding link
+++++++++++

.. method:: network.add_link(node_origin,node_final,weight=0,index_origin=False,index_final=False,iout=False)

   :arg node_origin: source of the link.
   :type node_origin: int or str.
   :arg node_final: destination of the link.
   :type node_final: int or str.
   :arg float weight: weight of link.
   :arg bool index_origin: True if node_origin takes the integer value of the index in the network.
   :arg bool index_final: True if node_final takes the integer value of the index in the network.
   :arg bool iout: Returning the indexes of nodes linked.
   :return: Index of linked nodes if iout=True.
   :prints out: None

Info
++++

.. method:: network.info(update=True,verbose=True)

   :arg bool update: Computing and updating general variables as network.weight, network.num_nodes, network.k ...
   :arg bool verbose: Prints out general information.
   :return: None
   :prints out: General information of the network if verbose.


Merging network
+++++++++++++++

.. method:: network.merge_net(net=None,traj=None,verbose=True)

   :arg net: network to be merged or included.
   :type net: network class.
   :arg traj: None.
   :arg bool verbose:  Prints out general information.
   :prints out: General information of the network if verbose.

Extracting network
++++++++++++++++++

.. method:: network.extract_net(nodes=None,verbose=True)

   :arg nodes: List of nodes with which a new network is built.
   :type nodes: list[int]
   :arg bool verbose:  Prints out general information.
   :return: network class.
   :prints out: General information of the new network if verbose.


Writting network
++++++++++++++++

.. method:: network.write_net(name_file=None,format='text',pymol=False,with_index=True,with_cluster=False)

   :arg str name_file: Name of new file.
   :arg str format: Format of output file: 'text' or 'native'.
   :arg bool pymol: True if the file will be read by pymol.
   :arg bool with_index: The network will be written with node indexes if True.
   :arg bool with_cluster: The hnodes will be written with a cluster indentification.

Writting labels
+++++++++++++++

.. method:: network.write_labels(name_file=None,format='text')

   :arg str name_file: Name of new file
   :arg str format: Format of output file: 'text'


Clustering
==========

MCL
+++

Gradient Clusters
+++++++++++++++++

Free Energy Profiles
====================

CFEPs
+++++

.. method:: network.cfep(mode='pfold',A=0,B=0,num_bins=0,num_iter=20000,KbT=((0.0019872*300.0)))

   :arg str mode: Version of CFEP: 'pfold' or 'mfpt'.  
   :arg int A: Index of node with Pfold=1.0 if mode='pfold', or with mfpt=0.0 if
               mode='mfpt'.  
   :arg int B: Index of node with Pfold=0.0 if
               mode='pfold'. Not required if mode='mfpt'.  
   :arg int num_bins: Number of bins to discretize the X coord (Za/Z). If num_bins=0, a
                       point [X,Y] is computed every ranked node according to the Pfold or mfpt value.
   :arg int num_iter: Number of iterations to compute the Pfold or mfpt values for each node.
   :arg float KbT: Value of KbT. (Kb=0.0019872 kcal/mol/K or Kb=0.0083144 kJ/mol/K)
   :returns: 
             1. xcfep: X coord of CFEP (Za/Z).
	     2. ycfep: Y coord of CFEP (-KbT*log(Zab/Z))
	     3. zcfep: pfold or mfpt corresponding to X.
	     4. node_vals: values of pfold or mfpt per node.
	     5. node_x: location in xcfep array per node.
   :rtype:
	  1. xcfep: numpy.array[num_bins](float) or numpy.array[network.num_nodes](float) if num_bins=0
	  2. ycfep: numpy.array[num_bins](float) or numpy.array[network.num_nodes](float) if num_bins=0
	  3. zcfep: numpy.array[num_bins](float) or numpy.array[network.num_nodes](float) if num_bins=0
	  4. node_vals: numpy.array[network.num_nodes](float)
	  5. node_x: numpy.array[networ.num_nodes](int)
   :prints out: None
   :example: xcfep, ycfep, zcfep, node_vals, node_x = network.cfep(mode='pfold', A=10, B=30)
