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
