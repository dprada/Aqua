
Networks
*********

Description of the network object and attributes.
This is a test to include the comments from the code.

nada

Network
=======

.. class:: network(file_net=None,file_labels=None,net_format='text',labels_format='text',directed=True,kinetic=False,verbose=True)

   :arg str file_net: Input file with network (see formats)
   :arg str file_labels: Input file with labels for nodes (see formats)
   :arg str net_format: Format of input file: 'text',... (see formats)
   :arg str labels_format: Format of labels file: 'text',... (see formats)
   :arg bool directed: Network with directed links.
   :arg bool kinetic: Kinetic network. (see)
   :arg bool verbose: Printing out information about the loaded network.
   :return: returns the object intialized.

   :Example: 

	    (followed by a blank line)


.. attribute:: network.node

   List of nodes. Each node is an object of node class. (list)

.. attribute:: network.num_nodes

   Total number of nodes in the network. (int)

.. attribute:: network.cluster

   List of clusters. Each cluster is an object of cluster class. (list)

.. attribute:: network.num_clusters

   Total number of clusters detected in the network

.. attribute:: network.component

   List of components. Each component is an object of cluster class. (list)

.. attribute:: network.num_components

   Total number of components in the network

.. attribute:: network.k_total

   Total number of links. Degree of the network. (int)

.. attribute:: network.k_max

   Maximum number of links of a node. Maximum degree of a node. (int)

.. attribute:: network.k_out

   Total number of links out of nodes if the network is directed. (int)


.. attribute:: network.k_in

   Total number of links in of nodes if the network is directed. (int)

.. attribute:: network.weight

   Total weight of nodes if they are weighted. (int)

.. attribute:: network.labels

   Dictionary of with labels and indexes of nodes. (dict[str]=int)

.. attribute:: network.directed

   True if the links are directed. (bool)

.. attribute:: network.kinetic 

   True if it has the properties of a Conformational Space Network, or kinetic network.

.. attribute:: network.file_net

   Name of input file with the network. (str)

.. attribute:: network.file_labels

   Name of input file with labels network. (str)




.. seealso:: blabla

Node
++++

Nodes are represented by objects with the common attributes of class node:

.. attribute:: node.label

   Label of node. (str)

.. attribute:: node.link

   Dictionary of links with index of the destination node as key, and weight of link as value. (dict[int]=float)

.. attribute:: node.k

   Degree or number of links of node. (int)

.. attribute:: node.k_out

   Degree or number of links out of this node. If the network is directed. (int)

.. attribute:: node.k_in

   Degree or number of links in this node. If the network is directed. (int)

.. attribute:: node.weight

   Weight of node. If the network is weighted. (float)

.. attribute:: node.cluster

   Index of cluster to which the node belongs.

.. attribute:: node.component

   Index of component to which the node belongs.

.. attribute:: node.coors

   List of coordinates of this node for a graphical representation. (list[float])

.. attribute:: node.color

   Color of node in a graphical representation. (str,int,float...)

.. attribute:: node.size

   Size of network in a graphical representation. (int, float)


Cluster
+++++++

Clusters are represented by objects with the common attributes of class cluster:

.. attribute:: cluster.num_nodes

   Number of nodes in the cluster. (int)

.. attribute:: cluster.nodes

   List of indexes of nodes which belong to the cluster. (list[int])

.. attribute:: cluster.weight

   Total weight of cluster if the nodes are weighted. (int)

.. attribute:: cluster.label

   Label of the cluster. Label of the most weighted node in the cluster -if the network is weighted-. (str)

.. attribute:: cluster.link

   Dictionary of links with the index of the destination cluster as key, and its weight as value. (dict[int]=float)

.. attribute:: cluster.k

    Degree or number of links of the cluster. (int)

.. attribute:: cluster.k_out

   Degree or number of links out of the cluster if the network is directed. (int)

.. attribute:: cluster.k_in

   Degree or number of links in the cluster if the network is directed. (int)

Component
+++++++++

Components are represented by objects with the common attributes of a cluster:

.. attribute:: component.num_nodes

   Number of nodes in the component. (int)

.. attribute:: component.nodes

   List of indexes of nodes which belong to the component. (list[int])

.. attribute:: component.weight

   Total weight of component if the nodes are weighted. (int)

.. attribute:: component.label

   Label of the component. Label of the most weighted node in the component -if the network is weighted-. (str)




