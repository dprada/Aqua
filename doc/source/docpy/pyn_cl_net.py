class cl_node():
    """Fundamental unit to constitute a network together with links.

       :var label: label or key
       :type label: str
       :var weight: weight
       :type weight: float
       :var link: linked nodes, {node_index, weight_link}.
       :type link: dict [int, float]
       :ivar k_out: degree or connectivity with directed links out.
       :type k_out: int
       :cvar k_in: degree or connectivity with directed links in.
       :type k_in: int
       :var k: total degree or connectivity.
       :type k: int
       :var cluster: index of cluster the node belongs to
       :type cluster: int
       :var component: index of disconnected component the node belongs to
       :type component: int
       :var coors: spatial coordinates for representation 
       :type coors: N-array [float] 
    """

    def most_weighted_links(self,length=1):
        """Ranked indexes of connected nodes according to the weight of the links.
           
           :param length: N number of links requested
           :type length: int
           :return: ranked node indexes
           :rtype: N-dim list [int]
        """
        pass

class cl_cluster():

    def __init_(self):
        """Attributes of a set of nodes representing a cluster, community or disconnected component.
           :var nodes: indexes of nodes grouped in the set
           :type nodes: list [int]
           :var num_nodes: number of nodes grouped in the set
           :type num_nodes: int
           :var label: label or key
           :type label: str
           :var weight: total weight of the set
           :type weight: int 
           :var link: linked set, {set index, weight_link}
           :type link: dict [int, float]
           :var k_out: degree or connectivity with directed links out.
           :type k_out: int
           :var k_in: degree or connectivity with directed links in.
           :type k_in: int
           :var k: total degree or connectivity.
           :type k: int
"""

    pass

class net():
    """ Supra-structure composed by nodes"""
    def __init__(self,file_net=None,file_keys=None,directed=True,labels_format='2columns',verbose=True):
        """Attributes of a network 
           :var directed: network made by directed of links.
           :type directed: logic
           :var num_nodes: number of nodes in the network.
           :type num_nodes: int
           :var weight: total weight of the network.
           :type weight: int
           :var num_links: total number of links.
           :type num_links: int
           :var num_clusters: number of clusters in the network.
           :type num_clusters: int
           :var num_components: number of disconnected components in the network.
           :type num_components: int
           :var node: list of nodes
           :type node: list [class cl_node]
           :var cluster: list of cluster or communities.
           :type cluster: list [class cl_cluster]
           :var component: list of disconnected components.
           :type component: list [class cl_cluster]
           :var k_total: total number of links.
           :type k_total: int
           :var k_out: total number of links out.
           :type k_out: int
           :var k_in: total number of links in.
           :type k_in: int
           :var k_out_max: maximum number of links out of a node.
           :type k_out_max: int
           :var k_in_max: maximum number of links in of a node.
           :type k_in_max: int
           :var labels:  node index of each label, {label, index of node}
           :type labels: dict [str, int]
           :var clustering_method: name of clustering method if used
           :type clustering_method: str

           :var file_net: name of loaded file with the network.
           :type file_net: str
           :var file_keys: name of loaded file with labels.
           :type file_keys: str

           :var Ts: auxiliary arrays to optimize analysis.
           :type Ts: logic
           :var T_ind: auxiliary compact array with node indexes.
           :type T_ind: numpy array (format: Fortran).
           :var T_start: auxiliary array to handle Ts.
           :type T_start: numpy array (format: Fortran).
           :var T_weight: auxiliary compact array with link weights.
           :type T_weight: numpy array (format: Fortran).

           :param file_net: name of file to load a network. (see read_net)
           :type file_net: str, *opt*.
           :param file_keys: name of file to load a the labels. (see read_keys)
           :type file_keys: str, *opt*.
           :param labels_format: format of labels in file. (see read_keys)
           :type labels_format: str, ['2columns'].
           :param verbose: return network info. (see info)
           :type verbose: logic
        """

    def info(self):
        """ Print network general variables: num_nodes, k_total, weight."""    
        pass
        
    def add_node(self,new_node,weight=0):
        """Adding a new node to the network.

           :param new_node: key or label of the new node
           :type new_node: str,int,float...
           :param weight: weight of new node
           :type weight: int, *opt*
        """
        pass
