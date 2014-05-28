from liboptim import glob as fs

def load_net(net):
    fs.load_net(net.T_start,net.T_ind,net.T_wl,net.num_nodes,net.k_total)

    pass

def random_pos():
    fs.random_distribution(1,3)

    pass


