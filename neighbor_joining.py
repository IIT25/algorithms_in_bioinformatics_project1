import numpy as np

class Node:
    def __init__(self, name, children = []):
        self.name = name
        self.children = children
        self.distance = None
        self.parent = None
    
    def write_in_newick(self, taxa):
        if self.name < len(taxa):
            name = taxa[self.name]
        else:
            name = ''
        if len(self.children) == 3:
            return f"({self.children[0].write_in_newick(taxa)},{self.children[1].write_in_newick(taxa)},{self.children[2].write_in_newick(taxa)}){name}" 
        elif len(self.children) == 2:
            return f"({self.children[0].write_in_newick(taxa)},{self.children[1].write_in_newick(taxa)}){name}:{round(self.distance,2)}"
        elif len(self.children) == 0:
            return f"{name}:{round(self.distance,2)}"


def read_distance_matrix(filename):
    with open(filename) as f:
        lines = f.readlines()
        n = int(lines[0])
        taxa = []
        distance_matrix = []
        for i in range(1, n+1):
            line = lines[i].split(' ')
            taxa.append(line[0])
            distance_matrix.append(list(map(float, line[1:])))
    return taxa, np.array(distance_matrix)

def join(nodes,i,j,distances,r,last_node):
    mask = [n for n in range(len(nodes)) if n not in [i,j]]

    new_node_name = last_node
    i_node = nodes[i]
    j_node = nodes[j]
    nodes = nodes[mask]

    d_ij = distances[i,j]
    distance_i = (d_ij+r[i]-r[j])/2
    distance_j = d_ij - distance_i
    i_node.distance = distance_i
    j_node.distance = distance_j
    
    nodes = np.append(nodes, Node(new_node_name, [i_node, j_node]))
    i_node.parent = nodes[-1]
    j_node.parent = nodes[-1]

    new_distances = np.array([(distances[i,m] + distances[j,m] - d_ij)/2 for m in range(len(distances))])
    new_distances = new_distances[mask]
    new_distances = np.expand_dims(new_distances, axis = 0)
    distances = distances[mask]
    distances = distances[:, mask]
    distances = np.append(distances, new_distances, axis = 0)
    new_distances = np.append(new_distances, np.array([[0]]), axis = 1)
    distances = np.append(distances, new_distances.transpose(), axis = 1)
    return nodes, distances

def join_last_nodes(nodes, distances,last_node):
    root = Node(last_node, nodes)
    i_node = nodes[0]
    j_node = nodes[1]
    m_node = nodes[2]

    d_root_i = (distances[0,1] + distances[0,2] - distances[1,2])/2
    d_root_j = (distances[0,1] + distances[1,2] - distances[0,2])/2
    d_root_m = (distances[0,2] + distances[1,2] - distances[0,1])/2

    i_node.distance = d_root_i
    j_node.distance = d_root_j
    m_node.distance = d_root_m
    return root

def neighbour_joining(distances, taxa):
    #Initialize nodes for every taxa
    nodes = np.array([Node(i) for i in range(len(taxa))])
    n_nodes = len(taxa)

    while len(nodes) > 3:
        n_taxa = distances.shape[0]
        r = distances.sum(axis = 0)/(n_taxa-2)
        n = np.subtract(distances, r)
        n = np.subtract(n.transpose(), r).transpose()
        np.fill_diagonal(n,np.inf)
        min_pair = (n.argmin()//distances.shape[0], n.argmin()%distances.shape[0])
        nodes, distances = join(nodes, min_pair[0], min_pair[1], distances, r, n_nodes)
        n_nodes += 1

    root = join_last_nodes(nodes, distances, n_nodes)
    return root

def write_newick_to_file(filename, tree):
    with open(filename, "w") as f:
        f.write(tree)


if __name__ == '__main__':
    taxa, distances = read_distance_matrix("example_slide4.phy")
    root = neighbour_joining(distances, taxa)
    newick_format = root.write_in_newick(taxa)
    write_newick_to_file("example_slide.nwk", newick_format)