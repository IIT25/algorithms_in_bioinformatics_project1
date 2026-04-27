import numpy as np

def write_newick_to_file(filename, tree):
    """
    Helper function to write the newick tree to a .nwk file
    """
    with open(filename, "w") as f:
        f.write(tree)

def read_distance_matrix(filename):
    """
    Helper function to read the distance matrix and the
    names of the taxa, recieves:
        - Filename: The name of the file where the distance matrix is 
            stored
    """
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


class Node:
    """
    Class to store the nodes and tree created by the neighbor
    joining algorithm
    """
    def __init__(self, name, children = None):
        """
        Constructor of the class, recieves:
            - name: name of the node, related to the taxon
            - children: List of children of the new node created
        
        In addition, the following attributes are initialized:
            - distance: Distance from the node to the parent
            - parent: Parent of the node
        Both attributes are initialized as None, since at first, the
        node doesn't have a parent.
        """
        self.name = name
        self.children = children if children is not None else []
        self.distance = None
        self.parent = None
    
    def write_in_newick(self, taxa):
        """
        Function to write the tree in newick format
        Does a recursive traverse of the tree to write the tree
        """
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


def join(nodes,i,j,distances,r,last_node):
    """
    One of the main functions of the neighbor joining algorithm
    Performs one step of the algorithm joining the two nodes selected
    by the algorithm and updating the nodes list and distance matrix.

    Recieves:
        - nodes: List of current nodes 
        - i: index of the first node to be joined
        - j: index of the second node to be joined
        - distances: Current distance matrix
        - r: r vector
        - last_node: name of the last created node, to continue with the sequence
    """
    mask = np.ones(len(nodes), dtype=bool)
    mask[[i, j]] = False

    new_node_name = last_node
    i_node = nodes[i]
    j_node = nodes[j]
    for idx in sorted([i, j], reverse=True):
        nodes.pop(idx)

    d_ij = distances[i,j]
    distance_i = (d_ij+r[i]-r[j])/2
    distance_j = d_ij - distance_i
    i_node.distance = distance_i
    j_node.distance = distance_j
    
    nodes.append(Node(new_node_name, [i_node, j_node]))
    i_node.parent = nodes[-1]
    j_node.parent = nodes[-1]

    new_distances = (distances[i] + distances[j] - d_ij) / 2
    new_distances = new_distances[mask]
    new_distances = np.expand_dims(new_distances, axis = 0)

    new_size = len(distances) - 1
    new_matrix = np.zeros((new_size, new_size))
    new_matrix[:-1, :-1] = distances[mask][:, mask]
    new_matrix[-1, :-1] = new_distances
    new_matrix[:-1, -1] = new_distances

    return nodes, new_matrix

def join_last_nodes(nodes, distances,last_node):
    """
    Last function in the neighbor joining procedure
    Joins the last 3 nodes of the tree, assigning the right distances
    and creating the root of the tree.
    Recieves:
        - nodes: 
    """
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
    """
    Neighbor joining procedure, iteratively runs the algorithm
    selecting with the r vector and the Q (also called n) matrix
    the nodes to be joined until there are only 3 remaining.

    Recieves:
        - distances: Initial distance matrix
        - taxa: List with the name of the taxa
    """
    
    nodes = [Node(i) for i in range(len(taxa))]
    n_nodes = len(taxa)

    while len(nodes) > 3:
        n_taxa = distances.shape[0]
        r = distances.sum(axis = 0)/(n_taxa-2)
        n = distances - r[:, None] - r[None, :]
        np.fill_diagonal(n,np.inf)
        idx = n.argmin()
        min_pair = (idx // n.shape[0], idx % n.shape[0])
        nodes, distances = join(nodes, min_pair[0], min_pair[1], distances, r, n_nodes)
        n_nodes += 1

    root = join_last_nodes(nodes, distances, n_nodes)
    return root

if __name__ == '__main__':
    taxa, distances = read_distance_matrix("example_slide4.phy")
    root = neighbour_joining(distances, taxa)
    newick_format = root.write_in_newick(taxa)
    write_newick_to_file("example_slide.nwk", newick_format)