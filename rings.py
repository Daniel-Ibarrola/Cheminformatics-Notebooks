import networkx as nx
import collections

def find_rings(molecular_graph):
    """Find all rings in a molecular graph"""
    rings = []
    for node in list(molecular_graph.nodes()):
        ringset = get_ring(molecular_graph, node)
        if ringset and ringset not in rings:
            rings.append(ringset)
    
    return rings

def get_ring(molecular_graph, root):
    """Find a ringset in a molecular graph starting from a root node. The output depends on the root node"""
    ringset = []
    queue = collections.deque([root, root])
    path = [None] * (len(molecular_graph.nodes())) #Initialize paths to null
    path[root] = [root]

    while queue:

        # Dequeue a node from front of queue
        front_node = queue.popleft()
        source = queue.popleft()

        for neighbor in molecular_graph.adj[front_node]:

            if neighbor == source:
                continue

            if path[neighbor] is None:
                path[neighbor] = path[front_node] + [neighbor]
                path[neighbor].sort()
                queue.append(neighbor)
                queue.append(front_node)

            elif path[neighbor] is not None:
                # Compute intersection
                intersect = intersection(path[front_node], path[neighbor])
                if len(intersect) == 1: # If the intersection is a singleton there is a valid ring closure
                    ringset = Union(path[front_node], path[neighbor])
                    break

        if ringset:
            break
    
    return ringset
    
def intersection(lst1, lst2):
    """Compute the intersection of two lists"""
    return list(set(lst1) & set(lst2))

def Union(lst1, lst2):
    """Union of two lists"""
    final_list = list(set(lst1) | set(lst2))
    return final_list