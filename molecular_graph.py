import molsysmt as msm
import networkx as nx
import numpy as np
import numpy.linalg as LA
from numpy import copy

def molecular_graph(molecular_system, selection='all', reset_index=True):
    """
       Transform a molecular system to a Networkx multigraph.
       Nodes are named by atom indices and they contain atom names as data.
       Edges represent bonds and they contain bond lenghts
    """
    Graph = nx.Graph()

    indices, types, names, bonded_atoms, coords = msm.get(molecular_system, selection=selection, target="atom", 
                                           atom_index=True, atom_type=True, name=True, inner_bonded_atoms=True, coordinates=True)
    
    cero_indx = True
    if not selection == 'all' or not selection == 'molecule_index==0':
        cero_indx = False #If selection isn't 'all' the indices array will not start at 0, causing an error in the distance function
    
    distances = get_bond_lenght(bonded_atoms, coords, indices, cero_indx=cero_indx)
    
    if reset_index:
        indx_dict = {value:index for index,value in enumerate(indices)} # dictionary with keys as the array value and values as the array index
        for k, v in indx_dict.items(): bonded_atoms[bonded_atoms==k] = v
        indices = np.arange(indices.shape[0]) # Indices now starts at cero
        
    
    for i in range(len(indices)):
        Graph.add_node(indices[i], atom_name=names[i], atom_type=types[i], coordinates=coords[0, i, :])
        Graph.add_edge(bonded_atoms[i,0], bonded_atoms[i,1], bond_lenght=distances[i][0])
    
    return Graph

def get_bond_lenght(bonded_atoms, coordinates, indices, cero_indx=False):
    """Computes the distance between bonded atoms"""
    distance = np.zeros(bonded_atoms.shape[0]) # initialize distance array
    coordinates = np.array(coordinates[0,:,:]) # convert to array to avoid error in numpy norm function
    
    bonded_atoms_copy = copy(bonded_atoms) # create a coppy of the array so it's not modified
    
    if not cero_indx:
       
        indx_dict = {value:index for index,value in enumerate(indices)} # dictionary with keys as the array value and values as the array index
        for k, v in indx_dict.items(): bonded_atoms_copy[bonded_atoms==k] = v
    
    coords_between_bonded_atoms = coordinates[bonded_atoms_copy] #matrix with pairs of coordiantes between bonded atoms
    bond_vectors = np.diff(coords_between_bonded_atoms, axis=1) #compute vectors for each pair of coordinates

    distance = LA.norm(bond_vectors, axis=-1)
    
    return distance