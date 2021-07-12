import networkx as nx
import numpy as np
import numpy.linalg as LA
import math
from numpy import copy

def mean_bond_angles(mol_graph, verbose=False):
    
    """Compute mean bond angle for each central atom"""
    
    n_centroid_atoms = count_atoms_with_n_bonds(1, mol_graph) # Count number of central atoms
    centroid_mean_angle = np.zeros((n_centroid_atoms,2)) # Initialize array

    i = 0
    for atom, neighbors in mol_graph.adjacency(): #adjacency returns an iterator for the neighbors

        atom_names = [mol_graph.nodes[atom]['atom_name']]  # List with atom names
        centroid = np.array(mol_graph.nodes[atom]['coordinates'])
        n_neighbors = len(list(neighbors))

        if n_neighbors > 1:

            vertices = np.zeros((n_neighbors, 3))
            for j, nbr_atom in enumerate(list(neighbors)):
                vertices[j, :] = np.array(mol_graph.nodes[nbr_atom]['coordinates'])
                atom_names.append(mol_graph.nodes[nbr_atom]['atom_name'])

            angles = get_bond_angles(centroid, vertices, atom_names, verbose=verbose)
            centroid_mean_angle[i, 0] = atom
            centroid_mean_angle[i, 1] = np.mean(angles)

        else:
            continue

        i += 1
        
    return centroid_mean_angle

def get_bond_angles(centroid, vertices, atom_names=[], verbose=False):
    
    """Compute bond angle between a central atom and two adjacent atoms
        Bond angles are calculated for each combination of central atom and adjacent atoms
    """
    
    bond_vectors = vertices - centroid
    bond_lenght = LA.norm(bond_vectors, axis=-1)

    # print(bond_vectors, '\n')

    ii = bond_vectors.shape[0]
    
    m = nCr(ii, 2)
    angles = np.zeros((m,1))
    
    k = 0
    # Calculate angle between each different pair of vectors
    for i in range(ii - 1):
        for j in range(i + 1, ii):
            angle = np.degrees(np.arccos(np.dot(bond_vectors[i, :], bond_vectors[j, :]) / (bond_lenght[i] * bond_lenght[j])))
            angles[k] = angle
            k += 1
            if verbose and len(atom_names) > 0:
                print("Angle between atom {}, atom {} and atom {} is: {:.2f}".format(atom_names[i+1],
                                                             atom_names[0],
                                                             atom_names[j+1],              
                                                             angle))
    k += 1
    
    return angles

def count_atoms_with_n_bonds(n, mol_graph):
    """Count atoms with more than n bonds"""
    count = 0
    for atom, n_bonds in mol_graph.degree():
        if n_bonds > n:
            count += 1
        else:
            continue
    return count

def nCr(n,r):
    """Compute number of combinations"""
    f = math.factorial
    return f(n) // f(r) // f(n-r)