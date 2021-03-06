from rdkit import RDConfig, Chem
from rdkit.Chem import ChemicalFeatures
from openpharmacophore.errors import NoConformersError, PointTypeError
from openpharmacophore.utils.centroid import feature_centroid
from openpharmacophore.utils.rdkit_to_point import rdkit_to_point
import numpy as np
import os

def conformer_pharmacophoric_points(ligand, conformer_idx=0, feat_list=None, feat_def="rdkit",
                                    point_type="spheres", radius=None):
    """
        Get pharmacophoric points for a single conformer. 
    """

    point_type_list = ["spheres", "spheres_vectors", "gaussian", "shapelet"]
    if point_type not in point_type_list:
        raise PointTypeError(point_type)

    if point_type == "spheres" or point_type == "spheres_vectors" and radius is None:
        raise ValueError("Radius cannot be null if point type is spheres or spheres and vectors") 

    if feat_def == 'rdkit':
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    else:
        raise NotImplementedError
    
    if not feat_list:
        feat_list = ['Acceptor', 'Aromatic', 'Donor', 'Hydrophobe', 'PosIonizable', 'NegIonizable']
    
    n_conformers = ligand.GetNumConformers()
    if n_conformers == 0:
        raise NoConformersError(n_conformers)
    
    feats = factory.GetFeaturesForMol(ligand)
    points = []
        
    for f in feats:
        
        feat_name = f.GetFamily()
        if feat_name not in feat_list:
            continue
        atom_idxs = f.GetAtomIds()
                 
        if len(atom_idxs) > 1: # Get the centroid of that feature
            coords = feature_centroid(ligand, atom_idxs, conformer_idx)
        else:
            position = ligand.GetConformer(conformer_idx).GetAtomPosition(atom_idxs[0])
            coords = np.zeros((3,))
            coords[0] = position.x
            coords[1] = position.y
            coords[2] = position.z

        point = rdkit_to_point(feat_name, coords, radius=radius, point_type=point_type)
        points.append(point)
    
    return points