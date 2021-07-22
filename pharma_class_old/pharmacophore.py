from dbscan import get_feature_clusters
from visualize import view_pharmacophore
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
from centroid import feature_centroid
import numpy as np
import os

class Pharmacophore():

    def __init__(self, ligands):
        self.features = {}
        self.num_features = len(self.features)
        self.ligands = ligands # List of rdkit molecules
    
    def common_pharmacophore(self, method, feat_list=None):
        """
            Computes the common pharmacophore.
            Returns a dictionary which keys are features and
            items are the coordinates of the features.
        """
        method = method.lower()

        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        
        if not feat_list:
            feat_list = ['Acceptor', 'Aromatic', 'Donor', 'Hydrophobe', 'PosIonizable', 'NegIonizable']
        
        feat_coords = {}
            
        for feature in feat_list:
            feat_coords[feature] = []
            
            for ligand in self.ligands:
            
                feats = factory.GetFeaturesForMol(ligand, includeOnly=feature)
                n_conformers = ligand.GetNumConformers()
                if n_conformers == 0: 
                    n_conformers == 1
                    
                for f in feats:

                    atom_idxs = f.GetAtomIds()
                    
                    for conformer in range(n_conformers):
                    
                        if len(atom_idxs) > 1: # Get the centroid of that feature
                            coords = feature_centroid(ligand, atom_idxs, conformer)
                        else:
                            position = ligand.GetConformer(conformer).GetAtomPosition(atom_idxs[0])
                            coords = np.zeros((3,))
                            coords[0] = position.x
                            coords[1] = position.y
                            coords[2] = position.z

                        feat_coords[feature].append(coords.tolist())
                
            feat_coords[feature] = np.array(feat_coords[feature])
            

            if method == "dbscan":
                self.features = get_feature_clusters(feat_coords, eps=2, min_samples=0.75*len(self.ligands))
            else:
                raise NotImplementedError
    
    def visualize(self, show_ligands=True):
       view = view_pharmacophore(self.features, self.ligands, show_ligands=show_ligands)
       return view