import numpy as np
import os
from sklearn.cluster import DBSCAN
from rdkit import RDConfig, Chem
from rdkit.Chem import AllChem, ChemicalFeatures, Draw
from rdkit.Chem import rdMolAlign

def get_feature_clusters(feat_coords, eps, min_samples):
    
    clusters = {}
    for feat, coords in feat_coords.items():
        db_scan = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
        
        labels = db_scan.labels_
        core_samples_mask = np.zeros_like(labels, dtype=bool)
        core_samples_mask[db_scan.core_sample_indices_] = True
        
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = list(labels).count(-1)

        # print(f'Estimated number of clusters: {n_clusters} for {feat} feature' )
        # print(f'Estimated number of noise points: {n_noise}\n')
        
        centroids = []
        unique_labels = set(labels)
        for k in unique_labels:
            if k == -1:
                continue
            class_member_mask = (labels == k)
            cluster = feat_coords[feat][class_member_mask & core_samples_mask]
            cluster_centroid = cluster.mean(axis=0)
            centroids.append(cluster_centroid.tolist())
        
        clusters[feat] = centroids
    
    return clusters