from openpharmacophore import pharmacophoric_elements as phe

import numpy as np
from plip.structure.preparation import PDBComplex
import pyunitwizard as puw

def protein_ligand_interactions(file_name):
    
    mol_system = PDBComplex()
    mol_system.load_pdb(file_name)
    mol_system.analyze()
    
    # Dictionary with interactions for each small molecule in the protein
    all_interactions = mol_system.interaction_sets 

    # Filter ions and molecules with less than 10 atoms
    for ligand in mol_system.ligands:
        n_atoms = len(ligand.atomorder)
        if ligand.type != "SMALLMOLECULE" or n_atoms <= 10:
            ligand_id = ":".join([ligand.hetid, ligand.chain, str(ligand.position)])
            all_interactions.pop(ligand_id)
    
    return all_interactions

def sb_pharmacophore_points(interactions, radius):
    
    radius = puw.quantity(radius, "angstroms")

    # List with all interactions
    interactions = interactions.all_itypes
    # list of pharmacophoric_elements
    points = [] 
    for interaction in interactions:
        interaction_name = type(interaction).__name__
        
        if interaction_name == "pistack":
            ligand_center = np.array(interaction.ligandring.center)
            protein_center = np.array(interaction.proteinring.center)
            direction = protein_center - ligand_center
            aromatic = phe.AromaticRingSphereAndVector(
                center=puw.quantity(ligand_center, "angstroms"),
                radius=radius,
                direction=direction
            )
            points.append(aromatic)

        elif interaction_name == "hydroph_interaction":
            center = puw.quantity(interaction.ligatom.coords, "angstroms")
            hydrophobic = phe.HydrophobicSphere(
                center=center,
                radius=radius
            )
            points.append(hydrophobic)
        
        elif interaction_name == "saltbridge":
            if interaction.protispos:
                # The ligand has a negative charge
                center = puw.quantity(interaction.negative.center, "angstroms")
                charge_sphere = phe.NegativeChargeSphere(
                    center=center,
                    radius=radius
                ) 
            else:
                # The ligand has a positive charge
                center = puw.quantity(interaction.positive.center, "angstroms")
                charge_sphere = phe.PositiveChargeSphere(
                    center=center,
                    radius=radius
                ) 
            points.append(charge_sphere)
        
        elif interaction_name == "hbond":
            if interaction.protisdon:
                # The ligand has an acceptor atom
                ligand_acceptor_center = np.array(interaction.a.coords)
                protein_donor_center = np.array(interaction.d.coords)
                direction = ligand_acceptor_center - protein_donor_center 
                acceptor = phe.HBAcceptorSphereAndVector(
                    center=puw.quantity(ligand_acceptor_center, "angstroms"),
                    radius=radius,
                    direction=direction
                )
                points.append(acceptor)
            else:
                # The ligand has a donor atom
                ligand_donor_center = np.array(interaction.d.coords)
                protein_acceptor_center = np.array(interaction.a.coords)
                direction = protein_acceptor_center - ligand_donor_center
                donor = phe.HBDonorSphereAndVector(
                    center=puw.quantity(ligand_donor_center, "angstroms"),
                    radius=radius,
                    direction=direction
                )
                points.append(donor)
        else:
            # TODO: Incorporate other features such as halogenbonds, waterbridges and metal complexes
            continue

    # Remove duplicate points
    points_filtered = []
    for p in points:
        if p not in points_filtered:
            points_filtered.append(p)
    
    return points_filtered