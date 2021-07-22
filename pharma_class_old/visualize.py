import nglview as nv
from rdkit import Chem

def view_ligands(molecules):
    """
    Generate a view of the ligand molecules.

    Parameters
    -----------
    molecules: an rdkit.Chem.rdchem.Mol or a list of rdkit.Chem.rdchem.Mol

    Returns
    ----------
    nglview.widget.NGLWidget
    """
    if not isinstance(molecules, list):
        molecules = [molecules]
    
    view = nv.NGLWidget()
    for molecule in molecules:
        component = view.add_component(molecule)
        component.clear()
        component.add_ball_and_stick(multipleBond=True)
    return view

def view_pharmacophore(pharmacophore, molecules, show_ligands=True):
    """
    Generate a view of a pharmacophore.

    Parameters
    -----------
    pharmacophore: a dictionary of features with coordinates
    molecules: a list of rdkit.Chem.rdchem.Mol

    Returns
    ----------
    nglview.widget.NGLWidget
    """
    
    feature_colors = {
    "Acceptor": (0.90, 0.30, 0.24),  # Red
    "Aromatic": (1, 0.9, 0),  # Yellow
    "Donor": (0.13, 0.56, 0.30), # Green
    "Hydrophobic": (1, 0.9, 0),  # Yellow,
    "PosIonizable": (0.12, 0.36, 0.52), # Blue
    "NegIonizable": (0.90, 0.30, 0.24),  # Red
    }

    if show_ligands:
        view = view_ligands(molecules)
    else:
        view = nv.NGLWidget()
    
    sphere_radius = 1
    
    for featname, points in pharmacophore.items():
        for i, coords in enumerate(points):
            feature_color = feature_colors[featname]
            label = f"{featname}_{i}"
            view.shape.add_sphere(coords, feature_color, sphere_radius, label)
    
    return view


def view_conformers(molecule):
    """
    Generate a view of the conformers of a molecule.

    Parameters
    -----------
    molecule: an rdkit.Chem.rdchem.Mol

    Returns
    ----------
    nglview.widget.NGLWidget
    """
    view = nv.NGLWidget()
    for conformer in range(molecule.GetNumConformers()):
        mol_string = Chem.MolToMolBlock(molecule, confId=conformer)
        temp_mol = Chem.MolFromMolBlock(mol_string, removeHs=False)
        component = view.add_component(temp_mol)
        component.clear()
        component.add_ball_and_stick(multipleBond=True)
    return view