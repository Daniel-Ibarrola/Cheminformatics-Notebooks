from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

def confgen(input, output, prunermsthresh, numconf, add_ref):
    mol = Chem.AddHs(Chem.MolFromMolFile(input), addCoords=True)
    refmol = Chem.AddHs(Chem.Mol(mol))
    param = rdDistGeom.ETKDGv2()
    param.pruneRmsThresh = prunermsthresh
    cids = rdDistGeom.EmbedMultipleConfs(mol, numconf, param)
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant='MMFF94s')
    w = Chem.SDWriter(output)
    if add_ref:
        refmol.SetProp('CID', '-1')
        refmol.SetProp('Energy', '')
        w.write(refmol)
    res = []

    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))
    sorted_res = sorted(res, key=lambda x:x[1])
    rdMolAlign.AlignMolConformers(mol)
    for cid, e in sorted_res:
        mol.SetProp('CID', str(cid))
        mol.SetProp('Energy', str(e))
        w.write(mol, confId=cid)
    w.close()

if __name__=='__main__':
    confgen()