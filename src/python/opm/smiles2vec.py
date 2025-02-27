import datamol as dm

from rdkit import DataStructs
from rdkit import Chem

import numpy as np

def is_valid_config(config:dict)->bool:
    try:
        import rdkit
    except ImportError:
        raise ImportError("rdkit is not installed")

    return True

def get_embedding(smiles:str,config:dict)->list:

    mol = dm.to_mol(smiles)
    fp = dm.to_fp(mol,fp_type='maccs')
    print(fp.shape)
    return str(fp.tolist())

if __name__ == "__main__":
    smiles_a = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"
    smiles_b = "CC(=O)O[C@H]1C=C[C@H]2[C@H]3CC4=C5[C@]2([C@H]1OC5=C(C=C4)OC(=O)C)CCN3C"

    emb_a = get_embedding('O', {})
    emb_b = get_embedding(smiles_b, {})
    print(emb_a)
    # # 
    # factory = Gobbi_Pharm2D.factory
    # fp1 = Generate.Gen2DFingerprint(mol_a, factory)
    # fp2 = Generate.Gen2DFingerprint(mol_b, factory)
    # tmp = np.zeros(fp1.GetNumBits(), dtype=int)
    # on_bits = np.array(fp1.GetOnBits())
    # tmp[on_bits] = 1
    # print(tmp)

    # simi = DataStructs.TanimotoSimilarity(fp1, fp2)
    # print(simi)
    