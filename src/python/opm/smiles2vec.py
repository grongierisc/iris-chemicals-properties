import datamol as dm

def is_valid_config(config:dict)->bool:
    try:
        import rdkit
    except ImportError:
        raise ImportError("rdkit is not installed")

    return True

def get_embedding(smiles:str,config:dict)->list:

    mol = dm.to_mol(smiles)
    fp = dm.to_fp(mol,fp_type='maccs')

    return str(fp.tolist())