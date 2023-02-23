from rdkit import Chem

# read SD file
mols = Chem.SDMolSupplier("misc/data/ibuprofen-3D-structure-CT1078642946.sdf")

# get the smiles
sm = Chem.MolToSmiles(mols[0])

