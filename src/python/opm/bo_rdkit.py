from grongier.pex import BusinessOperation

import requests

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
#from .pka.predictor import PkaPredictor

class RDKitOperation(BusinessOperation):
    """
    This class is used to 
        get the properties of a molecule
        get the image of a molecule
    """
    def _get_properties(self, mol:Chem.Mol) -> dict:
        """
        Returns a dictionary of properties for the molecule.

        :param mol: The molecule to calculate properties for.
        :return: A dictionary of properties for the molecule.
        """

        properties = {}

        properties["IUPAC name"] = self._calculate_iupac_name(mol)
        properties["Formula"] = self._calculate_mol_formula(mol)
        properties["MW"] = self._calculate_mol_weight(mol)
        properties["SMILES"] = self._calculate_smiles(mol)
        properties["ClogP"] = self._calculate_clogp(mol)
        properties["ClogD"] = self._calculate_clogd(mol)
        properties["TPSA"] = self._calculate_tpsa(mol)
        properties["pKa"] = self._calculate_pka(mol)
        properties["H donor"] = self._calculate_num_h_donor(mol)
        properties["H acceptor"] = self._calculate_num_h_acceptor(mol)
        properties["Heavy atom count"] = self._calculate_heavy_atom_count(mol)
        properties["Rotatable bonds"] = self._calculate_num_rotatable_bonds(mol)

        return properties

    def _calculate_mol_formula(self, mol:Chem.Mol) -> str:
        """
        Returns the molecular formula of the molecule.

        :param mol: The molecule to calculate the molecular formula for.
        :return: The molecular formula of the molecule.
        """
        return Descriptors.rdMolDescriptors.CalcMolFormula(mol)

    def _calculate_mol_weight(self, mol:Chem.Mol) -> float:
        """
        Returns the molecular weight of the molecule.

        :param mol: The molecule to calculate the molecular weight for.
        :return: The molecular weight of the molecule.
        """
        return Descriptors.MolWt(mol)

    def _calculate_smiles(self, mol:Chem.Mol) -> str:
        """
        Returns the SMILES string of the molecule.

        :param mol: The molecule to calculate the SMILES string for.
        :return: The SMILES string of the molecule.
        """
        return Chem.MolToSmiles(mol,isomericSmiles=True)

    def _calculate_clogp(self, mol:Chem.Mol) -> float:
        """
        Returns the ClogP value of the molecule.

        :param mol: The molecule to calculate the ClogP value for.
        :return: The ClogP value of the molecule.
        """
        return Descriptors.MolLogP(mol)

    def _calculate_clogd(self, mol:Chem.Mol) -> float:
        """
        Returns the ClogD value of the molecule.

        :param mol: The molecule to calculate the ClogD value for.
        :return: The ClogD value of the molecule.
        """
        return Descriptors.MolLogP(mol,True)

    def _calculate_tpsa(self, mol:Chem.Mol) -> float:
        """
        Returns the TPSA value of the molecule.

        :param mol: The molecule to calculate the TPSA value for.
        :return: The TPSA value of the molecule.
        """
        return Descriptors.TPSA(mol)

    def _calculate_pka(self, mol:Chem.Mol) -> list:
        """
        Returns a list of pKa values for the molecule.

        :param mol: The molecule to calculate pKa values for.
        :return: A list of pKa values for the molecule.
        """
        # pka_predictor = PkaPredictor()
        # pka = pka_predictor.predict(mol)
        return None

    def _calculate_num_h_donor(self, mol:Chem.Mol) -> int:
        """
        Returns the number of hydrogen donors in the molecule.

        :param mol: The molecule to calculate the number of hydrogen donors for.
        :return: The number of hydrogen donors in the molecule.
        """
        return Descriptors.NumHDonors(mol)

    def _calculate_num_h_acceptor(self, mol:Chem.Mol) -> int:
        """
        Returns the number of hydrogen acceptors in the molecule.

        :param mol: The molecule to calculate the number of hydrogen acceptors for.
        :return: The number of hydrogen acceptors in the molecule.
        """
        return Descriptors.NumHAcceptors(mol)

    def _calculate_heavy_atom_count(self, mol:Chem.Mol) -> int:
        """
        Returns the number of heavy atoms in the molecule.

        :param mol: The molecule to calculate the number of heavy atoms for.
        :return: The number of heavy atoms in the molecule.
        """
        return Descriptors.HeavyAtomCount(mol)

    def _calculate_num_rotatable_bonds(self, mol:Chem.Mol) -> int:
        """
        Returns the number of rotatable bonds in the molecule.

        :param mol: The molecule to calculate the number of rotatable bonds for.
        :return: The number of rotatable bonds in the molecule.
        """
        return Descriptors.NumRotatableBonds(mol)

    def _draw_molecule(mol, size=(300,300)):
        """
        _draw_molecule takes a molecule and returns the image of the molecule in bytes.
        mol: molecule to be drawn
        size: size of the image to be drawn
        """
        if mol is None:
            return None
        try:
            img = Draw.MolToImage(mol, size=size)
            return img.tobytes()
        except:
            return None

    def _calculate_iupac_name(self,smiles):
        # This code uses the CACTUS web service to convert a SMILES string
        # to an IUPAC name. It returns the name as a string.
        CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
        rep = "iupac_name"
        url = CACTUS.format(smiles, rep)
        response = requests.get(url,timeout=10)
        rsp = None
        if (response.status_code == 200):
            rsp = response.text
        return rsp

if __name__ == "__main__":
    rd = RDKitOperation()
    mols = Chem.SDMolSupplier("/Users/grongier/git/iris-chemicals-properties/misc/data/ibuprofen-3D-structure-CT1078642946.sdf")
    mol = mols[0]
    pp = mol.GetPropsAsDict()
    print(pp)

    properties = rd._get_properties(mol)
    print(properties)

    smi = 'CC(C)Cc1ccc(cc1)C(C)C(O)=O'
    mol = Chem.MolFromSmiles(smi)
    properties = rd._get_properties(mol)
    print(properties)