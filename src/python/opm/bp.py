from dataclasses import dataclass
import requests

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors


from grongier.pex import BusinessProcess,Message

@dataclass
class SmilesRequest(Message):
    smiles:str = None

@dataclass
class SmilesResponse(Message):
    image:bytes = None
    properties:dict = None

class ChemBase:
    """
    This class is used to 
        get the properties of a molecule
        get the image of a molecule
    """

    def get_properties(self, mol:Chem.Mol) -> dict:
        properties = {}
        # Calculate SMILES
        smiles = Chem.MolToSmiles(mol)

        # Calculate IUPAC name
        properties["IUPAC name"] = self.smiles_to_iupac(smiles)

        # Calculate formula
        properties["Formula"] = Descriptors.rdMolDescriptors.CalcMolFormula(mol)

        # Calculate molecular weight
        properties["MW"] = Descriptors.MolWt(mol)

        # Calculate SMILES
        properties["SMILES"] = smiles

        # Calculate ClogP the calculated octanol/water partition coefficient
        properties["ClogP"] = Descriptors.MolLogP(mol)

        # Calculate ClogD the calculated distribution coefficient
        # To discuss with the team
        properties["ClogD"] = Descriptors.MolLogP(mol, True)

        # Calculate TPSA
        properties["TPSA"] = Descriptors.TPSA(mol)

        # Calculate pKa Basic 1st strongest
        # To discuss with the team
        # properties["pKa Basic 1st strongest"] = self.get_strongest_basic_pka(mol)

        # Calculate H donor
        properties["H donor"] = Descriptors.NumHDonors(mol)

        # Calculate H acceptor
        properties["H acceptor"] = Descriptors.NumHAcceptors(mol)

        # Calculate heavy atom count
        properties["Heavy atom count"] = Descriptors.HeavyAtomCount(mol)

        # Calculate rotatable bonds
        properties["Rotatable bonds"] = Descriptors.NumRotatableBonds(mol)

        # Calculate MPO (Molecular Property Optimizer)
        mpo = self.get_mpo(mol)
        properties["mpo"] = mpo

        # Calculate MPO_ClogP
        properties["mpo_clogP"] = mpo.GetDescriptors(['MolLogP'])[0]

        # Calculate MPO_ClogD
        properties["mpo_clogD"] = Descriptors.MolProperty(mol, Descriptors.MolLogD)

        # Calculate MPO_TPSA
        properties["mpo_TPSA"] = Descriptors.MolProperty(mol, Descriptors.TPSA)

        # Calculate MPO_pKa
        properties["mpo_pKa"] = Descriptors.MolProperty(mol, Descriptors.pKa)

        # Calculate MPO_HBD
        properties["mpo_HBD"] = Descriptors.MolProperty(mol, Descriptors.NumHDonors)

        # Calculate MPO_MW
        properties["mpo_MW"] = Descriptors.MolProperty(mol, Descriptors.MolWt)

        return properties

    def get_mpo(self, mol:Chem.Mol) -> dict:
        """
        This function is used to get the MPO of a molecule
        """
        # Calculate the MPO
        

        return Chem.rdMolDescriptors.MolPropertyCalculator(mol)


    def get_pka(self, mol:Chem.Mol) -> dict:
        """
        This function is used to get the pKa of a molecule
        """
        # Calculate the number of acidic and basic hydrogens
        num_acidic_h = Chem.rdMolDescriptors.CalcNumAcidicH(mol)
        num_basic_h = Chem.rdMolDescriptors.CalcNumBasicH(mol)

        # Calculate the pKa using the Henderson-Hasselbalch equation
        pKa = 7.0 + (num_basic_h - num_acidic_h)
        return pKa

    def get_strongest_basic_pka(self,mol):
        pKa_values = rdMolStandardize.Normalizer().normalize(mol)

        # Create a list of ionizable groups and their predicted pKa values
        ionizable_groups = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
                ionizable_groups.append(('N', atom.GetIdx(), pKa_values['pKa'][atom.GetIdx()]))
            elif atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
                ionizable_groups.append(('O', atom.GetIdx(), pKa_values['pKa'][atom.GetIdx()]))

        # Sort the ionizable groups by predicted pKa value in descending order
        ionizable_groups = sorted(ionizable_groups, key=lambda x: -x[2])

        # Print the first (most basic) ionizable group
        if len(ionizable_groups) > 0:
            print(f'The first (most basic) ionizable group is a {ionizable_groups[0][0]} atom at index {ionizable_groups[0][1]} with a predicted pKa of {ionizable_groups[0][2]:.2f}')
        else:
            print('No ionizable groups found in the molecule')


    def smiles_to_iupac(self,smiles):
        # This code uses the CACTUS web service to convert a SMILES string
        # to an IUPAC name. It returns the name as a string.
        CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
        rep = "iupac_name"
        url = CACTUS.format(smiles, rep)
        response = requests.get(url)
        response.raise_for_status()
        return response.text


    def get_image(self, mol:Chem.Mol) -> bytes:
        """
        This function is used to get the image of a molecule
        """
        img = Draw.MolToImage(mol)
        # save the image on disk
        img.save("mol.png")
        # return the image as a byte array
        return img.tobytes()

class SmilesProcess(BusinessProcess, ChemBase):
    """
    This class is used to process Smiles data
    It generate the image and the properties of the molecule
    """
    def on_message(self, request:SmilesRequest) -> SmilesResponse:
        """
        Main function of this process
        Make create the image as a byte array
        Make the dictionary of properties
        """
        # sm data
        sm = request.smiles
        # sm to mol
        mol = Chem.MolFromSmiles(sm)
        # mol to image
        img = self.get_image(mol)
        # mol to properties
        prop = self.get_properties(mol)
        # return the response
        return SmilesResponse(image=img, properties=prop)

class SdFileProcess(BusinessProcess, ChemBase):
    """
    This class is used to process SD file
    It generate the image and the properties of the molecule
    """
    def on_message(self, request:SmilesRequest) -> SmilesResponse:
        """
        Main function of this process
        Make create the image as a byte array
        Make the dictionary of properties
        """
        # read SD file
        mols = Chem.SDMolSupplier(request.smiles)
        # mol to image
        img = self.get_image(mols[0])
        # mol to properties
        prop = self.get_properties(mols[0])
        # return the response
        return SmilesResponse(image=img, properties=prop)

if __name__ == "__main__":
    # test SmilesProcess
    sp = SmilesProcess()
    rsp = sp.on_message(SmilesRequest(smiles="CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"))
    print(rsp.properties)
    # test SdFileProcess
    # sp = SdFileProcess()
    # sp.on_message(SmilesRequest(smiles="data/1.sdf"))