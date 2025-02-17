from io import BytesIO
import codecs
from rdkit import Chem
from rdkit.Chem import Draw

from msg import SmilesRequest, SmilesResponse, CreateImageRequest, CreateImage2SdfRequest, CreateImage2SmilesRequest

from grongier.pex import BusinessOperation

import requests

import iris

class IUPACOperation(BusinessOperation):
    """
    IUPACOperation is a chemical operation that returns the IUPAC name of the molecule.
    """
    def process(self, request:SmilesRequest) -> SmilesResponse:
        """
        Processes the molecule and returns the IUPAC name of the molecule.

        :param mol: The molecule to process.
        :return: The IUPAC name of the molecule.
        """
        rsp = SmilesResponse()
        mol = Chem.MolFromSmiles(request.smiles)
        rsp.properties = {}
        rsp.properties['iupac_name'] = self._calculate_iupac_name(mol)
        return rsp

    def _calculate_iupac_name(self,mol):
        # This code uses the CACTUS web service to convert a SMILES string
        # to an IUPAC name. It returns the name as a string.
        CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
        rep = "iupac_name"
        smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
        url = CACTUS.format(smiles, rep)
        # try catch block to handle timeout
        try:
            response = requests.get(url,timeout=10)
        except requests.exceptions.Timeout:
            self.log_warning("Timeout occurred")
            return None
        rsp = None
        if (response.status_code == 200):
            rsp = response.text
        else:
            self.log_warning("Error occurred with status code: " + str(response.status_code))
        return rsp

class GenerateImageOperation(BusinessOperation):
    """
    GenerateImageOperation is a chemical operation that generates an image of the molecule.
    """

    def process_2_sdf(self, request:CreateImage2SdfRequest):
        """
        Processes the molecule and returns an image of the molecule.

        :param mol: The molecule to process.
        :return: An image of the molecule.
        """
        image = None
        if request.filename_a is not None and request.filename_b is not None:
            image = self._draw_molecule_2_sdf_diff(request.filename_a,request.filename_b)
        else:
            return None

        resp = iris.cls('Opm.ImageDisplay')._New()
       
        # Converting the image into a binary format and then writing it into the 
        # BinaryImage field of the response.
        output = BytesIO()
        image.save(output, format="png")
        binary = output.getvalue()
        base64 = codecs.encode(binary, 'base64').decode()
        buffer = 3600
        chunks = [binary[i:i+buffer] for i in range(0, len(binary), buffer)]
        for chunk in chunks:
            resp.BinaryImage.Write(chunk)
        base64_chunks = [base64[i:i+buffer] for i in range(0, len(base64), buffer)]
        for base64_chunk in base64_chunks:
            resp.Base64Image.Write(base64_chunk)

        return resp
    
    def process_2_smiles(self, request:CreateImage2SmilesRequest):
        """
        Processes the molecule and returns an image of the molecule.

        :param mol: The molecule to process.
        :return: An image of the molecule.
        """
        image = None
        if request.smiles_a is not None and request.smiles_b is not None:
            image = self._draw_molecule_2_smiles_diff(request.smiles_a,request.smiles_b)
        else:
            return None

        resp = iris.cls('Opm.ImageDisplay')._New()
       
        # Converting the image into a binary format and then writing it into the 
        # BinaryImage field of the response.
        output = BytesIO()
        image.save(output, format="png")
        binary = output.getvalue()
        base64 = codecs.encode(binary, 'base64').decode()
        buffer = 3600
        chunks = [binary[i:i+buffer] for i in range(0, len(binary), buffer)]
        for chunk in chunks:
            resp.BinaryImage.Write(chunk)
        base64_chunks = [base64[i:i+buffer] for i in range(0, len(base64), buffer)]
        for base64_chunk in base64_chunks:
            resp.Base64Image.Write(base64_chunk)

        return resp


    def process(self, request:CreateImageRequest):
        """
        Processes the molecule and returns an image of the molecule.

        :param mol: The molecule to process.
        :return: An image of the molecule.
        """
        image = None
        if request.filename is not None and request.smiles is not None:
            image = self._draw_molecule_diff(request.smiles,request.filename)
        elif request.smiles is not None:
            image = self._draw_molecule_smiles(request.smiles)
        elif request.filename is not None:
            image = self._draw_molecule_sdf(request.filename)
        else:
            return None

        resp = iris.cls('Opm.ImageDisplay')._New()
       
        # Converting the image into a binary format and then writing it into the 
        # BinaryImage field of the response.
        output = BytesIO()
        image.save(output, format="png")
        binary = output.getvalue()
        base64 = codecs.encode(binary, 'base64').decode()
        buffer = 3600
        chunks = [binary[i:i+buffer] for i in range(0, len(binary), buffer)]
        for chunk in chunks:
            resp.BinaryImage.Write(chunk)
        base64_chunks = [base64[i:i+buffer] for i in range(0, len(base64), buffer)]
        for base64_chunk in base64_chunks:
            resp.Base64Image.Write(base64_chunk)

        return resp
    
    def _draw_molecule_2_smiles_diff(self,smiles_a,smiles_b):
        """
        _draw_molecule_2_smiles_diff takes two smiles strings and returns the image of the difference between the two molecules in bytes.
        smiles_a: smiles string of the first molecule to be drawn
        smiles_b: smiles string of the second molecule to be drawn
        size: size of the image to be drawn
        """
        mol_a = Chem.MolFromSmiles(smiles_a)
        mol_b = Chem.MolFromSmiles(smiles_b)

        return Draw.MolsToGridImage([mol_a,mol_b],
                                    legends=['Smiles A','Smiles B'],
                                    highlightAtomLists=[mol_a.GetSubstructMatch(mol_b),mol_b.GetSubstructMatch(mol_a)],
                                    molsPerRow=1)
    
    def _draw_molecule_2_sdf_diff(self,filename_a,filename_b):
        """
        _draw_molecule_2_sdf_diff takes two sdf files and returns the image of the difference between the two molecules in bytes.
        filename_a: sdf file of the first molecule to be drawn
        filename_b: sdf file of the second molecule to be drawn
        size: size of the image to be drawn
        """
        mols_a = Chem.SDMolSupplier(filename_a)
        mols_b = Chem.SDMolSupplier(filename_b)
        mol_a = mols_a[0]
        mol_b = mols_b[0]
        # Highlight the differences between the two molecules
        mol_a = Chem.MolFromSmiles(Chem.MolToSmiles(mol_a,isomericSmiles=True))
        mol_b = Chem.MolFromSmiles(Chem.MolToSmiles(mol_b,isomericSmiles=True))
        return Draw.MolsToGridImage([mol_a,mol_b],molsPerRow=1,legends=['SDF A','SDF B'],highlightAtomLists=[mol_a.GetSubstructMatch(mol_b),mol_b.GetSubstructMatch(mol_a)])

    def _draw_molecule_diff(self,smiles,filename):
        """
        _draw_molecule_diff takes a smiles string and sdf file and returns the image of the difference between the two molecules in bytes.
        smiles: smiles string of the molecule to be drawn
        filename: sdf file of the molecule to be drawn
        size: size of the image to be drawn
        """
        mol_smiles = Chem.MolFromSmiles(smiles)
        mols = Chem.SDMolSupplier(filename)
        mol_sdf = mols[0]
        # Highlight the differences between the two molecules
        mol_sdf = Chem.MolFromSmiles(Chem.MolToSmiles(mol_sdf,isomericSmiles=True))
        mol_smiles = Chem.MolFromSmiles(Chem.MolToSmiles(mol_smiles,isomericSmiles=True))
        return Draw.MolsToGridImage([mol_smiles,mol_sdf],molsPerRow=1,legends=['Smiles','SDF'],highlightAtomLists=[mol_smiles.GetSubstructMatch(mol_sdf),mol_sdf.GetSubstructMatch(mol_smiles)])

    def _draw_molecule_smiles(self,smiles):
        """
        _draw_molecule_smiles takes a smiles string and returns the image of the molecule in bytes.
        smiles: smiles string of the molecule to be drawn
        size: size of the image to be drawn
        """
        mol = Chem.MolFromSmiles(smiles)
        return Draw.MolToImage(mol)

    def _draw_molecule_sdf(self,filename):
        """
        _draw_molecule_sdf takes a sdf file and returns the image of the molecule in bytes.
        filename: sdf file of the molecule to be drawn
        size: size of the image to be drawn
        """
        mols = Chem.SDMolSupplier(filename)
        mol = mols[0]
        return Draw.MolToImage(mol)


if __name__ == '__main__':
    bo = GenerateImageOperation()
    msg = CreateImageRequest(smiles='CC(C)Cc1ccc(C(C)C(=O)O)cc1',filename='test.png')
    bo.process(msg)