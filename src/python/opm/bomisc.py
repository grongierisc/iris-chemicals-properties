from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw

from msg import SmilesRequest, SmilesResponse, CreateImageRequest

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
            return None
        rsp = None
        if (response.status_code == 200):
            rsp = response.text
        return rsp

class GenerateImageOperation(BusinessOperation):
    """
    GenerateImageOperation is a chemical operation that generates an image of the molecule.
    """
    def process(self, request:CreateImageRequest):
        """
        Processes the molecule and returns an image of the molecule.

        :param mol: The molecule to process.
        :return: An image of the molecule.
        """
        mol = Chem.MolFromSmiles(request.smiles)
        resp = iris.cls('Opm.ImageDisplay')._New()

        image= self._draw_molecule(mol)
        # Converting the image into a binary format and then writing it into the 
        # BinaryImage field of the response.
        output = BytesIO()
        image.save(request.filename, format="png")
        image.save(output, format="png")
        binary = output.getvalue()
        buffer = 3600
        chunks = [binary[i:i+buffer] for i in range(0, len(binary), buffer)]
        for chunk in chunks:
            resp.BinaryImage.Write(chunk)

        return resp

    def _draw_molecule(self,mol):
        """
        _draw_molecule takes a molecule and returns the image of the molecule in bytes.
        mol: molecule to be drawn
        size: size of the image to be drawn
        """
        return Draw.MolToImage(mol)


if __name__ == '__main__':
    bo = GenerateImageOperation()
    msg = CreateImageRequest(smiles='CC(C)Cc1ccc(C(C)C(=O)O)cc1',filename='test.png')
    bo.process(msg)