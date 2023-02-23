from grongier.pex import BusinessOperation

from .pka.predictor import PkaPredictor

from rdkit import Chem

from msg import PkaRequest, PkaResponse

class PkaPredictorOperation(BusinessOperation):

    def get_pka_from_smiles(self, request: PkaRequest) -> dict:
        """
        Returns a list of pKa values for the molecule.

        :param smiles: The smiles of the molecule to calculate pKa values for.
        :return: A list of pKa values for the molecule.
        """
        mol = Chem.MolFromSmiles(request.smiles)
        pka_dict = self._calculate_pka(mol)
        pka = list(pka_dict.values())[0]
        pka_type = list(pka_dict.keys())[0]
        return PkaResponse(pka=pka, pka_type=pka_type)


    def _calculate_pka(self, mol:Chem.Mol) -> dict:
        """
        Returns a list of pKa values for the molecule.

        :param mol: The molecule to calculate pKa values for.
        :return: A list of pKa values for the molecule.
        """
        pka_predictor = PkaPredictor()
        pka = pka_predictor.predict([mol])
        return pka[0]