import opm.bopka as bopka
from rdkit import Chem

class TestBoPka:

    bo = bopka.PkaPredictorOperation()

    def test_bo_pka(self):
        mol = Chem.MolFromSmiles('CC(=O)O')
        rsp = self.bo._calculate_pka(mol)
        assert float(rsp['acidic']) == 4.072120666503906