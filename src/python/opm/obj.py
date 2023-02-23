from dataclasses import dataclass

@dataclass
class MolProperties:
    iupac_name: str = None
    formula: str = None
    mw: float = None
    smiles: str = None
    clogp: float = None
    clogd: float = None
    tpsa: float = None
    pka: float = None
    pka_type: str = None
    h_donor: int = None
    h_acceptor: int = None
    heavy_atom_count: int = None
    rotatable_bonds: int = None