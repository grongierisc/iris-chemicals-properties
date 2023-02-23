from grongier.pex import Message
from dataclasses import dataclass
from obj import MolProperties

@dataclass
class SmilesRequest(Message):
    smiles:str = None

@dataclass
class SmilesResponse(Message):
    image:bytes = None
    smiles:str = None
    properties:MolProperties = None

@dataclass
class PkaRequest(Message):
    smiles:str = None

@dataclass
class PkaResponse(Message):
    pka:float = None
    pka_type:str = None

@dataclass
class GenerateSdfRequest(Message):
    smiles:str = None
    filename:str = None

@dataclass
class GenerateSdfResponse(Message):
    filename:str = None

@dataclass
class CompareRequest(Message):
    smiles:str = None
    filename:str = None

@dataclass
class CompareResponse(Message):
    prop_smiles:MolProperties = None
    prop_sdf:MolProperties = None
    diff_prop:MolProperties = None