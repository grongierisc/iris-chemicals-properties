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