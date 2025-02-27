from decimal import Decimal
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
class Compare2SmilesRequest(Message):
    smiles_a:str = None
    smiles_b:str = None

@dataclass
class Compare2SmilesResponse(Message):
    prop_a:MolProperties = None
    prop_b:MolProperties = None
    diff_prop:MolProperties = None

@dataclass
class Compare2SdfRequest(Message):
    filename_a:str = None
    filename_b:str = None

@dataclass
class Compare2SdfResponse(Message):
    prop_a:MolProperties = None
    prop_b:MolProperties = None
    diff_prop:MolProperties = None

@dataclass
class CompareResponse(Message):
    prop_smiles:MolProperties = None
    prop_sdf:MolProperties = None
    diff_prop:MolProperties = None

@dataclass
class SdfExtractorRequest(Message):
    filename:str = None

@dataclass
class SdfExtractorResponse(Message):
    properties:MolProperties = None

@dataclass
class CreateSdfRequest(Message):
    properties:MolProperties = None
    filename:str = None

@dataclass
class CreateSdfResponse(Message):
    filename:str = None

@dataclass
class CreateImageRequest(Message):
    smiles:str = None
    filename:str = None

@dataclass
class CreateImage2SmilesRequest(Message):
    smiles_a:str = None
    smiles_b:str = None

@dataclass
class CreateImage2SdfRequest(Message):
    filename_a:str = None
    filename_b:str = None

@dataclass
class CreatePersistenceRequest(Message):
    smiles:str = None

@dataclass
class CreatePersistenceResponse(Message):
    smiles:str = None
    embedding:list[int] = None
    embedding_random:list[float] = None

@dataclass
class DeletePersistenceRequest(Message):
    smiles:str = None

@dataclass
class DeletePersistenceResponse(Message):
    smiles:str = None

@dataclass
class SmilesVectorCosineRequest(Message):
    smiles:str = None

@dataclass
class SmilesVectorCosine:
    smiles:str = None
    embedding:list[int] = None
    embedding_random:list[float] = None
    cosine:float = None
    cosine_random:float = None

@dataclass
class SmilesVectorCosineResponse(Message):
    result: list[SmilesVectorCosine] = None

@dataclass
class AllPersistenceRequest(Message):
    pass

@dataclass
class AllPersistenceResponse(Message):
    result: list[CreatePersistenceResponse] = None