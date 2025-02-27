import codecs
import os

from fastapi import FastAPI, UploadFile
from fastapi.responses import Response, JSONResponse
from fastapi.middleware.cors import CORSMiddleware

from msg import (
    SmilesResponse, SdfExtractorResponse, Compare2SdfResponse, 
    Compare2SmilesResponse, CompareResponse, SmilesRequest,
    SdfExtractorRequest, Compare2SdfRequest, Compare2SmilesRequest, 
    CreateImage2SdfRequest, CreateImage2SmilesRequest, CreateImageRequest, 
    CompareRequest,
    CreatePersistenceRequest, CreatePersistenceResponse, 
    DeletePersistenceRequest, DeletePersistenceResponse,
    SmilesVectorCosineRequest, AllPersistenceRequest,
    SmilesVectorCosineResponse, AllPersistenceResponse,
)
from iop import Director

WORKING_DIR = "/usr/irissys/mgr/IRISAPP_DATA/"
director = Director()
app = FastAPI()

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    expose_headers=["X-Session-Id"]
)

def get_business_service():
    bs = director.get_business_service('Python.bs.Rest', True)
    return bs, bs.iris_handle._SessionId

def save_uploaded_file(file: UploadFile) -> str:
    filepath = os.path.join(WORKING_DIR, file.filename)
    with open(filepath, "wb") as buffer:
        buffer.write(file.file.read())
    return filepath

def create_image_response(img_base64: str, session_id: str) -> Response:
    img = codecs.decode(img_base64.encode(), 'base64')
    return Response(img, media_type="image/png", headers={"X-Session-Id": session_id})

# SMILES endpoints
@app.post("/smiles", response_model=SmilesResponse)
def smiles(request: SmilesRequest, response: Response):
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_smiles_process(request)

@app.post("/smiles/img", responses={200: {"content": {"image/png": {}}}})
def smiles_img(request: SmilesRequest):
    bs, session_id = get_business_service()
    img_base64 = bs.on_smiles_img(request)
    return create_image_response(img_base64, session_id)

# SDF endpoints
@app.post("/sdf", response_model=SdfExtractorResponse)
def sdf(file: UploadFile, response: Response):
    filepath = save_uploaded_file(file)
    request = SdfExtractorRequest(filename=filepath)
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_sdf_process(request)

@app.post("/sdf/img", responses={200: {"content": {"image/png": {}}}})
def sdf_img(file: UploadFile):
    filepath = save_uploaded_file(file)
    request = SdfExtractorRequest(filename=filepath)
    bs, session_id = get_business_service()
    img_base64 = bs.on_sdf_img(request)
    return create_image_response(img_base64, session_id)

# Compare endpoints
@app.post("/compare/smiles", response_model=Compare2SmilesResponse)
def compare_smiles(request: Compare2SmilesRequest, response: Response):
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_compare_smiles(request)

@app.post("/compare/smiles/img", responses={200: {"content": {"image/png": {}}}})
def compare_smiles_img(request: CreateImage2SmilesRequest):
    bs, session_id = get_business_service()
    img_base64 = bs.on_compare_smiles_img(request)
    return create_image_response(img_base64, session_id)

@app.post("/compare/sdf", response_model=Compare2SdfResponse)
def compare_sdf(file_a: UploadFile, file_b: UploadFile, response: Response):
    filepath_a = save_uploaded_file(file_a)
    filepath_b = save_uploaded_file(file_b)
    request = Compare2SdfRequest(filename_a=filepath_a, filename_b=filepath_b)
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_compare_sdf(request)

@app.post("/compare/sdf/img", responses={200: {"content": {"image/png": {}}}})
def compare_sdf_img(file_a: UploadFile, file_b: UploadFile):
    filepath_a = save_uploaded_file(file_a)
    filepath_b = save_uploaded_file(file_b)
    request = CreateImage2SdfRequest(filename_a=filepath_a, filename_b=filepath_b)
    bs, session_id = get_business_service()
    img_base64 = bs.on_compare_sdf_img(request)
    return create_image_response(img_base64, session_id)

@app.post("/compare", response_model=CompareResponse)
def compare(smiles: str, file: UploadFile, response: Response):
    filepath = save_uploaded_file(file)
    request = CompareRequest(smiles=smiles, filename=filepath)
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_compare(request)

@app.post("/compare/img", responses={200: {"content": {"image/png": {}}}})
def compare_img(smiles: str, file: UploadFile):
    filepath = save_uploaded_file(file)
    request = CreateImageRequest(smiles=smiles, filename=filepath)
    bs, session_id = get_business_service()
    img_base64 = bs.on_compare_img(request)
    return create_image_response(img_base64, session_id)

# Persistence endpoints
@app.post("/persistence", response_model=CreatePersistenceResponse)
def create_persistence(request: CreatePersistenceRequest, response: Response):
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_query(request)

@app.delete("/persistence", response_model=DeletePersistenceResponse)
def delete_persistence(request: DeletePersistenceRequest, response: Response):
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_query(request)

@app.post("/persistence/cosine", response_model=SmilesVectorCosineResponse)
def smiles_vector_cosine(request: SmilesVectorCosineRequest, response: Response):
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_query(request)

@app.get("/persistence/all", response_model=AllPersistenceResponse)
def all_persistence(response: Response):
    bs, session_id = get_business_service()
    response.headers["X-Session-Id"] = session_id
    return bs.on_query(AllPersistenceRequest())

@app.get("/iupac/{smiles}/iupac_name")
def iupac_name(smiles: str):
    msg = SmilesRequest(smiles=smiles)
    bs, _ = get_business_service()
    return bs.on_iupac(msg)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=53773)
