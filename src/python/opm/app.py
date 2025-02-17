import codecs

from fastapi import FastAPI, UploadFile
from starlette.datastructures import FormData
from fastapi.responses import Response, FileResponse
from fastapi.middleware.cors import CORSMiddleware

from msg import SmilesResponse,SdfExtractorResponse, Compare2SdfResponse, Compare2SmilesResponse, CompareResponse
from msg import SmilesRequest,SdfExtractorRequest,Compare2SdfRequest, Compare2SmilesRequest, CreateImage2SdfRequest, CreateImage2SmilesRequest, CreateImageRequest, CompareRequest

from iop import Director

director = Director()

app = FastAPI()

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/smiles", response_model=SmilesResponse)
def smiles(request:SmilesRequest):
    bs=director.get_business_service('Python.bs.Rest')
    return bs.on_smiles_process(request)

@app.post("/sdf", response_model=SdfExtractorResponse)
def sdf(file:UploadFile):
    # write the file to disk
    with open(file.filename, "wb") as buffer:
        buffer.write(file.file.read())
    request = SdfExtractorRequest(filename=file.filename)
    bs=director.get_business_service('Python.bs.Rest')
    return bs.on_sdf_process(request)

@app.post("/smiles/img", responses={200: {"content": {"image/png": {}}}})
def smiles_img(request:SmilesRequest):
    bs=director.get_business_service('Python.bs.Rest')
    img_base64 = bs.on_smiles_img(request)
    # convert base64 to binary
    img = codecs.decode(img_base64.encode(), 'base64')
    return Response(img, media_type="image/png")


@app.post("/sdf/img" , responses={200: {"content": {"image/png": {}}}})
def sdf_img(file:UploadFile):
    # write the file to disk
    with open(file.filename, "wb") as buffer:
        buffer.write(file.file.read())
    request = SdfExtractorRequest(filename=file.filename)
    bs=director.get_business_service('Python.bs.Rest')
    img_base64 = bs.on_sdf_img(request)
    # convert base64 to binary
    img = codecs.decode(img_base64.encode(), 'base64')
    return Response(img, media_type="image/png")

@app.post("/compare/smiles", response_model=Compare2SmilesResponse)
def compare_smiles(request:Compare2SmilesRequest):
    bs=director.get_business_service('Python.bs.Rest')
    return bs.on_compare_smiles(request)

@app.post("/compare/smiles/img", responses={200: {"content": {"image/png": {}}}})
def compare_smiles_img(request:CreateImage2SmilesRequest):
    bs=director.get_business_service('Python.bs.Rest')
    img_base64 = bs.on_compare_smiles_img(request)
    # convert base64 to binary
    img = codecs.decode(img_base64.encode(), 'base64')
    return Response(img, media_type="image/png")

@app.post("/compare/sdf/img", responses={200: {"content": {"image/png": {}}}})
def compare_sdf_img(file_a: UploadFile, file_b: UploadFile):
    request = CreateImage2SdfRequest(filename_a=file_a.filename, filename_b=file_b.filename)
    bs=director.get_business_service('Python.bs.Rest')
    img_base64 = bs.on_compare_sdf_img(request)
    # convert base64 to binary
    img = codecs.decode(img_base64.encode(), 'base64')
    return Response(img, media_type="image/png")

@app.post("/compare/sdf", response_model=Compare2SdfResponse)
def compare_sdf(file_a: UploadFile, file_b: UploadFile):
    request = Compare2SdfRequest(filename_a=file_a.filename, filename_b=file_b.filename)
    bs=director.get_business_service('Python.bs.Rest')
    return bs.on_compare_sdf(request)

@app.post("/compare", response_model=CompareResponse)
def compare(smiles: str, file: UploadFile):
    # write the file to disk
    req = CompareRequest(smiles=smiles, filename=file.filename)
    bs=director.get_business_service('Python.bs.Rest')
    return bs.on_compare(req)

@app.post("/compare/img", responses={200: {"content": {"image/png": {}}}})
def compare_img(smiles: str, file: UploadFile):
    # write the file to disk
    req = CreateImageRequest(smiles=smiles, filename=file.filename)
        
    bs=director.get_business_service('Python.bs.Rest')
    img_base64 = bs.on_compare_img(req)
    # convert base64 to binary
    img = codecs.decode(img_base64.encode(), 'base64')
    return Response(img, media_type="image/png")

if __name__ == "__main__":
    import uvicorn
    # bind port 53773
    uvicorn.run(app, host="0.0.0.0", port=53773)
