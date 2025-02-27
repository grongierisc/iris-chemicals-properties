from iop import BusinessService
from msg import SmilesRequest, CreateImageRequest

import codecs

class Rest(BusinessService):
    def on_smiles_process(self, request):
        return self.send_request_sync('Python.bp.SmilesProcess', request)
    
    def on_sdf_process(self, request):
        return self.send_request_sync('Python.bp.SdfProcess', request)
    
    def on_smiles_img(self, request):
        msg = CreateImageRequest(smiles=request.smiles)
        rsp = self.send_request_sync('Python.bomisc.GenerateImageOperation', msg)
        img_stream = ''
        while not rsp.Base64Image.AtEnd:
            img_stream += rsp.Base64Image.Read(1024)
        
        return img_stream
    
    def on_sdf_img(self, request):
        msg = CreateImageRequest(filename=request.filename)
        rsp = self.send_request_sync('Python.bomisc.GenerateImageOperation', msg)
        img_stream = ''
        while not rsp.Base64Image.AtEnd:
            img_stream += rsp.Base64Image.Read(1024)
        
        return img_stream
    
    def on_compare_smiles(self, request):
        return self.send_request_sync('Python.bp.CompareProcess', request)

    def on_compare_sdf(self, request):
        return self.send_request_sync('Python.bp.CompareProcess', request)
    
    def on_compare(self, request):
        return self.send_request_sync('Python.bp.CompareProcess', request)
    
    def on_compare_smiles_img(self, msg):
        rsp = self.send_request_sync('Python.bomisc.GenerateImageOperation', msg)
        img_stream = ''
        while not rsp.Base64Image.AtEnd:
            img_stream += rsp.Base64Image.Read(1024)
        
        return img_stream
    
    def on_compare_sdf_img(self, msg):
        rsp = self.send_request_sync('Python.bomisc.GenerateImageOperation', msg)
        img_stream = ''
        while not rsp.Base64Image.AtEnd:
            img_stream += rsp.Base64Image.Read(1024)
        
        return img_stream
    
    def on_compare_img(self, msg):
        rsp = self.send_request_sync('Python.bomisc.GenerateImageOperation', msg)
        img_stream = ''
        while not rsp.Base64Image.AtEnd:
            img_stream += rsp.Base64Image.Read(1024)
        
        return img_stream
    
    def on_query(self, msg):
        return self.send_request_sync('Python.bopersist.Persist', msg)
    
    def on_iupac(self, msg):
        rsp= self.send_request_sync('Python.bomisc.IUPACOperation', msg)
        return rsp.properties.iupac_name
