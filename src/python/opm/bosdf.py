from grongier.pex import BusinessOperation

from msg import SdfExtractorRequest, SdfExtractorResponse, CreateSdfRequest, CreateSdfResponse, CreateImageRequest

from rdkit import Chem

class SdfOperation(BusinessOperation):

    def extract_sdf(self, request: SdfExtractorRequest) -> SdfExtractorResponse:
        """
        Extract the properties from the sdf file
        """
        # get the properties from the sdf file
        prop_sdf = self.extract_sdf_properties(request.filename)

        #get the drawing from the sdf file
        self.send_request_sync("Python.bomisc.GenerateImageOperation", CreateImageRequest(smiles=None, filename=request.filename))

        return SdfExtractorResponse(properties=prop_sdf)

    def extract_sdf_properties(self, filename):
        """
        Extract the properties from the sdf file
        """
        # get the properties from the sdf file
        mols = Chem.SDMolSupplier(filename)
        mol = mols[0]
        # get all the properties as a dictionary
        properties = mol.GetPropsAsDict()
        # change all the properties name in lower case
        properties = {k.lower(): v for k, v in properties.items()}

        return properties

    def create_sdf(self, request: CreateSdfRequest) -> CreateSdfRequest:
        """
        Create the sdf file
        """
        # create an sdf from the proprties
        mol = Chem.MolFromSmiles(request.properties.smiles)
        # for each property not null, add it to the mol as str
        for k, v in request.properties.__dict__.items():
            if v:
                mol.SetProp(k.upper(), str(v))

        w = Chem.SDWriter(request.filename)
        w.write(mol)
        w.close()

        return CreateSdfResponse(filename=request.filename)
