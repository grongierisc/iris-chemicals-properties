from grongier.pex import BusinessOperation

from msg import SdfExtractorRequest, SdfExtractorResponse, CreateSdfRequest, CreateSdfResponse

from rdkit import Chem

class SdfOperation(BusinessOperation):

    def extract_sdf(self, request: SdfExtractorRequest) -> SdfExtractorResponse:
        """
        Extract the properties from the sdf file
        """
        # get the properties from the sdf file
        prop_sdf = self.extract_sdf_properties(request.filename)

        return SdfExtractorResponse(prop_sdf)

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
        mol.SetProp("IUPAC_NAME", request.properties.iupac_name)
        mol.SetProp("FORMULA", request.properties.formula)
        mol.SetProp("MW", str(request.properties.mw))
        mol.SetProp("SMILES", request.properties.smiles)
        mol.SetProp("CLOGP", str(request.properties.clogp))
        mol.SetProp("CLOGD", str(request.properties.clogd))
        mol.SetProp("TPSA", str(request.properties.tpsa))
        mol.SetProp("PKA", str(request.properties.pka))
        mol.SetProp("PKA_TYPE", str(request.properties.pka_type))
        mol.SetProp("H_DONOR", str(request.properties.h_donor))
        mol.SetProp("H_ACCEPTOR", str(request.properties.h_acceptor))
        mol.SetProp("HEAVY_ATOM_COUNT", str(request.properties.heavy_atom_count))
        mol.SetProp("ROTATABLE_BONDS", str(request.properties.rotatable_bonds))

        w = Chem.SDWriter(request.filename)
        w.write(mol)
        w.close()

        return CreateSdfResponse(filename=request.filename)
