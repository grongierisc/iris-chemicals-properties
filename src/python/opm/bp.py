from rdkit import Chem

from grongier.pex import BusinessProcess

from msg import GenerateSdfRequest, GenerateSdfResponse, SmilesRequest, SmilesResponse, CompareRequest, CompareResponse

class GenerateSdFileProcess(BusinessProcess):
    """

    """
    def on_message(self, request:GenerateSdfRequest) -> GenerateSdfResponse:
        """

        """
        rsp = self.send_request_sync("Python.bordkit.RDKitOperation", request=request)
        # from bo_rdkit import RDKitOperation
        # bo = RDKitOperation()
        # rsp = bo.get_properties_from_smiles(SmilesRequest(smiles=request.smiles))
        
        # create an sdf from the proprties
        mol = Chem.MolFromSmiles(rsp.smiles)
        mol.SetProp("IUPAC_NAME", rsp.properties.iupac_name)
        mol.SetProp("FORMULA", rsp.properties.formula)
        mol.SetProp("MW", str(rsp.properties.mw))
        mol.SetProp("SMILES", rsp.properties.smiles)
        mol.SetProp("CLOGP", str(rsp.properties.clogp))
        mol.SetProp("CLOGD", str(rsp.properties.clogd))
        mol.SetProp("TPSA", str(rsp.properties.tpsa))
        mol.SetProp("PKA", str(rsp.properties.pka))
        mol.SetProp("PKA_TYPE", str(rsp.properties.pka_type))
        mol.SetProp("H_DONOR", str(rsp.properties.h_donor))
        mol.SetProp("H_ACCEPTOR", str(rsp.properties.h_acceptor))
        mol.SetProp("HEAVY_ATOM_COUNT", str(rsp.properties.heavy_atom_count))
        mol.SetProp("ROTATABLE_BONDS", str(rsp.properties.rotatable_bonds))
        mol.SetProp("IMAGE", str(rsp.image))

        w = Chem.SDWriter(request.filename)
        w.write(mol)
        w.close()
        
        return GenerateSdfResponse(
            filename=request.filename
        )

class SmilesProcess(BusinessProcess):
    """

    """
    def get_smiles(self, request:SmilesRequest) -> SmilesResponse:
        """

        """
        rsp = self.send_request_sync("Python.bordkit.RDKitOperation", request=request)
        # from bo_rdkit import RDKitOperation
        # bo = RDKitOperation()
        # rsp = bo.get_properties_from_smiles(SmilesRequest(smiles=request.smiles))

        return rsp

class CompareProcess(BusinessProcess):
    """

    """
    def on_message(self, request:CompareRequest) -> CompareResponse:
        """

        """
        # get the properties from the sdf file
        prop_sdf = self.extract_sdf_properties(request.filename)
        # get the properties from the smiles
        prop_smiles = self.extract_smiles_properties(request.smiles)
        # diff the properties
        diff_prop = self.diff_properties(prop_smiles, prop_sdf)

        return CompareResponse(
            prop_smiles=prop_smiles,
            prop_sdf=prop_sdf,
            diff_prop=diff_prop
        )

    def diff_properties(self, prop_smiles, prop_sdf):
        """
        Diff the properties
        """
        # diff the properties name and value
        diff_prop = {k: (v, prop_sdf[k]) for k, v in prop_smiles.items() if k in prop_sdf and v != prop_sdf[k]}
        # diff the properties name
        diff_prop_name = {k: (v, prop_sdf[k]) for k, v in prop_smiles.items() if k not in prop_sdf}
        diff_prop_name.update({k: (v, prop_sdf[k]) for k, v in prop_sdf.items() if k not in prop_smiles})
        # add the diff properties name
        diff_prop.update(diff_prop_name)

        return diff_prop

    def extract_smiles_properties(self, smiles):
        """
        Extract the properties from the smiles
        """
        # create the process
        process = SmilesProcess()
        # start the process
        msg = SmilesRequest(smiles=smiles)
        rsp = process.on_message(msg)
        # get the properties
        properties = rsp.properties.__dict__
        # change all the properties name in lower case
        properties = {k.lower(): v for k, v in properties.items()}

        return properties

    def extract_sdf_properties(self,filename):
        """
        Extract the properties from the sdf file
        """
        mols = Chem.SDMolSupplier(filename)
        mol = mols[0]
        # get all the properties as a dictionary
        properties = mol.GetPropsAsDict()
        # change all the properties name in lower case
        properties = {k.lower(): v for k, v in properties.items()}

        return properties

if __name__ == "__main__":
    # create the process
    process = CompareProcess()
    # start the process
    filename = '/Users/grongier/git/iris-chemicals-properties/misc/test.sdf'
    msg = CompareRequest(filename=filename, smiles="CC(=O)Nc1ccc(cc1)C(=O)O")
    rsp = process.on_message(msg)
    # print the result
    print(rsp)
    # process = GenerateSdFileProcess()
    # msg = GenerateSdfRequest(smiles="CC(=O)Nc1ccc(cc1)C(=O)O", filename="/Users/grongier/git/iris-chemicals-properties/misc/test.sdf")
    # process.on_message(msg)