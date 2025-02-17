msg.SmilesRequest
{
"smiles":"CC(=O)Nc1ccc(cc1)C(=O)O"
}

msg.GenerateSdfRequest
{
"smiles":"CC(=O)Nc1ccc(cc1)C(=O)O",
"filename":"/irisdev/app/misc/test4.sdf"
}

msg.SdfExtractorRequest
{
"filename":"/irisdev/app/misc/test4.sdf"
}

msg.CompareRequest
{
"smiles":"C(=O)O",
"filename":"/irisdev/app/misc/test.sdf"
}

msg.CreateImageRequest
{
"smiles":"CC(=O)Nc1ccc(cc1)C(=O)O",
"filename":"/irisdev/app/misc/test.png"
}


iop --test Python.bp.CompareProcess --classname msg.CompareRequest --body '{ "smiles":"C(=O)O", "filename":"/irisdev/app/misc/test.sdf" }'