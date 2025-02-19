Class Opm.Embedding Extends %Embedding.Interface
{

ClassMethod Embedding(
	input,
	configuration) As %Vector
{
	set config = {}.%FromJSON(configuration)
	set transformers = ##class(%SYS.Python).Import("transformers")

	set args = { "pretrained_model_name_or_path": (config.%Get("modelName")), "cache_dir": (config.%Get("hfCachePath")) }
	set model = transformers.RobertaModel."from_pretrained"(args...)
	set tokenizer = transformers.RobertaTokenizer."from_pretrained"(config.%Get("modelName"))

	set args = { "add_special_tokens": 1, "return_tensors": "pt" }
	set inputs = tokenizer."encode_plus"(input, args...)

	set args = { "input_ids": (inputs."input_ids"), "attention_mask": (inputs."attention_mask") }
	set embeddings = model.forward(args...)

	set builtins = ##class(%SYS.Python).Import("builtins")

	set strEmb = builtins.str(embeddings."last_hidden_state".tolist()."__getitem__"(0)."__getitem__"(0))

	set embeddings = [].%FromJSON(strEmb)

	// Convert from dynamic array to vector
	set iter = embeddings.%GetIterator()
	while iter.%GetNext(.key,.value, .type){
		set $vector(v, $increment(i), "xf32") = value
	}

	return v
}

/// Description
ClassMethod Test() As %Status
{
	Set sc = $$$OK
	set config = {
		"modelName": "seyonec/PubChem10M_SMILES_BPE_450k",
		"hfCachePath": "/usr/irissys/hfCache"
	}
	set configuration = {
		"config": (config)
	}
	set input = "CCC"
	set embeddings = ..Embedding(input, configuration)
	write embeddings
	Return sc
}

ClassMethod IsValidConfig(
	config As %DynamicObject,
	ByRef errorMsg As %String) As %Boolean
{
	if config.%Get("modelName") = "" {
        set errorMsg = $$$Text("'modelName' not set", "%SQL.VECTOR")
        return 0
    }
    if config.%Get("hfCachePath") = "" {
        set errorMsg = $$$Text("'hfCachePath' not set", "%SQL.VECTOR")
        return 0
    }
	try {
        do ..CheckInstall(config.%Get("pythonPath", ""))
    } catch e {
        set errorMsg = e.Data
        return 0
    }
    try {
        do ..DownloadModel(config.%Get("modelName"), config.%Get("hfCachePath"), config.%Get("hfToken",""),config.%Get("pythonPath", ""))
    } catch e {
        set errorMsg = $$$FormatText($$$Text("%1. Error downloading model","%SQL.VECTOR"),e.Data)
        return 0
    }
	return 1
}

ClassMethod CheckInstall(pythonPath As %String)
{
	set transformers = ##class(%SYS.Python).Import("transformers")
	quit
}

ClassMethod DownloadModel(
	modelName As %String,
	hfCachePath As %String,
	hfToken As %String,
	pythonPath As %String)
{
	set transformers = ##class(%SYS.Python).Import("transformers")
	//model = RobertaModel.from_pretrained(pretrained_path, num_labels=2, add_pooling_layer=True, cache_dir="hfCache")
	set args = { "pretrained_model_name_or_path": (modelName), "cache_dir": (hfCachePath) }
	set model = transformers.RobertaModel."from_pretrained"(args...)
	quit
}

}
