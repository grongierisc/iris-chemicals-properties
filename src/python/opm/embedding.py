
def is_valid_config(config:dict)->bool:
    if config.get("modelName") is None:
        raise ValueError("modelName is required in the config")
    if config.get("hfCachePath") is None:
        raise ValueError("hfCachePath is required in the config")
    try:
        import transformers
    except ImportError:
        raise ImportError("transformers is not installed")
    transformers.RobertaModel.from_pretrained(config.get("modelName"), num_labels=2, 
                                 add_pooling_layer=True, 
                                 cache_dir=config.get("hfCachePath"))
    return True

def get_embedding(smiles:str,config)->list:
    import transformers
    pretrained_path = config.get("modelName")
    tokenizer = transformers.AutoTokenizer.from_pretrained(pretrained_path)
    model = transformers.RobertaModel.from_pretrained(pretrained_path, num_labels=2, 
                                 add_pooling_layer=True, 
                                 cache_dir=config.get("hfCachePath"))
    inputs = tokenizer.encode_plus(smiles, return_tensors="pt")
    outputs = model.forward(input_ids=inputs["input_ids"])
    embeddings = outputs.last_hidden_state
    return str(embeddings.tolist()[0][0])
