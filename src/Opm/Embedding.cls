Class Opm.Embedding Extends %Embedding.Interface
{

ClassMethod GetPythonModule() As %SYS.Python
{
	Do ##class(Opm.Embedding).SetPythonPath("/irisdev/app/src/python/opm")
	set importlib = ##class(%SYS.Python).Import("importlib")
	quit importlib."import_module"("embedding")
}

ClassMethod SetPythonPath(pClasspaths)
{
    set sys = ##class(%SYS.Python).Import("sys")

    for i=0:1:(sys.path."__len__"()-1) {
        Try {
            if sys.path."__getitem__"(i) = pClasspaths {
                do sys.path."__delitem__"(i)
            }
        }
        Catch ex {
            // do nothing
        }

    }
    do sys.path.insert(0, pClasspaths)
}

ClassMethod Embedding(
	input,
	configuration) As %Vector
{
	set json = ##class(%SYS.Python).Import("json")
	set configuration = json."loads"(configuration)
	set module = ##class(Opm.Embedding).GetPythonModule()
	set embeddings = module."get_embedding"(input, configuration)
	return ##class(%Library.Vector).DisplayToLogical(embeddings)
	#; // Convert from dynamic array to vector
	#; set iter = embeddings.%GetIterator()
	#; while iter.%GetNext(.key,.value, .type){
	#; 	set $vector(v, $increment(i), "xf32") = value
	#; }

	#; return v
}

ClassMethod IsValidConfig(
	config As %DynamicObject,
	ByRef errorMsg As %String) As %Boolean
{
	set json = ##class(%SYS.Python).Import("json")
	set configuration = json."loads"(config.%ToJSON())
	set module = ##class(Opm.Embedding).GetPythonModule()
	quit module."is_valid_config"(configuration)
}

}
