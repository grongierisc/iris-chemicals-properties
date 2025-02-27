Class Opm.Smiles2Vec Extends %Embedding.Interface
{

ClassMethod GetPythonModule() As %SYS.Python
{
	Do ..SetPythonPath("/irisdev/app/src/python/opm")
	set importlib = ##class(%SYS.Python).Import("importlib")
	quit importlib."import_module"("smiles2vec")
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

ClassMethod Test() As %String
{
    set smiles = "CCO"
    set configuration = {"smiles": "smiles"}.%ToJSON()
    set embeddings = ..Embedding(smiles, configuration)
    return embeddings
}

ClassMethod Embedding(
	input,
	configuration) As %Vector
{
	set json = ##class(%SYS.Python).Import("json")
	set conf = json."loads"(configuration)
	set module = ..GetPythonModule()

	set embeddings = module."get_embedding"(input, conf)
    set strippedInput = $zstrip(embeddings, "*W","[]")
    set vec = $vop("fromstring", strippedInput, "integer", ",", 3)
    return vec
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
	set module = ..GetPythonModule()
	quit module."is_valid_config"(configuration)
}

ClassMethod JACCARD(
	vec1 As %Library.Vector = "",
	vec2 As %Library.Vector = "") As %Library.Double [ SqlProc ]
{
        set return=""
        set v1len=$vop("length",vec1),v2len=$vop("length",vec2)
        throw:v1len'=v2len ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARD on vectors of different lengths "_v1len_" and "_v2len)
        set v1type=$vop("type",vec1),v2type=$vop("type",vec2)
        throw:",integer,double,decimal,xf32,"'[(","_v1type_",") ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARD because Argument #1 is of type "_v1type_", not a numeric datatype")
        throw:",integer,double,decimal,xf32,"'[(","_v2type_",") ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARD because Argument #2 is of type "_v2type_", not a numeric datatype")
        throw:v1type'=v2type ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARD on vectors of different datatypes "_v1type_" and "_v2type)
        try {
            set intersection=0,union=0
            for i=1:1:v1len {
                set v1=$VECTOR(vec1,i)
                for j=1:1:v2len {
                    set v2=$VECTOR(vec2,j)
                    if v1=v2 { set intersection=intersection+1 }
                }
                set union=v1len+v2len-intersection
            }
            set return=intersection/union
        }
        catch ex {
            if ex.Name="<ILLEGAL VALUE>" { throw ##class(%Exception.General).%New(,,,"Overflow error occured on call to $SYSTEM.SQL.Functions.JACCARD")
            } else { throw ##class(%Exception.General).%New(,,,"Unexpected error occurred on call to $SYSTEM.SQL.Functions.JACCARD")}
        }
        RETURN return
}

ClassMethod JACCARDBIT(
	vec1 As %Library.Vector = "",
	vec2 As %Library.Vector = "") As %Library.Double [ SqlProc ]
{
        set return=""
        set v1len=$vop("length",vec1),v2len=$vop("length",vec2)
        throw:v1len'=v2len ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARDBIT on vectors of different lengths "_v1len_" and "_v2len)
        set v1type=$vop("type",vec1),v2type=$vop("type",vec2)
        throw:",integer,"'[(","_v1type_",") ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARDBIT because Argument #1 is of type "_v1type_", not a numeric datatype")
        throw:",integer,"'[(","_v2type_",") ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARDBIT because Argument #2 is of type "_v2type_", not a numeric datatype")
        throw:v1type'=v2type ##class(%Exception.General).%New(,,,"Cannot perform vector function JACCARDBIT on vectors of different datatypes "_v1type_" and "_v2type)
        try {
            set M11=0,M10=0,M01=0,M00=0
            for i=1:1:v1len {
                set v1=$VECTOR(vec1,i)
                set v2=$VECTOR(vec2,i)
                if ((v1=1) && (v2=1)) { set M11=M11+1 }
                if ((v1=1) && (v2=0)) { set M10=M10+1 }
                if ((v1=0) && (v2=1)) { set M01=M01+1 }
                if ((v1=0) && (v2=0)) { set M00=M00+1 }
            }
            set return=M11/(M11+M10+M01)
        }
        catch ex {
            if ex.Name="<ILLEGAL VALUE>" { throw ##class(%Exception.General).%New(,,,"Overflow error occured on call to $SYSTEM.SQL.Functions.JACCARDBIT")
            } else { throw ##class(%Exception.General).%New(,,,"Unexpected error occurred on call to $SYSTEM.SQL.Functions.JACCARDBIT")}
        }
        RETURN return
}

}
