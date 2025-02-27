from bomisc import IUPACOperation, GenerateImageOperation
from bp import CompareProcess, SdfProcess, SmilesProcess
from bopka import PkaPredictorOperation
from bordkit import RDKitOperation
from bosdf import SdfOperation
from bopersist import Persist
from bs import Rest

# preload models from persistence
bo = Persist()
bo.log_info("Preloading models...", to_console=True)
try:
    bo.on_init()
except Exception as e:
    bo.log_error(f"Error preloading models: {e}", to_console=True)
bo.log_info("Models preloaded", to_console=True)

CLASSES = {
    'Python.bomisc.IUPACOperation': IUPACOperation,
    'Python.bomisc.GenerateImageOperation': GenerateImageOperation,
    'Python.bp.CompareProcess': CompareProcess,
    'Python.bp.SdfProcess': SdfProcess,
    'Python.bp.SmilesProcess': SmilesProcess,
    'Python.bopka.PkaPredictorOperation': PkaPredictorOperation,
    'Python.bordkit.RDKitOperation': RDKitOperation,
    'Python.bosdf.SdfOperation': SdfOperation,
    'Python.bopersist.Persist': Persist,
    'Python.bs.Rest': Rest
}

PRODUCTIONS = [
    {
    "Opm.Production": {
        "@Name": "Opm.Production",
        "@TestingEnabled": "true",
        "@LogGeneralTraceEvents": "true",
        "Description": "",
        "ActorPoolSize": "2",
        "Item": [
            {
                "@Name": "Python.bs.Rest",
                "@Category": "",
                "@ClassName": "Python.bs.Rest",
            },
            {
                "@Name": "Python.bp.CompareProcess",
                "@Category": "",
                "@ClassName": "Python.bp.CompareProcess",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bp.SdfProcess",
                "@Category": "",
                "@ClassName": "Python.bp.SdfProcess",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bp.SmilesProcess",
                "@Category": "",
                "@ClassName": "Python.bp.SmilesProcess",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bopka.PkaPredictorOperation",
                "@Category": "",
                "@ClassName": "Python.bopka.PkaPredictorOperation",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bordkit.RDKitOperation",
                "@Category": "",
                "@ClassName": "Python.bordkit.RDKitOperation",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bosdf.SdfOperation",
                "@Category": "",
                "@ClassName": "Python.bosdf.SdfOperation",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bomisc.IUPACOperation",
                "@Category": "",
                "@ClassName": "Python.bomisc.IUPACOperation",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bomisc.GenerateImageOperation",
                "@Category": "",
                "@ClassName": "Python.bomisc.GenerateImageOperation",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            },
            {
                "@Name": "Python.bopersist.Persist",
                "@Category": "",
                "@ClassName": "Python.bopersist.Persist",
                "@PoolSize": "1",
                "@Enabled": "true",
                "@Foreground": "false",
                "@Comment": "",
                "@LogTraceEvents": "false",
                "@Schedule": ""
            }
        ]
    }
}
]