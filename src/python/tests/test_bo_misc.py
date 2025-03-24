from opm import bomisc
from opm import msg

from iop._serialization import (
    SerializationError,
    serialize_message,
    deserialize_message,
    serialize_pickle_message,
    deserialize_pickle_message,
)

class TestBoMisc:

    def test_iupac(self):
        bo = bomisc.IUPACOperation()
        bo.on_init()
        rsp = bo.process(msg.SmilesRequest(smiles='CC1=CC=CC=C1'))
        ser = serialize_message(rsp)
        des = deserialize_message(ser)
        assert des.properties.iupac_name == rsp.properties['iupac_name']