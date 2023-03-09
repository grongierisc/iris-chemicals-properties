from grongier.pex import BusinessOperation

import pandas as pd

from sqlalchemy import create_engine,types

from msg import CreatePersistenceRequest
from obj import MolProperties

import iris

class Persist(BusinessOperation):
    def on_init(self):
        pass
        

    def on_message(self, msg:CreatePersistenceRequest):
        # check is msg is a CreatePersistenceRequest
        if not isinstance(msg, CreatePersistenceRequest):
            raise Exception('msg is not a CreatePersistenceRequest')
        if not isinstance(msg.properties, MolProperties):
            raise Exception('msg.properties is not a MolProperties')
        # cast msg to a dict
        dict = msg.__dict__
        # flatten msg['properties'] for each property defined in MolProperties
        for key, value in msg.properties.__dict__.items():
            if key in MolProperties.__annotations__:
                dict[key] = value
        del dict['properties']

        df = pd.DataFrame([dict])
        db = create_engine('iris+emb:///',execution_options={
            "isolation_level": "AUTOCOMMIT",
        })

        df.to_sql('Mol', db, if_exists='append', dtype= {
            'smiles': types.VARCHAR(length=255),
            'filename': types.VARCHAR(length=255),
            'formula': types.VARCHAR(length=255),
            'clogd': types.Float,
            'clogp': types.Float,
            'h_acceptor': types.Integer,
            'h_donor': types.Integer,
            'heavy_atom_count': types.Integer,
            'mw': types.Float,
            'iupac_name': types.VARCHAR(length=255),
            'pka': types.Float,
            'pka_type': types.VARCHAR(length=255),
            'rotatable_bonds': types.Integer,
            'tpsa': types.Float,
        })
        db.dispose()
        return iris.cls('Ens.Response')._New()

    def on_tear_down(self):
        pass

if __name__ == '__main__':
    bo = Persist()
    bo.on_init()
    msg = CreatePersistenceRequest()
    msg.smiles = 'C1=CC=CC=C1'
    msg.filename = 'test.sdf'
    msg.properties = MolProperties(
        smiles='C1=CC=CC=C1',
        formula='C6H6',
        clogd=0.0,
        clogp=0.0,
        h_acceptor=0,
        h_donor=0,
        heavy_atom_count=6,
        mw=78.11,
        iupac_name='benzene',
        pka=6.0,
        pka_type='pKa',
        rotatable_bonds=0,
        tpsa=0.0,
        )
    bo.on_message(msg)
    print('Done')
