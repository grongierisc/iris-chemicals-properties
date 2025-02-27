from grongier.pex import BusinessOperation

from sqlalchemy import create_engine,text

from msg import (
    CreatePersistenceRequest, CreatePersistenceResponse, 
    DeletePersistenceRequest, DeletePersistenceResponse,
    SmilesVectorCosineRequest, AllPersistenceRequest,
    SmilesVectorCosineResponse, AllPersistenceResponse,
    SmilesVectorCosine
)

INIT_EMBEDDING = """
INSERT INTO %Embedding.Config (Name, Configuration, EmbeddingClass, VectorLength, Description)
  VALUES ('opm-transformers-config',
          '{"modelName":"smilestovec"}',
          'Opm.Smiles2Vec',
          167,
          'embedding model')
"""

INIT_EMBEDDING_RANDOM = """
INSERT INTO %Embedding.Config (Name, Configuration, EmbeddingClass, Description)
  VALUES ('sentence-transformers-config',
          '{"modelName":"sentence-transformers/all-MiniLM-L6-v2",
            "hfCachePath":"/usr/irissys/hfCache",
            "maxTokens": 256,
            "checkTokenCount": true}',
          '%Embedding.SentenceTransformers',
          'a small SentenceTransformers embedding model')
"""

CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS
Opm.VectorTable (
  Smiles VARCHAR(250) PRIMARY KEY UNIQUE,
  SmilesEmbedding EMBEDDING('opm-transformers-config','Smiles'),
  RandomEmbedding EMBEDDING('sentence-transformers-config','Smiles')
  )
"""

SELECT_ALL = """
SELECT Smiles,SmilesEmbedding, RandomEmbedding FROM Opm.VectorTable
"""


SELECT_VECTOR = """
SELECT Smiles,SmilesEmbedding, RandomEmbedding FROM Opm.VectorTable where Smiles = :smiles
"""

INSERT_VECTOR = """
INSERT INTO Opm.VectorTable (Smiles) VALUES (:smiles)
"""

DELETE_VECTOR = """
DELETE FROM Opm.VectorTable where Smiles = :smiles
"""

SEARCH_VECTOR = """
SELECT Smiles,SmilesEmbedding, RandomEmbedding, Opm.Smiles2Vec_JACCARDBIT(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')),VECTOR_COSINE(RandomEmbedding,EMBEDDING(:smiles,'sentence-transformers-config')) FROM Opm.VectorTable
ORDER BY Opm.Smiles2Vec_JACCARDBIT(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')) DESC
"""

class Persist(BusinessOperation):

    def on_init(self):
        self.engine = create_engine('iris+emb://IRISAPP', isolation_level="AUTOCOMMIT")
        self.str_to_emb_int = lambda x: [int(i) for i in x.split(",")]
        self.str_to_emb_float = lambda x: [float(i) for i in x.split(",")]

        self.log_info("Creating model for ChemBert embedding", to_console=True)
        # connect with autocommit=True
        with self.engine.connect() as connection:
            try:
                connection.execute(
                    text(INIT_EMBEDDING)
                )
            except Exception as e:
                self.log_warning(f"Error initializing embedding: {e}")

            self.log_info("Creating model for SentenceTransformers embedding", to_console=True)
            try:
                connection.execute(
                    text(INIT_EMBEDDING_RANDOM),
                )
            except Exception as e:
                self.log_warning(f"Error initializing random embedding: {e}")

            self.log_info("Creating table for SMILES persistence", to_console=True)

            connection.execute(
                text(CREATE_TABLE)
            )

            self.log_info("Inserting examples", to_console=True)
            for smiles in [
                'O',
                'CCO',
                'CO',
                'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
                'CC(Cc1ccc(cc1)C(C(=O)O)C)C',
                'Cc1ccccc1',
                'CC1=CC=CC=C1',
                'CC(=O)OC1=CC=CC=C1C(=O)O',
                'CC(=O)NC1=CC=C(C=C1)O',
                'CC(=O)O[C@H]1C=C[C@H]2[C@H]3CC4=C5[C@]2([C@H]1OC5=C(C=C4)OC(=O)C)CCN3C',
                'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O'
            ]:
                try:
                    self.on_create(CreatePersistenceRequest(smiles=smiles))
                except Exception as e:
                    self.log_warning(f"Error inserting examples:{e}")


    def on_create(self, msg:CreatePersistenceRequest):
        with self.engine.connect() as connection:
            new_row = connection.execute(
                text(SELECT_VECTOR),
                {"smiles":msg.smiles}
            ).fetchone()
            if new_row is None:
                connection.execute(
                    text(INSERT_VECTOR),
                    {"smiles":msg.smiles}
                )
                new_row = connection.execute(
                    text(SELECT_VECTOR),
                    {"smiles":msg.smiles}
                ).fetchone()
            else:
                new_row = (msg.smiles, new_row[1], new_row[2])
            connection.commit()

        return CreatePersistenceResponse(smiles=new_row[0], embedding=self.str_to_emb_int(new_row[1]), embedding_random=self.str_to_emb_float(new_row[2]))
    
    def on_delete(self, msg:DeletePersistenceRequest):
        with self.engine.connect() as connection:
            connection.execute(
                text(DELETE_VECTOR),
                {"smiles":msg.smiles}
            )
            connection.commit()
        return DeletePersistenceResponse(smiles=msg.smiles)
    
    def on_smiles_vector_cosine(self, msg:SmilesVectorCosineRequest):
        with self.engine.connect() as connection:
            result = connection.execute(
                text(SEARCH_VECTOR),
                {"smiles":msg.smiles}
            ).fetchall()
        return SmilesVectorCosineResponse(
            result=[SmilesVectorCosine(
                smiles=row[0],
                embedding=self.str_to_emb_int(row[1]),
                embedding_random=self.str_to_emb_float(row[2]),
                cosine=float(row[3]),
                cosine_random=float(row[4])
            ) for row in result]
        )


    def on_all(self, msg:AllPersistenceRequest):
        with self.engine.connect() as connection:
            all_rows = connection.execute(
                text(SELECT_ALL)
            ).fetchall()
        return AllPersistenceResponse(
            result=[SmilesVectorCosine(smiles=row[0], embedding=self.str_to_emb_int(row[1]), embedding_random=self.str_to_emb_float(row[2]) ) for row in all_rows]
        )
        
