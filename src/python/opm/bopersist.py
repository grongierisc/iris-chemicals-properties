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
          '{"modelName":"seyonec/PubChem10M_SMILES_BPE_450k",
          "hfCachePath":"/usr/irissys/hfCache"}',
          'Opm.Embedding',
          768,
          'embedding model')
"""

CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS
Opm.VectorTable (
  Smiles VARCHAR(250) PRIMARY KEY UNIQUE,
  SmilesEmbedding EMBEDDING('opm-transformers-config','Smiles')
  )
"""

SELECT_VECTOR = """
SELECT Smiles, SmilesEmbedding FROM Opm.VectorTable where Smiles = :smiles
"""

INSERT_VECTOR = """
INSERT INTO Opm.VectorTable (Smiles) VALUES (:smiles)
"""

DELETE_VECTOR = """
DELETE FROM Opm.VectorTable where Smiles = :smiles
"""

SEARCH_VECTOR = """
SELECT TOP 5 Smiles, SmilesEmbedding, VECTOR_COSINE(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')) FROM Opm.VectorTable
ORDER BY VECTOR_COSINE(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')) DESC
"""

class Persist(BusinessOperation):

    def on_init(self):
        self.engine = create_engine('iris+emb://IRISAPP')
        self.str_to_emb = lambda x: [float(i) for i in x.split(",")]
        with self.engine.connect() as connection:
            try:
                connection.execute(
                    text(INIT_EMBEDDING)
                )
            except Exception as e:
                self.log_warning(f"Error initializing embedding: {e}")
            connection.execute(
                text(CREATE_TABLE)
            )
            connection.commit()

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
            connection.commit()

        return CreatePersistenceResponse(smiles=new_row[0], embedding=self.str_to_emb(new_row[1]))
    
    def on_delete(self, msg:DeletePersistenceRequest):
        with self.engine.connect() as connection:
            connection.execute(
                text(DELETE_VECTOR),
                {"smiles":msg.smiles}
            )
        return DeletePersistenceResponse(smiles=msg.smiles, embedding=None)
    
    def on_smiles_vector_cosine(self, msg:SmilesVectorCosineRequest):
        with self.engine.connect() as connection:
            result = connection.execute(
                text(SEARCH_VECTOR),
                {"smiles":msg.smiles}
            ).fetchall()
        return SmilesVectorCosineResponse(
            result=[SmilesVectorCosine(smiles=row[0], embedding=self.str_to_emb(row[1]), cosine=float(row[2])) for row in result]
        )


    def on_all(self, msg:AllPersistenceRequest):
        with self.engine.connect() as connection:
            all_rows = connection.execute(
                text("SELECT Smiles, SmilesEmbedding FROM Opm.VectorTable")
            ).fetchall()
        return AllPersistenceResponse(
            result=[SmilesVectorCosine(smiles=row[0], embedding=self.str_to_emb(row[1])) for row in all_rows]
        )
        
