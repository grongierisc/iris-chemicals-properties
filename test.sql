INSERT INTO %Embedding.Config (Name, Configuration, EmbeddingClass, VectorLength, Description)
  VALUES ('opm-transformers-config',
          '{"modelName":"seyonec/PubChem10M_SMILES_BPE_450k",
          "hfCachePath":"/usr/irissys/hfCache"}',
          'Opm.Embedding',
          768,
          'embedding model')
          ---
INSERT INTO %Embedding.Config (Name, Configuration, EmbeddingClass, Description)
  VALUES ('sentence-transformers-config',
          '{"modelName":"sentence-transformers/all-MiniLM-L6-v2",
            "hfCachePath":"/usr/irissys/hfCache",
            "maxTokens": 256,
            "checkTokenCount": true}',
          '%Embedding.SentenceTransformers',
          'a small SentenceTransformers embedding model')
          ---

CREATE TABLE Opm.VectorTable (
  Smiles VARCHAR(250),
  SmilesEmbedding EMBEDDING('opm-transformers-config','Smiles')
)
---
SELECT TOP 5 Smiles, SmilesEmbedding, VECTOR_COSINE(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')),VECTOR_COSINE(SmilesEmbedding,EMBEDDING(:smiles,'sentence-transformers-config')) FROM Opm.VectorTable
ORDER BY VECTOR_COSINE(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')) DESC
          ---
INSERT INTO Opm.VectorTable (Smiles)
  VALUES ('CCO')
          ---
SELECT * FROM Opm.VectorTable