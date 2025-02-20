import React from 'react';
import './Overview.css';

const VectorSearch: React.FC = () => {
  return (
    <div>
      <h2>Vector Search Overview</h2>
      <div className="overview-section">
      <p>
        Store the SMILES representation of the chemicals in a vector database<br />
        Using a ChemBert model to extract the vector representation of the chemicals<br />
        Similarities between the vectors are calculated using the cosine similarity<br />
        <br />
        Example of a query: <br />
        <p className="sql-highlight">
          SELECT TOP 5 Smiles, SmilesEmbedding, VECTOR_COSINE(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')) FROM Opm.VectorTable
          ORDER BY VECTOR_COSINE(SmilesEmbedding,EMBEDDING(:smiles,'opm-transformers-config')) DESC
        </p>
        The EMBEDDING uses two parameters: <br />
        <ul>
        <li>the SMILES string and the model name</li> 
        <li>the model name is configured like this</li>
        </ul>
        <p className="sql-highlight">
          INSERT INTO %Embedding.Config (Name, Configuration, EmbeddingClass, VectorLength, Description)<br />
          VALUES (&apos;opm-transformers-config&apos;,<br />
            &apos;&#123;&quot;modelName&quot;:&quot;seyonec/PubChem10M_SMILES_BPE_450k&quot;,<br />
            &quot;hfCachePath&quot;:&quot;/usr/irissys/hfCache&quot;&#125;&apos;,<br />
            &apos;Opm.Embedding&apos;,<br />
            768,<br />
            &apos;embedding model&apos;<br />
          )
        </p>
      </p>
    </div>
    </div>
  );
};

export default VectorSearch;
