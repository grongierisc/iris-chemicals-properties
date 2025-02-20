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
          <span className="sql-keyword">SELECT</span> <span className="sql-keyword">TOP</span> 5 Smiles, SmilesEmbedding, 
          <span className="sql-function"> VECTOR_COSINE</span>(SmilesEmbedding,<span className="sql-function">EMBEDDING</span>(:smiles,&apos;opm-transformers-config&apos;)) 
          <span className="sql-keyword"> FROM</span> Opm.VectorTable
          <span className="sql-keyword"> ORDER BY </span> <span className="sql-function">VECTOR_COSINE</span>(SmilesEmbedding,<span className="sql-function">EMBEDDING</span>(:smiles,&apos;opm-transformers-config&apos;)) <span className="sql-keyword">DESC</span>
        </p>
        The EMBEDDING uses two parameters: <br />
        <ul>
        <li>the SMILES string and the model name</li> 
        <li>the model name is configured like this</li>
        </ul>
        <p className="sql-highlight">
          <span className="sql-keyword">INSERT INTO</span> %Embedding.Config (Name, Configuration, EmbeddingClass, VectorLength, Description)<br />
          <span className="sql-keyword">VALUES</span> (&apos;opm-transformers-config&apos;,<br />
            <span className="json-content">&apos;&#123;&quot;modelName&quot;:&quot;seyonec/PubChem10M_SMILES_BPE_450k&quot;,<br />
            &quot;hfCachePath&quot;:&quot;/usr/irissys/hfCache&quot;&#125;&apos;</span>,<br />
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
