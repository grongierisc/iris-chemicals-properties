import React from 'react';
import './Overview.css';
import cosine from '../../imgs/cosine.png';
import jaccard from '../../imgs/jaccard.png';

const VectorSearch: React.FC = () => {
  return (
    <div>
      <h2>Vector Search Overview</h2>
      <div className="overview-section">
      <p>
        Store the SMILES representation of the chemicals in a vector database<br />
        Using a MiniLM model to extract the vector representation of the chemicals<br />
        and MACCS fingerprints for similarity search<br />
        <br />
        Similarities between the vectors are calculated using the cosine similarity<br />
        and the Jaccard similarity for the MACCS fingerprints<br />
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
        <span className="sql-keyword">INSERT INTO</span> %Embedding.Config (Name, Configuration, EmbeddingClass, Description)
          <span className="sql-keyword">VALUES</span> (&apos;sentence-transformers-config&apos;,
          &apos;&#123;&quot;modelName&quot;:&quot;sentence-transformers/all-MiniLM-L6-v2&quot;,
            &quot;hfCachePath&quot;:&quot;/usr/irissys/hfCache&quot;,
            &quot;maxTokens&quot;: 256,
            &quot;checkTokenCount&quot;: true&#125;&apos;,
          &apos;%Embedding.SentenceTransformers&apos;,
          &apos;a small SentenceTransformers embedding model&apos;)
        </p>
        The aim of this demo is to show how to use the vector search capabilities of IRIS<br />
        Then to explain which embedding model is best suited for the task<br />
        Same for the similarity search algorithm<br />
        <br />
        COSINE distance is defined as the cosine of the angle between two vectors<br />
        <div className="text-center mt-4">
        <img 
          src={cosine} 
          alt="Chemical Properties System Architecture" 
          style={{ maxWidth: '100%', height: 'auto' }}
        />
        </div>
        JACCARD distance is defined as the intersection over the union of the two sets<br />
        <div className="text-center mt-4">
        <img 
          src={jaccard} 
          alt="Chemical Properties System Architecture" 
          style={{ maxWidth: '100%', height: 'auto' }}
        />
        </div>
        <br />
        The MACCS fingerprints are a set of 166 bits that represent the presence or absence of substructures in the molecule<br />
        The MiniLM model is a transformer model that is trained on a large corpus of text, basically this model is not suited for the task<br />
        <br />
      </p>
    </div>
    </div>
  );
};

export default VectorSearch;
