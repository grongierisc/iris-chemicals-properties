import React from 'react';
import chemSchema from '../../imgs/chem-schema.png';

const Architecture: React.FC = () => {
  return (
    <div className="container mt-4">
      <h2>System Architecture</h2>
      <div className='overview-section'>
      <p>
      This project aim to showcase the use of the <a href="http://localhost:53795/csp/irisapp/EnsPortal.ProductionConfig.zen?$NAMESPACE=IRISAPP&">IRIS</a> to 
      extract chemicals properties and use them in a vector database.<br/>
      <br/>
      The project is composed of:<br/>
      <ul>
      <li>
        <strong>Smiles processing</strong>
        <ul>
        <li>RDKit is used to extract the properties of the chemicals</li>
        <li>Extracted the IUAPC name of the chemicals</li>
        <li>Estimate the molecular pKa of the chemicals</li>
        </ul>
      </li>
      <li>
        <strong>SDF processing</strong>
        <ul>
        <li>Simple extraction of the properties of the chemicals present in the SDF file</li>
        </ul>
      </li>
      <li>
        <strong>Compare the properties of the chemicals</strong>
        <ul>
        <li>Based on the SMILES processing or the SDF processing</li>
        </ul>
      </li>
      <li>
        <strong>Store the SMILES representation of the chemicals in a vector database</strong>
        <ul>
        <li>Using a ChemBert model to extract the vector representation of the chemicals</li>
        </ul>
      </li>
      </ul>
      <div className="text-center mt-4">
      <img 
        src={chemSchema} 
        alt="Chemical Properties System Architecture" 
        style={{ maxWidth: '100%', height: 'auto' }}
      />
      </div>
      </p>
    </div>
    </div>
  );
};

export default Architecture;
