import React from 'react';

const SmilesExtractor: React.FC = () => {
  return (
    <div className="overview-section">
      <h2>SMILES Extractor Overview</h2>
      <p>The SMILES extraction is done in 3 steps:</p>
      
      <ul>
      <li>Extract properties with RDKit</li>
      <li>Extract the IUAPC name</li>
      <li>Estimate the molecular pKa</li>
      </ul>

      <h3>RDKit</h3>
      <p>The SMILES extraction is done using the RDKit library. The properties extracted are:</p>
      
      <ul>
      <li>Molecular weight</li>
      <li>Number of atoms</li>
      <li>Number of bonds</li>
      <li>Number of H-bonds</li>
      <li>And other standard properties</li>
      </ul>

      <h3>IUAPC name</h3>
      <p>The IUAPC name is extracted using the service provided by cactus.nci.nih.gov.</p>

      <h3>Molecular pKa</h3>
      <p>The molecular pKa is estimated using two models:</p>
      
      <ul>
      <li>one for the classification of the pKa
        <ul>
        <li>acidic</li>
        <li>basic</li>
        <li>neutral</li>
        </ul>
      </li>
      <li>one predicting the pKa value</li>
      </ul>
    </div>
  );
};

export default SmilesExtractor;
