import React from 'react';

const SmilesExtractor: React.FC = () => {
  return (
    <div className="overview-section">
      <h2>SMILES Extractor Overview</h2>
      <ul>
        <li>Extract properties with RDKit</li>
        <li>Extract the IUAPC name</li>
        <li>Estimate the molecular pKa</li>
      </ul>
    </div>
  );
};

export default SmilesExtractor;
