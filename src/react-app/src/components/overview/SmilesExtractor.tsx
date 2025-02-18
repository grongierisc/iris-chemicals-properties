import React from 'react';

const SmilesExtractor: React.FC = () => {
  return (
    <div className="overview-section">
      <h2>SMILES Extractor Overview</h2>
      <p>The SMILES Extractor allows you to:</p>
      <ul>
        <li>Input SMILES strings for chemical compounds</li>
        <li>Extract and visualize molecular properties</li>
        <li>Generate 2D structure representations</li>
        <li>Calculate chemical descriptors</li>
      </ul>
    </div>
  );
};

export default SmilesExtractor;
