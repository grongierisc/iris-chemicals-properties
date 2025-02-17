import React, { useState } from 'react';
import { compareSMILES, compareSMILESImage } from '../services/api';
import './CompareSmiles.css';

interface Properties {
  iupac_name: string;
  formula: string;
  mw: number;
  smiles: string;
  clogp: number;
  clogd: number;
  tpsa: number;
  pka: number;
  pka_type: string;
  h_donor: number;
  h_acceptor: number;
  heavy_atom_count: number;
  rotatable_bonds: number;
}

interface ComparisonResult {
  prop_a: Properties;
  prop_b: Properties;
  diff_prop: Record<keyof Properties, [any, any] | null>;
}

const CompareSmiles: React.FC = () => {
  const [smiles1, setSmiles1] = useState('');
  const [smiles2, setSmiles2] = useState('');
  const [result, setResult] = useState<ComparisonResult | null>(null);
  const [image, setImage] = useState<string | null>(null);

  const handleCompare = async () => {
    try {
      const data = await compareSMILES(smiles1, smiles2);
      setResult(data);

      const imgBlob = await compareSMILESImage(smiles1, smiles2);
      setImage(URL.createObjectURL(imgBlob));
    } catch (error) {
      console.error('Error:', error);
    }
  };

  const PropertyDisplay = ({ properties, title }: { properties: Properties; title: string }) => (
    <div className="property-section">
      <h3>{title}</h3>
      <div className="property-grid">
        {Object.entries(properties).map(([key, value]) => (
          <div key={key} className="property-item">
            <span className="property-label">{key.replace(/_/g, ' ')}:</span>
            <span className="property-value">{value}</span>
          </div>
        ))}
      </div>
    </div>
  );

  const DifferenceDisplay = ({ differences }: { differences: ComparisonResult['diff_prop'] }) => (
    <div className="difference-section">
      <h3>Differences</h3>
      <div className="property-grid">
        {Object.entries(differences).map(([key, value]) => {
          if (value === null) return null;
          return (
            <div key={key} className="property-item">
              <span className="property-label">{key.replace(/_/g, ' ')}:</span>
              <span className="property-value">{value.join(' → ')}</span>
            </div>
          );
        })}
      </div>
    </div>
  );

  return (
    <div className="compare-smiles-container">
      <div className="input-section">
        <input
          type="text"
          value={smiles1}
          onChange={(e) => setSmiles1(e.target.value)}
          placeholder="Enter first SMILES"
        />
        <input
          type="text"
          value={smiles2}
          onChange={(e) => setSmiles2(e.target.value)}
          placeholder="Enter second SMILES"
        />
        <button onClick={handleCompare}>Compare</button>
      </div>

      {result && (
        <div className="results-container">
          <div className="properties-row">
            <PropertyDisplay properties={result.prop_a} title="Properties A" />
            {image && (
              <div className="image-container">
                <img src={image} alt="Comparison visualization" />
              </div>
            )}
            <PropertyDisplay properties={result.prop_b} title="Properties B" />
          </div>
          <DifferenceDisplay differences={result.diff_prop} />
        </div>
      )}
    </div>
  );
};

export default CompareSmiles;
