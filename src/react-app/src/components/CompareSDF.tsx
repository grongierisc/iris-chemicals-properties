import React, { useState } from 'react';
import { compareSDF, compareSDFImage } from '../services/api';
import './CompareSDF.css';

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

const CompareSDF: React.FC = () => {
  const [file1, setFile1] = useState<File | null>(null);
  const [file2, setFile2] = useState<File | null>(null);
  const [result, setResult] = useState<ComparisonResult | null>(null);
  const [image, setImage] = useState<string | null>(null);

  const handleCompare = async () => {
    if (!file1 || !file2) return;

    try {
      const data = await compareSDF(file1, file2);
      setResult(data);

      const imgBlob = await compareSDFImage(file1, file2);
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
    <div className="compare-sdf-container">
      <div className="input-section">
        <div className="file-input-wrapper">
          <label htmlFor="file1">Select first SDF file</label>
          <input
            id="file1"
            type="file"
            onChange={(e) => setFile1(e.target.files?.[0] || null)}
            accept=".sdf"
          />
          {file1 && <span className="file-name">{file1.name}</span>}
        </div>
        <div className="file-input-wrapper">
          <label htmlFor="file2">Select second SDF file</label>
          <input
            id="file2"
            type="file"
            onChange={(e) => setFile2(e.target.files?.[0] || null)}
            accept=".sdf"
          />
          {file2 && <span className="file-name">{file2.name}</span>}
        </div>
        <button onClick={handleCompare} disabled={!file1 || !file2}>
          Compare
        </button>
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

export default CompareSDF;
