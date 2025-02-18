import React, { useState } from 'react';
import { compareMixed, compareMixedImage } from '../services/api';
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
  prop_smiles: Properties;
  prop_sdf: Properties;
  diff_prop: Record<keyof Properties, [any, any] | null>;
}

const CompareMixed: React.FC = () => {
  const [smiles, setSmiles] = useState('');
  const [file, setFile] = useState<File | null>(null);
  const [result, setResult] = useState<ComparisonResult | null>(null);
  const [image, setImage] = useState<string | null>(null);
  const [sessionId, setSessionId] = useState<string | null>(null);
  const [sessionIdImg, setSessionIdImg] = useState<string | null>(null);

  const handleCompare = async () => {
    if (!smiles || !file) return;

    try {
      const { data, sessionId: compareSessionId } = await compareMixed(smiles, file);
      setResult(data);
      setSessionId(compareSessionId);

      const { data: imgBlob, sessionId: imgSessionId } = await compareMixedImage(smiles, file);
      setSessionIdImg(imgSessionId);
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
          <label htmlFor="smiles">Enter SMILES string</label>
          <input
            type="text"
            id="smiles"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="Enter SMILES"
            className="smiles-input"
          />
        </div>
        <div className="file-input-wrapper">
          <label htmlFor="file">Select SDF file</label>
          <input
            id="file"
            type="file"
            onChange={(e) => setFile(e.target.files?.[0] || null)}
            accept=".sdf"
          />
          {file && <span className="file-name">{file.name}</span>}
        </div>
        <button onClick={handleCompare} disabled={!smiles || !file}>
          Compare
        </button>
      </div>

      {result && (
        <div className="results-container">
          <div className="properties-row">
            <PropertyDisplay properties={result.prop_smiles} title="SMILES Properties" />
            {image && (
              <div className="image-container">
                <img src={image} alt="Comparison visualization" />
              </div>
            )}
            <PropertyDisplay properties={result.prop_sdf} title="SDF Properties" />
          </div>
          <DifferenceDisplay differences={result.diff_prop} />
        </div>
      )}

      {sessionId && (
        <div className="session-link">
          <a 
            href={`http://localhost:53795/csp/irisapp/EnsPortal.VisualTrace.zen?SESSIONID=${sessionId}`}
            target="_blank"
            rel="noopener noreferrer"
          >
            View Message Trace for Mixed Comparison
          </a>
        </div>
      )}

      {sessionIdImg && (
        <div className="session-link">
          <a 
            href={`http://localhost:53795/csp/irisapp/EnsPortal.VisualTrace.zen?SESSIONID=${sessionIdImg}`}
            target="_blank"
            rel="noopener noreferrer"
          >
            View Message Trace for Mixed Comparison Image
          </a>
        </div>
      )}
    </div>
  );
};

export default CompareMixed;
