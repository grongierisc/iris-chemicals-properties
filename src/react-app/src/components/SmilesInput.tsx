import React, { useState } from 'react';
import { processSMILES, getSMILESImage } from '../services/api';
import './CompareSmiles.css'; // Reuse the same CSS

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

const SmilesInput: React.FC = () => {
    const [smiles, setSmiles] = useState('');
    const [result, setResult] = useState<Properties | null>(null);
    const [imageUrl, setImageUrl] = useState<string | null>(null);
    const [sessionId, setSessionId] = useState<string | null>(null);
    const [sessionIdImg, setSessionIdImg] = useState<string | null>(null);
    const [isLoading, setIsLoading] = useState(false);

    const handleSubmit = async (e: React.FormEvent) => {
        e.preventDefault();
        setIsLoading(true);
        try {
            const { data: data, sessionId: processSessionId } = await processSMILES(smiles);
            setResult(data.properties);
            setSessionId(processSessionId);
            
            const { data: imageBlob, sessionId: imgSessionId } = await getSMILESImage(smiles);
            setSessionIdImg(imgSessionId)
            setImageUrl(URL.createObjectURL(imageBlob));
        } catch (error) {
            console.error('Error:', error);
        } finally {
            setIsLoading(false);
        }
    };

    const PropertyDisplay = ({ properties }: { properties: Properties }) => (
        <div className="property-section">
            <h3>Properties</h3>
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

    return (
        <div className="smiles-input compare-smiles-container">
            {isLoading && (
                <div className="loading-overlay">
                    <div className="loading-spinner"></div>
                </div>
            )}
            <form onSubmit={handleSubmit} className="input-section">
                <input
                    type="text"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="Enter SMILES string"
                />
                <button type="submit">Process</button>
            </form>

            {result && (
                <div className="results-container">
                    <div className="properties-row">
                        {imageUrl && (
                            <div className="image-container">
                                <img src={imageUrl} alt="Molecule visualization" />
                            </div>
                        )}
                        <PropertyDisplay properties={result} />
                    </div>
                </div>
            )}

            {sessionId && (
                <div className="session-link">
                    <a 
                        href={`http://localhost:53795/csp/irisapp/EnsPortal.VisualTrace.zen?SESSIONID=${sessionId}`}
                        target="_blank"
                        rel="noopener noreferrer"
                    >
                        View Message Trace for Smiles Extractor
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
                        View Message Trace for Img Smiles
                    </a>
                </div>
            )}
        </div>
    );
};

export default SmilesInput;
