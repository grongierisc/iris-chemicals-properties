import React, { useState } from 'react';
import { uploadSDF, getSDFImage } from '../services/api';
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

const SDFUpload: React.FC = () => {
    const [file, setFile] = useState<File | null>(null);
    const [result, setResult] = useState<Properties | null>(null);
    const [imageUrl, setImageUrl] = useState<string | null>(null);
    const [sessionId, setSessionId] = useState<string | null>(null);
    const [sessionIdImg, setSessionIdImg] = useState<string | null>(null);

    const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (e.target.files && e.target.files[0]) {
            setFile(e.target.files[0]);
        }
    };

    const handleSubmit = async (e: React.FormEvent) => {
        e.preventDefault();
        if (!file) return;

        try {
            const { data, sessionId: processSessionId } = await uploadSDF(file);
            setResult(data.properties);
            setSessionId(processSessionId);
            
            const { data: imageBlob, sessionId: imgSessionId } = await getSDFImage(file);
            setSessionIdImg(imgSessionId);
            setImageUrl(URL.createObjectURL(imageBlob));
        } catch (error) {
            console.error('Error:', error);
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
        <div className="sdf-upload compare-smiles-container">
            <form onSubmit={handleSubmit} className="input-section">
                <input
                    type="file"
                    accept=".sdf"
                    onChange={handleFileChange}
                />
                <button type="submit" disabled={!file}>
                    Upload and Process
                </button>
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
                        View Message Trace for SDF Extractor
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
                        View Message Trace for Img SDF
                    </a>
                </div>
            )}
        </div>
    );
};

export default SDFUpload;
