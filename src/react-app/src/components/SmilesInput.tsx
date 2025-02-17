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

    const handleSubmit = async (e: React.FormEvent) => {
        e.preventDefault();
        try {
            const data = await processSMILES(smiles);
            setResult(data.properties);
            
            const imageBlob = await getSMILESImage(smiles);
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
        <div className="smiles-input compare-smiles-container">
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
        </div>
    );
};

export default SmilesInput;
