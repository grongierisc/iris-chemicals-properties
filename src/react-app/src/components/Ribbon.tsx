import React from 'react';
import '../styles/Ribbon.css';

interface RibbonProps {
  activeTab: string;
  onTabChange: (tab: string) => void;
}

const Ribbon: React.FC<RibbonProps> = ({ activeTab, onTabChange }) => {
  const tabs = [
    { id: 'smiles', label: 'SMILES' },
    { id: 'sdf', label: 'SDF' },
    { id: 'compareSmiles', label: 'Compare SMILES' },
    { id: 'compareSDF', label: 'Compare SDF' },
    { id: 'compareMixed', label: 'Compare SMILES vs SDF' }
  ];

  return (
    <div className="ribbon">
      {tabs.map(tab => (
        <button
          key={tab.id}
          className={activeTab === tab.id ? 'active' : ''}
          onClick={() => onTabChange(tab.id)}
        >
          {tab.label}
        </button>
      ))}
    </div>
  );
};

export default Ribbon;
