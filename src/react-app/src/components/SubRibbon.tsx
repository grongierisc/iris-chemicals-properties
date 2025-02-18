import React from 'react';
import '../styles/Ribbon.css';

interface SubRibbonProps {
  activeMainTab: string;
  activeSubTab: string;
  onSubTabChange: (tab: string) => void;
}

const SubRibbon: React.FC<SubRibbonProps> = ({ activeMainTab, activeSubTab, onSubTabChange }) => {
  const subTabs = {
    single: [
      { id: 'smiles', label: 'SMILES' },
      { id: 'sdf', label: 'SDF' }
    ],
    compare: [
      { id: 'compareSmiles', label: 'Compare SMILES' },
      { id: 'compareSDF', label: 'Compare SDF' },
      { id: 'compareMixed', label: 'Compare SMILES vs SDF' }
    ]
  };

  return (
    <div className="ribbon sub-ribbon">
      {subTabs[activeMainTab as keyof typeof subTabs]?.map(tab => (
        <button
          key={tab.id}
          className={activeSubTab === tab.id ? 'active' : ''}
          onClick={() => onSubTabChange(tab.id)}
        >
          {tab.label}
        </button>
      ))}
    </div>
  );
};

export default SubRibbon;
