import React from 'react';
import '../styles/Ribbon.css';

interface RibbonProps {
  activeMainTab: string;
  onMainTabChange: (tab: string) => void;
}

const Ribbon: React.FC<RibbonProps> = ({ activeMainTab, onMainTabChange }) => {
  const mainTabs = [
    { id: 'overview', label: 'Overview' },
    { id: 'single', label: 'Single Molecule' },
    { id: 'compare', label: 'Compare Molecules' }
  ];

  return (
    <div className="ribbon main-ribbon">
      {mainTabs.map(tab => (
        <button
          key={tab.id}
          className={activeMainTab === tab.id ? 'active' : ''}
          onClick={() => onMainTabChange(tab.id)}
        >
          {tab.label}
        </button>
      ))}
    </div>
  );
};

export default Ribbon;
