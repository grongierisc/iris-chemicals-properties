import React from 'react';
import '../styles/Ribbon.css';

interface SubRibbonProps {
  activeMainTab: string;
  activeSubTab: string;
  onSubTabChange: (tab: string) => void;
}

const SubRibbon: React.FC<SubRibbonProps> = ({ activeMainTab, activeSubTab, onSubTabChange }) => {
  const renderTabs = () => {
    switch (activeMainTab) {
      case 'overview':
        return (
          <>
            <button
              className={activeSubTab === 'architecture' ? 'active' : ''}
              onClick={() => onSubTabChange('architecture')}
            >
              Architecture
            </button>
            <button 
              className={activeSubTab === 'smilesExtractor' ? 'active' : ''}
              onClick={() => onSubTabChange('smilesExtractor')}
            >
              SMILES Extractor
            </button>
            <button 
              className={activeSubTab === 'sdfExtractor' ? 'active' : ''}
              onClick={() => onSubTabChange('sdfExtractor')}
            >
              SDF Extractor
            </button>
            <button 
              className={activeSubTab === 'compareOverview' ? 'active' : ''}
              onClick={() => onSubTabChange('compareOverview')}
            >
              Compare
            </button>
            <button 
              className={activeSubTab === 'vectorSearch' ? 'active' : ''}
              onClick={() => onSubTabChange('vectorSearch')}
            >
              Vector Search
            </button>

          </>
        );
      case 'single':
        return (
          <>
            <button
              className={activeSubTab === 'smiles' ? 'active' : ''}
              onClick={() => onSubTabChange('smiles')}
            >
              SMILES
            </button>
            <button
              className={activeSubTab === 'sdf' ? 'active' : ''}
              onClick={() => onSubTabChange('sdf')}
            >
              SDF
            </button>
          </>
        );
      case 'compare':
        return (
          <>
            <button
              className={activeSubTab === 'compareSmiles' ? 'active' : ''}
              onClick={() => onSubTabChange('compareSmiles')}
            >
              Compare SMILES
            </button>
            <button
              className={activeSubTab === 'compareSDF' ? 'active' : ''}
              onClick={() => onSubTabChange('compareSDF')}
            >
              Compare SDF
            </button>
            <button
              className={activeSubTab === 'compareMixed' ? 'active' : ''}
              onClick={() => onSubTabChange('compareMixed')}
            >
              Compare SMILES vs SDF
            </button>
          </>
        );
      default:
        return null;
    }
  };

  return (
    <div className="ribbon sub-ribbon">
      {renderTabs()}
    </div>
  );
};

export default SubRibbon;
