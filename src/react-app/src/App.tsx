import React, { useState, useEffect } from 'react';
import './App.css';
import Ribbon from './components/Ribbon';
import SmilesInput from './components/SmilesInput';
import SDFUpload from './components/SDFUpload';
import CompareSmiles from './components/CompareSmiles';
import CompareSDF from './components/CompareSDF';
import CompareMixed from './components/CompareMixed';
import SubRibbon from './components/SubRibbon';
import SmilesExtractor from './components/overview/SmilesExtractor';
import SDFExtractor from './components/overview/SDFExtractor';
import CompareOverview from './components/overview/CompareOverview';
import VectorSearch from './components/overview/VectorSearch';
import Architecture from './components/overview/Architecture';
import SimilaritySearch from './components/vector/SimilaritySearch';

function App() {
  const [activeMainTab, setActiveMainTab] = useState('overview');
  const [activeSubTab, setActiveSubTab] = useState('smilesExtractor');

  useEffect(() => {
    // Update subtab when main tab changes
    if (activeMainTab === 'single') {
      setActiveSubTab('smiles');
    } else if (activeMainTab === 'compare') {
      setActiveSubTab('compareSmiles');
    } else if (activeMainTab === 'overview') {
      setActiveSubTab('smilesExtractor');
    } else if (activeMainTab === 'vector') {
      setActiveSubTab('similaritySearch');
    }
  }, [activeMainTab]);

  const renderContent = () => {
    switch (activeSubTab) {
      case 'smiles':
        return <SmilesInput />;
      case 'sdf':
        return <SDFUpload />;
      case 'compareSmiles':
        return <CompareSmiles />;
      case 'compareSDF':
        return <CompareSDF />;
      case 'compareMixed':
        return <CompareMixed />;
      case 'smilesExtractor':
        return <SmilesExtractor />;
      case 'sdfExtractor':
        return <SDFExtractor />;
      case 'compareOverview':
        return <CompareOverview />;
      case 'vectorSearch':
        return <VectorSearch />;
      case 'architecture':
        return <Architecture />;
      case 'similaritySearch':
        return <SimilaritySearch />;
      default:
        return <SmilesExtractor />;
    }
  };

  return (
    <div className="App">
      <header className="App-header">
        <h1>Chemical Properties Viewer</h1>
      </header>
      <Ribbon activeMainTab={activeMainTab} onMainTabChange={setActiveMainTab} />
      <SubRibbon 
        activeMainTab={activeMainTab}
        activeSubTab={activeSubTab}
        onSubTabChange={setActiveSubTab}
      />
      <main>
        {renderContent()}
      </main>
    </div>
  );
}

export default App;
