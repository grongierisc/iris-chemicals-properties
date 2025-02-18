import React, { useState, useEffect } from 'react';
import './App.css';
import Ribbon from './components/Ribbon';
import SmilesInput from './components/SmilesInput';
import SDFUpload from './components/SDFUpload';
import CompareSmiles from './components/CompareSmiles';
import CompareSDF from './components/CompareSDF';
import CompareMixed from './components/CompareMixed';
import SubRibbon from './components/SubRibbon';

function App() {
  const [activeMainTab, setActiveMainTab] = useState('single');
  const [activeSubTab, setActiveSubTab] = useState('smiles');

  useEffect(() => {
    // Update subtab when main tab changes
    if (activeMainTab === 'single') {
      setActiveSubTab('smiles');
    } else {
      setActiveSubTab('compareSmiles');
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
      default:
        return <SmilesInput />;
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
