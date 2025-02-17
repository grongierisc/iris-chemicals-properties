import React, { useState } from 'react';
import './App.css';
import Ribbon from './components/Ribbon';
import SmilesInput from './components/SmilesInput';
import SDFUpload from './components/SDFUpload';
import CompareSmiles from './components/CompareSmiles';
import CompareSDF from './components/CompareSDF';
import CompareMixed from './components/CompareMixed';

function App() {
  const [activeTab, setActiveTab] = useState('smiles');

  const renderContent = () => {
    switch (activeTab) {
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
      <Ribbon activeTab={activeTab} onTabChange={setActiveTab} />
      <main>
        {renderContent()}
      </main>
    </div>
  );
}

export default App;
