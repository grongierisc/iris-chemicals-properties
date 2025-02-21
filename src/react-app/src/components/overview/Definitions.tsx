import React from "react";

const Definitions: React.FC = () => {
  return (
    <div className="overview-section">
      <p>
      <h2>Definitions</h2>
      <h3>SMILES</h3>
        Simplified Molecular Input Line Entry System (SMILES) is a line notation for entering and representing molecules and reactions. <br />
        Examples: <br />
            <ul>
            <li>Water: <code>O</code></li>
            <li>Carbon dioxide: <code>O=C=O</code></li>
            <li>Caffeine: <code>CC1=CN(C(=O)N=C1N)C</code></li>
            <li>Aspirin: <code>CC(=O)OC1=CC=CC=C1C(=O)O</code></li>
            <li>Paracetamol: <code>CC(=O)NC1=CC=C(C=C1)O</code></li>
            <li>Acetaminophen: <code>CC(=O)Nc1ccc(C(=O)O)cc1</code></li>
            <li>Penicillin: <code>CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)[C@H](C3=CC=CC=C3)N)C(=O)O)C</code></li>
            </ul>
      <h3>SDF</h3>
        Structure Data File (SDF) is a file format used to store chemical structures and associated data. It is a standard file format used in the chemical industry
      <h3>RDKit</h3>
        RDKit is an open-source cheminformatics software toolkit that provides a wide range of functionality for working with chemical structures. It is written in C++ and Python and is widely used in the pharmaceutical industry and in academic research.
      <h3>ChemBert</h3>
        ChemBERT is a transformer-based model that is pre-trained on a large corpus of chemical data. It is trained to predict chemical properties and generate chemical representations. ChemBERT is a powerful tool for cheminformatics and can be used for a wide range of tasks, including chemical property prediction, molecular generation, and virtual screening.
      <h3>pKa</h3>
        The pKa of a compound is a measure of its acidity or basicity.<br />
        It's an empirical measure of the strength of an acid in solution.
    </p>    
    </div>
  );
}

export default Definitions;