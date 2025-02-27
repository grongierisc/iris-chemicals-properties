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
            <li>Caffeine: <code>CC1=CN(C(=O)N=C1N)C</code></li>
            <li>Aspirin: <code>CC(=O)OC1=CC=CC=C1C(=O)O</code></li>
            <li>Paracetamol: <code>CC(=O)NC1=CC=C(C=C1)O</code></li>
            <li>Ibuprofen: <code>CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O</code></li>
            <li>Nurofen: <code>CC(Cc1ccc(cc1)C(C(=O)O)C)C</code></li>
            <li>Codeine: <code>CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@H]3[C@H](C=C4)O</code></li>
            </ul>
      <h3>SDF</h3>
        Structure Data File (SDF) is a file format used to store chemical structures and associated data. It is a standard file format used in the chemical industry
      <h3>RDKit</h3>
        RDKit is an open-source cheminformatics software toolkit that provides a wide range of functionality for working with chemical structures. It is written in C++ and Python and is widely used in the pharmaceutical industry and in academic research.
      <h3>pKa</h3>
        The pKa of a compound is a measure of its acidity or basicity.<br />
        It's an empirical measure of the strength of an acid in solution.
    </p>    
    </div>
  );
}

export default Definitions;