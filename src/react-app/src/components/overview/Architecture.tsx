import React from 'react';
import chemSchema from '../../imgs/chem-schema.png';

const Architecture: React.FC = () => {
  return (
    <div className="container mt-4">
      <h2>System Architecture</h2>
      <div className="text-center mt-4">
        <img 
          src={chemSchema} 
          alt="Chemical Properties System Architecture" 
          style={{ maxWidth: '100%', height: 'auto' }}
        />
      </div>
    </div>
  );
};

export default Architecture;
