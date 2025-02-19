import React, { useState, useEffect } from 'react';
import { createPersistence, deletePersistence, getAllPersistence, getSimilarMolecules } from '../../services/api';
import { Molecule, SearchResult } from '../../types/molecule';

const SimilaritySearch: React.FC = () => {
  const [molecules, setMolecules] = useState<Molecule[]>([]);
  const [searchSmiles, setSearchSmiles] = useState('');
  const [newSmiles, setNewSmiles] = useState('');
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    loadMolecules();
  }, []);

  const loadMolecules = async () => {
    try {
      setLoading(true);
      const { data } = await getAllPersistence();
      setMolecules(data || []);
    } catch (error) {
      console.error('Error loading molecules:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleAddMolecule = async () => {
    if (!newSmiles.trim()) return;
    try {
      setLoading(true);
      await createPersistence(newSmiles);
      await loadMolecules();
      setNewSmiles('');
    } catch (error) {
      console.error('Error adding molecule:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleSearch = async () => {
    if (!searchSmiles.trim()) {
      await loadMolecules();
      return;
    }
    try {
      setLoading(true);
      const { data } = await getSimilarMolecules(searchSmiles, 0);
      const searchResults = data.result as SearchResult[];
      setMolecules(searchResults.map(result => ({
        smiles: result.smiles,
        embedding: result.embedding || [],
        cosine: result.cosine
      })));
    } catch (error) {
      console.error('Error searching molecules:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleDelete = async (smiles: string) => {
    try {
      setLoading(true);
      await deletePersistence(smiles);
      await loadMolecules();
    } catch (error) {
      console.error('Error deleting molecule:', error);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="p-4">
          <h1 className="text-2xl font-bold mb-6">Molecule Database</h1>

          <div className="grid grid-cols-2 gap-4 mb-6">
              <div className="p-4 border rounded">
              <h2 className="text-lg font-semibold mb-3">Add New Molecule</h2>
              <div className="flex gap-2">
                  <input
                      type="text"
                      value={newSmiles}
                      onChange={(e) => setNewSmiles(e.target.value)}
                      placeholder="Enter SMILES"
                      className="flex-1 p-2 border rounded" />
                  <button
                      onClick={handleAddMolecule}
                      disabled={loading || !newSmiles}
                      className="px-4 py-2 bg-blue-500 text-white rounded disabled:bg-gray-400"
                  >
                      Add
                  </button>
              </div>
          </div>

          <div className="p-4 border rounded">
              <h2 className="text-lg font-semibold mb-3">Search Molecules</h2>
              <div className="flex gap-2">
                  <input
                      type="text"
                      value={searchSmiles}
                      onChange={(e) => setSearchSmiles(e.target.value)}
                      placeholder="Search by SMILES"
                      className="flex-1 p-2 border rounded" />
                  <button
                      onClick={handleSearch}
                      disabled={loading}
                      className="px-4 py-2 bg-green-500 text-white rounded disabled:bg-gray-400"
                  >
                      Search
                  </button>
              </div>
          </div>
      </div><div className="overflow-x-auto">
              <table className="min-w-full border rounded">
                  <thead className="bg-gray-50">
                      <tr>
                          <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">SMILES</th>
                          <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Embedding</th>
                          {searchSmiles && <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Similarity</th>}
                          <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Actions</th>
                      </tr>
                  </thead>
                  <tbody className="bg-white divide-y divide-gray-200">
                      {molecules.map((molecule) => (
                          <tr key={molecule.smiles}>
                              <td className="px-6 py-4 whitespace-nowrap">
                                  <div className="text-sm text-gray-900">{molecule.smiles}</div>
                              </td>
                              <td className="px-6 py-4">
                                  <div className="text-sm text-gray-500 truncate max-w-md">
                                      [{molecule.embedding.slice(0, 3).join(', ')}...]
                                  </div>
                              </td>
                              {searchSmiles && (
                                  <td className="px-6 py-4 whitespace-nowrap">
                                      <div className="text-sm text-gray-900">
                                          {molecule.cosine?.toFixed(4) || '-'}
                                      </div>
                                  </td>
                              )}
                              <td className="px-6 py-4 whitespace-nowrap">
                                  <button
                                      onClick={() => handleDelete(molecule.smiles)}
                                      className="text-red-600 hover:text-red-900"
                                  >
                                      Delete
                                  </button>
                              </td>
                          </tr>
                      ))}
                  </tbody>
              </table>
          </div>
    </div>
  );
};

export default SimilaritySearch;
