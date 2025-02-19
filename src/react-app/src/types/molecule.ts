export interface Molecule {
  smiles: string;
  embedding: number[];
  cosine?: number;
}

export interface SearchResult {
  smiles: string;
  embedding: number[];
  cosine: number;
}
