export interface Molecule {
  smiles: string;
  embedding: number[];
  cosine?: number;
  cosine_random?: number;
  iupacName?: string;
}

