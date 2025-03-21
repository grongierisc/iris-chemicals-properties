export interface Molecule {
  smiles: string;
  embedding: number[];
  embedding_random: number[];
  cosine?: number;
  cosine_random?: number;
  iupacName?: string;
  img?: Blob;
}

