# iris-chemicales-properties

Schema :

![alt text](https://github.com/grongierisc/iris-chemicals-properties/blob/master/misc/chem-schema.png?raw=true)

## Description

This project aim to showcase the use of the [IRIS](https://www.intersystems.com/products/intersystems-iris/) to extract chemicals properties and use them in a vector database.

The project is composed of :

- Smiles processing
  - RDKit is used to extract the properties of the chemicals
  - Extracted the IUAPC name of the chemicals
  - Estimate the molecular pKa of the chemicals
- SDF processing
  - Simple extraction of the properties of the chemicals present in the SDF file
- Compare the properties of the chemicals
  - Based on the SMILES processing or the SDF processing
- Store the SMILES representation of the chemicals in a vector database
  - Using a ChemBert model to extract the vector representation of the chemicals

## SMILES Extraction

The SMILES extraction is done in 3 steps :

- Extract properties with RDKit
- Extract the IUAPC name
- Estimate the molecular pKa

### RDKit

The SMILES extraction is done using the RDKit library. The properties extracted are :

- Molecular weight
- Number of atoms
- Number of bonds
- Number of H-bonds
- And other standard properties

### IUAPC name

The IUAPC name is extracted using the service provided by cactus.nci.nih.gov.

### Molecular pKa

The molecular pKa is estimated using two models :

- one for the classification of the pKa
  - acidic
  - basic
  - neutral
- one predicting the pKa value

## SDF Extraction

The SDF extraction is a simple extraction of the properties of the chemicals present in the SDF file.

## Compare the properties of the chemicals

The properties of the chemicals can be compared based on the SMILES processing or the SDF processing.

## Store the SMILES representation of the chemicals in a vector database

The SMILES representation of the chemicals is stored in a vector database. The vector representation is extracted using a ChemBert model.