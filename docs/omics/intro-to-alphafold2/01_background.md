## Protein Organization

- **Primary Structure**: amino acid sequence
- **Secondary Structure**: amino acid sequences linked by hydrogen bonds
- **Tertiary Structure**: organization of secondary structures
- **Quaternary Structure**: organization of multiple amino acid chains

![](images/protein_org.jpg)

## The Importance of Protein Structure

- Can help determine what a protein does
- Often more conserved than the amino acid sequences that form them

![](images/different_prot_struct.jpg)

## Protein Structure Problem

As of now there are a few different ways to predict protein structure in the lab:

!!! info ""
    - X-ray Crystallography
    - Nuclear Magnetic Resonance (NMR) Spectroscopy
    - 3D Electron Microscopy
    
However there are 100,000,000 known distinct proteins, each with a unique structure that determines function. Determining protein structure is time consuming and as such only a small fraction of exact 3D structures are known.

## Levinthal's Paradox

- Finding the native folded state of a protein by random searching of all possible configurations would take an enormous amount of time
- However, proteins can often fold within seconds
- Meaning some process must be guiding this folding

![](image/levinthals_paradox.png)

## Using Sequence To Predict Structure

- Instead of laboratory experimentation, there have been massive efforts to use a protein’s sequence to determine structure
- In 1994, the Critical Assessment of Structure Protein (CASP) was established as a biennial assessment of methods to predict structure from sequence

![](images/seq_to_structure.png)


## Enter AlphaFold2 

- Google’s DeepMind team Entered AlphaFold 2 in CASP14 
- Achieved a median Global Distance Test Score of 92.4 
    - This score essentially says: How close is my predicited structure to the known structure?
- AlphaFold 2 works by:
    - starts with a user's query protein sequence
    - finding similar sequences to that query
    - aligns these sequences to create a mutliple sequence alignment
    - uses available structure data based on query sequence to create initial distances between residues
    - uses a neural network to iteratively update the distances between residues by using information from the sequence alignment
    - passes this to another neural network to determine how these residues are oriented in 3D space

![](images/af2_workflow.png)

## Our Data

Today we will be comparing AlphaFold2's structure prediction of Proliferating Cell Nuclear Antigen (PCNA) and DNA ligase 1 (LIG1)!

