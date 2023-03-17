## AlphaFold2 Accuracy Assessment

We can assess the accuracy of the AlphaFold prediction using:

- Predicted Local Distance Difference Test (pLDDT)
- Predicted Alignment Error

## Predicted Local Distance Difference Test (pLDDT)

- per-residue confidence metric  ranging from 0-100 (100 being the highest confidence)
- Regions below 50 could indicate disordered regions

![](images/plddt.png)

## Predicted Alignment Error (PAE)

- The Predicted Alignment Error (PAE) gives us an expected distance error based on each residue.
- If we are more confident that the distance between two residues is accurate, then the PAE is lower (darker green). If we are less confident that the distance between two residues is accurate, the PAE is higher (lighter green)

![](images/pae.png)
