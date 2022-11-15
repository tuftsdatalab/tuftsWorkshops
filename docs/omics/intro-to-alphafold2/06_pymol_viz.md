## Visualizing With Pymol

- In the previous slide we plotted our MSA alignment, the pLDDT scores, and the predicted alignement error. However, it is also useful to visualize the actual predicted protein structure and compare it to the known structure if there is one. Here we use a software called PyMOL to do just that:

![](images/pymolOverview.png)

- Here we see that PyMOL takes either the PDB ID or a PDB file and creates a vizualization for us to examine. If you have not done so already please [download PyMOL](https://pymol.org/2/) and open the app. You should see a window like the follwing:

![](images/pymolSession.png)

- Here we have a:
  - **History Window** with log of previous commands
  - **Command Interface** to enter PyMOL commands
  - **List of Objects Loaded** which list of objects/proteins that have been loaded into PyMOL
  - **Visualization Window** to visualize protiens loaded into PyMOL

- Let's try on our data!

## Download AlphaFold Output

- First we will need to download our predicted structure pdb files. To this go to Files > Home Directory:

![](images/homeDir.png)

- Then navigate to your AlphaFold Workshop directory, where you will note the two folders we examined earlier that have our AlphaFold Outputs:

![](images/pdbFolders.png)

- In each folder you will note a `visuals` folder and an ID number. Navigate to the ID number folder and download the file "ranked_0.pdb" - this structure is AlphaFold's best prediction of the protein's structure.

![](images/id.png)
![](images/ranked0.png)

- Make sure to add a prefix to each `ranked_0.pdb` file as to not confuse the two. You can rename the `ranked_0.pdb` in the `lig1` folder to `lig1_ranked_0.pdb` and the `ranked_0.pdb` in the `pcna` folder to `pcna_ranked_0.pdb`. 

## AlphaFold Output In PyMOL

- To visualize these protein structures in PyMOL go to File > Open - then choose your pdb files. 
- Loading two objects can make it difficult to see examine both individually, so enter the following command to hide lig1:

```
disable lig1_ranked_0
```
- Now you should only see `pnca_ranked_0`. It would be interesting to see how well this prediction lines up with the know structure for PCNA. We can load that structure with the following command

```
fetch 1axc
```
- To see how well our predicted structure lines up we can compare the two by aligning them:

```
align pcna_ranked_0, 1axc
```

![](images/align.png)

- In the history window you'll note that when we aligned our structures we were given an RMSD or root mean square deviation value. The smaller this value is, the better our two structures have aligned. 
- Now that we have aligned the predicted PCNA structure, repeat these steps to align LIG1 (the PDB ID for LIG1 id `1x9n`).

## AlphaFold2 Limitations

AlphaFold2 attempts to predict protein structures based on available structure data in the [PDB](https://www.rcsb.org/). But do these structures in the [PDB](https://www.rcsb.org/) reflect actual protein structures?


- PDB structures are usually created from experiments where the context of that structure is specific to the study question. 

!!! info ""

    For example there are lots of studies examining what a particular protein structure looks like when bound to ions, when it’s chemically modified, or when its in larger complexes
    
- Protein interactions/multimers might not be captured in the [PDB](https://www.rcsb.org/) database. Given this, AlphaFold2’s multimeric prediction might not be reflective of the true interaction structure.
- Proteins can also contain disordered regions (i.e. loops), which are difficult to crystallize and as such AlphaFold’s prediction of these disordered regions is bound to be poor.
        
AlphaFold2 is indeed a powerful tool but just be aware of what it is prediciting and if any of the items mentioned above interfere with the study question you are using AlphaFold2 to answer!

## References

!!! abstract "References"

    1. [PROTEIN](https://www.genome.gov/genetics-glossary/Protein)
    2. [Protein Function](https://www.nature.com/scitable/topicpage/protein-function-14123348/)
    3. [Analyzing Protein Structure and Function](https://www.ncbi.nlm.nih.gov/books/NBK26820/)
    4. [serial scanning 3D electron microscopy](https://directorsblog.nih.gov/tag/serial-scanning-3d-electron-microscopy/)
    5. [AI revolutions in biology](https://www.embopress.org/doi/full/10.15252/embr.202154046)
    6. [Methods for Determining Atomic Structures](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/methods-for-determining-structure)
    7. [X-ray crystallography](https://simple.wikipedia.org/wiki/X-ray_crystallography)
    8. [AlphaFold](https://deepmind.com/research/case-studies/alphafold)
    9. [Levinthal's paradox](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC48166/)
    10. [AlphaFold: a solution to a 50-year-old grand challenge in biology](https://deepmind.com/blog/article/alphafold-a-solution-to-a-50-year-old-grand-challenge-in-biology)
    11. [Protein Structure Prediction Center](https://predictioncenter.org/)
    12. [Neural network](https://en.wikipedia.org/wiki/Neural_network)
    13. [Understanding LSTMs](https://colah.github.io/posts/2015-08-Understanding-LSTMs/)
    14. [Attention Is All You Need](https://proceedings.neurips.cc/paper/2017/file/3f5ee243547dee91fbd053c1c4a845aa-Paper.pdf)
    15. [Transformer Neural Network: Step-By-Step Breakdown of the Beast](https://towardsdatascience.com/transformer-neural-network-step-by-step-breakdown-of-the-beast-b3e096dc857f)
    16. [FASTA Format](https://en.wikipedia.org/wiki/FASTA_format)
    17. [Multiple Sequence Alignment](https://en.wikipedia.org/wiki/Multiple_sequence_alignment)
    18. [Origins of coevolution between residues distant in protein 3D structures](https://www.pnas.org/content/114/34/9122)
    19. [AlphaFold 2 is here: what’s behind the structure prediction miracle](https://www.blopig.com/blog/2021/07/alphafold-2-is-here-whats-behind-the-structure-prediction-miracle/)
    20. [Alphafold Github](https://github.com/deepmind/alphafold)
    21. [Q9FX77](https://alphafold.com/entry/Q9FX77)
    22. [Limitations of AlphaFold](https://blog.cyclicarx.com/limitations-of-alphafold)


    

