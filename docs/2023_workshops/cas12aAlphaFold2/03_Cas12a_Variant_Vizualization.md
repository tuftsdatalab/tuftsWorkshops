## Cas12a Variant Visualization

- In the previous slide we plotted our MSA alignment, the pLDDT scores, and the predicted alignement error. However, it is also useful to visualize the actual predicted protein structure and compare it to the known structure if there is one. Here we use a software called PyMOL to do just that:

![](images/pymolOverview.png)

- Here we see that PyMOL takes either the PDB ID or a PDB file and creates a vizualization for us to examine. If you have not done so already please [download PyMOL](https://access.tufts.edu/pymol){:target="_blank" rel="noopener"} and open the app. You should see a window like the follwing:

![](images/pymolSession.png)

- Here we have a:
  - **History Window** with log of previous commands
  - **Command Interface** to enter PyMOL commands
  - **List of Objects Loaded** which list of objects/proteins that have been loaded into PyMOL
  - **Visualization Window** to visualize protiens loaded into PyMOL

- Let's try on our data!

## Download AlphaFold Output

- First we will need to download our predicted structure pdb file for the mutant, Cas12a mut2-CWF. To this go to Files > Home Directory:

![](images/homeDir.png)

- Now, click on cas12a_af2 and download the following file, `mut2cwf.pdb`:

IMAGES TO BE INSERTED HERE

![](images/)

## Visualization with PyMOL

- To visualize this protein structure in PyMOL, open PyMOL on your computer
- Go to File > Open - then choose your pdb file.
- We will now align this structure with the wild type of Cas12a(PDB: 5XUS). So in the PyMOL command prompt enter:

```bash
fetch 5xus
```

```bash
align 5xus, mut2cwf
```

- Now that we have aligned our structures, let's select residues on the Bridge Helix region on the mutan, Cas12a mut2-CWF, and Cas12a wild-type:

```bash
select resid 890+885+884
```

- With these residues selected we can color them to visualize them easier:

```bash
color sele, red
```

- Let's now zoom into this region:

```bash
zoom sele
```

??? question "Can you spot any disordered regions that AlphaFold2 may not have predicted well? Include an image displaying one of these regions. "

??? question "[Ma et al. 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9825149/) note that the Cas12a mut2-CW has a more open active site than the wild type Cas12a. Zoom into the catalytic site for Cas12a mut2c-W and include an image displaying this catalytic site (Hint: the catalytic site will contain the following residues: 863, 952, 965, and 1214). Do you agree that the catalytic site appears more open in comparison to the wild type structure (5XUS)?"

