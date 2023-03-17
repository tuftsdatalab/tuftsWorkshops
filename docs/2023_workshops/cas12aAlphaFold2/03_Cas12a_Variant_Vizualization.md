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

- First we will need to download our predicted structure pdb files. To this go to Files > Home Directory:

![](images/homeDir.png)

- Now, click on cas12a_af2 and download each of these files that end in `.pdb`:

IMAGES TO BE INSERTED HERE

![](images/)

## Visualization with PyMOL

- To visualize these protein structures in PyMOL, open PyMOL on your computer
- Go to File > Open - then choose your pdb files.
- We will now align each of these structures so we can visualize the region of interest. So in the PyMOL command prompt enter:

```bash
align mut2bw, mut2cw
```

```bash
align mut2cwf, mut2cw
```


- Now that we have aligned all of our structures, let's disable 
- Let's visualize the Bridge Helix region on these Cas12a mutants. 

```bash
select resid 890+885+884
```

```bash
color sele, red
```

```bash
zoom sele
```
