### Introduction to the Microbiome

- The microbiome are the set of genes belonging to the microbiota in a specimen. The term microbiome can also refer to the microbes themselves
- Disturbances in the microbiome are linked to obesity, inflammatory bowel disease, alcoholic and nonalcoholic fatty liver disease, and hepatocellular carcinoma

![](images/microbiome.png)

### Microbiome Variability

- Assessing a microbiome disturbance is not a trivial task as it is highly variable from person to person.
- Large sample sizes, hundreds of patients, are needed to overcome interindividual variability.

![](images/variability.png)

### Sample Collection

- Sample collection is also a difficult challenge and highly dependent on the study question.
- The microbiome can change in an individual over time, especially in diseases marked by flare ups like IBD.
- Samples might not be representative of the site in question. For example, a stool sample sits in the rectum – an environment that is undergoing dehydration and fermentation which might select for different bacteria than in the small intestine.

![](images/sample-collection.png)

### Confounding Factors

- When conducting a clinical experiment, it is pertinent to stratify accounting for age, gender, diet, etc.
- Sampling over time is incredibly valuable as you can better capture intrapatient variability.
- Additionally, the way the sample is processed can also confound your results

![](images/confounding.jpg)

## What is an Amplicon?

- Microbiome Amplicon sequencing involves sequencing a specific gene from microbial community

![](images/amplicon.png)

### Why Sequence One Gene?

- Genes can vary per organism and may not be well conserved across species. To assess the microbial community composition, we need to sequence a conserved gene across organisms of interest:
    - **16S ribosome DNA (rDNA)** for prokaryotes
    - **18S rDNA** and internal transcribed spacers (ITS) for eukaryotes 

- In the selected gene there are different levels of conservation across organisms. To circumvent this parts of the gene with high conservation (like the V4 region of 16S rRNA) are selected for

![](images/v4-region.png)

## DADA2

We will go through 16S Metabarcoding analysis using DADA2, however, other methods are available as well:

- [USEARCH](https://drive5.com/usearch/)
- [Mothur](https://mothur.org/)
- [QIIME](http://qiime.org/)