## Background

Sequencing data analysis typically focuses on either assessing DNA or RNA. As a reminder here is the interplay between DNA, RNA, and protein:

![](images/dna-rna.png)

### DNA Sequencing

- Fixed copy of a gene per cell 
- Analysis goal: Variant calling and interpretation

### RNA Sequencing

- Copy of a transcript per cell depends on gene expression
- Analysis goal: Differential expression and interpretation

!!! note
    Here we are working with DNA sequencing
    
## Next Generation Sequencing

Here we will analyze a DNA sequence using next generation sequencing data. Here are the steps to get that data:

- **Library Preparation:** DNA is fragmented and adapters are added to these fragments

![](images/lib-prep.png)

- **Cluster Amplification:** This library is loaded onto a flow cell, where the adapters help hybridize the fragments to the flow cell. Each fragment is then amplified to form a clonal cluster

![](images/cluster-amp.png)

- **Sequencing:** Fluorescently labelled nucleotides are added to this flow cell and each time a base in the fragment bonds a light signal is emmitted telling the sequencer which base is which in the sequence.

![](images/sequencing.png)

- **Alignment & Data Analysis:** These sequenced fragments, or **reads**, can then be aligned to a reference sequence to determine differences.

![](images/alignment-data-analysis.png)


## Singe End v. Paired End Data

- **single-end** sequence each DNA fragement from one end only
- **paired-end** sequence each DNA fragement from both sides. Paired-end data is useful when sequencing highly repetitive sequences.
        
![](images/single-paired.png)

## Variant Calling

![](images/variant-overview.png)

## Ploidy 

- When discussing variant calling it is worth mentioning an organism's **ploidy**. Ploidy is the number of copies of each chromosomes.

    - Humans cells are diploid for autosomal chromosome and haploid for sex chromosomes
    - Bacteria are haploid
    - Viruses and Yeast can by haploid or diploid

![](images/ploidy.png)

- Variant callers can use ploidy to improve specificity (avoid false positives) because there are expected variant frequencies, e.g. for a diploid:

    - **Homozygous**
    - both copies contain variant
    - fraction of the reads ~1

    - **Heterozygous**
    - one copy of variant
    - fraction of reads with variant  ~0.5

![](images/het-homo.jpg)

