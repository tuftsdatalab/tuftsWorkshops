## DADA2 Error Model

- The DADA2 Error Model asks: What is the error rate for an amplicon sequence read i that was produced from a 
sequence j over L aligned nucleotides with a quality score q?

- Basically, this is a product of error probabilities given some quality score e.g. p(A > G, 35)

- This p-value assess if sequence i is too abundant for it to be explained by errors in amplicon sequencing

![](images/error-model.png)

Here we will leverage the parametric model to learn error rates and then plot them:

```R
# Learn Error Rates

## dada2 uses a parametric model to learn the error rates
## for each amplicon data set
errForward <- learnErrors(filtForward)
errReverse <- learnErrors(filtReverse)
plotErrors(errForward)
```

![](images/error-plot.png)

## Inferring Amplicon Sequence Variants (ASVs)

- So far, we have assigned p-values for each sequence in each sample
- DADA2 then tries to determine which sequences are of biological origin (ASVs) and which aren’t by assessing which sequences are present in other samples
- If a sequence is present in another sample, it is more likely that it is a real biological sequence

![](images/infer-asv.png)

```R
# Sample Inference

## we will now run the dada2 algorithm 
## this algorithm delivers "true" sequence variants
## with information gathered from the error model 
## generated above
dadaForward <- dada(filtForward, err=errForward)
dadaReverse <- dada(filtReverse, err=errReverse)
dadaForward[[1]]
```

## ASVs vs. OTUs

Traditional 16S metagenomic approaches use OTUs or operational taxonomic units instead of ASVs. So why does DADA2 use ASVs? First let's cover what an OTU is:

- Methods that use OTUs, cluster sequences are clustered together by similarity 
- Those sequences with above a 97% identity threshold are clustered into an OTU
- These OTUs are then combined into a consensus sequence and mapped to a reference database to determine which species it is from

![](images/otu.png)

Originally, OTUs were used to mitigate possible sequence errors by clustering similar sequences and getting a consensus sequence. 
However, this method has been found to inflate the number of unique sequences. 
By contrast, ASV analysis derives an error term to assess the possibility of a sequencing error. 
These sequences are then mapped directly to the organism of interest - giving nucleotide resolution. 





