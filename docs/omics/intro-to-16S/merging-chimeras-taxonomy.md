## Merging Reads

- For paired-end data there is a good deal of overlap between the forward and reverse read
- To resolve this redundancy, these reads are collapsed into contigs

<figure markdown>
  ![](images/merging1.png){ width="500" }
</figure>


Let's do this with code now!

**Code Chunk 6**

![](images/r-markdown-header.png)

```R
# Merge Read Pairs

## so far we have "denoised", so to speak, 
## these sequence variants. We now need to merge the
## forward and reverse strands
mergers <- mergePairs(
  dadaForward,
  filtForward,
  dadaReverse, 
  filtReverse, 
  verbose=TRUE)
```

## ASV Table

Now that our sequences are merged we can create an ASV counts table, basically telling us how which samples contain which ASV's:

**Code Chunk 7**

![](images/r-markdown-header.png)

```R
# Making a Sequence Table

## now that we have merged sequences we can construct
## an Amplicon Sequence Variant (ASV) table
## 
seqtab <- makeSequenceTable(mergers)
```

## Chimera Removal

- During Sequencing microbial DNA is subjected to PCR to amplify DNA
- During PCR it is possible for two unrelated templates to form a non-biological hybrid sequence
- DADA2 finds these chimeras by aligning each sequence to more abundant sequences and seeing if there are any low abundant sequences that can be created by  mixing the left and and right sides of the more abundant sequences

<figure markdown>
  ![](images/chimera1.jpeg){ width="500" }
</figure>

Now in code:

**Code Chunk 8**

![](images/r-markdown-header.png)

```R
# Removing Chimeras

## Chimeric sequences occur as errors during PCR 
## when two unrelated templates for a hybrid sequence
## we will need to remove them before going forward

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE)
```

Now let's check if any chimeric sequences are removed:

**Code Chunk 9**

![](images/r-markdown-header.png)

```R
## check to see if the dimensions are different
## between the chimera filtered and unfiltered
## ASV tables

dim(seqtab)
dim(seqtab.nochim)
```

```
> dim(seqtab)
[1]  8 27
> dim(seqtab.nochim)
[1]  8 27
```

!!! info 
    We can see here no chimeric sequences were removed because our before and after sequence count matrices have the same dimensions.

## Pipeline Quality Control 

We will also take a moment to do some final QC:

**Code Chunk 10**

![](images/r-markdown-header.png)

```R
# Final QC

## we have performed quite a few steps 
## and it would be nice to get a final qc check 
## before assigning taxonomy
getN <- function(x) sum(getUniques(x))
finalQC <- cbind(
  out, 
  sapply(dadaForward, getN),
  sapply(dadaReverse, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim))
colnames(finalQC) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(finalQC) <- sampleNames
finalQC
```

```
           input filtered denoisedF denoisedR merged nonchim
SRR5690809  1000      929       788       743    528     528
SRR5690810  1000      937       766       708    293     293
SRR5690811  1000      930       807       735    360     360
SRR5690812  1000      915       768       756    437     437
SRR5690819  1000      943       761       708    337     337
SRR5690820  1000      911       769       683    349     349
SRR5690821  1000      925       811       700    349     349
SRR5690822  1000      940       802       669    300     300
```

!!! info
    Here we see that we start with 1000 sequences per sample, end up with around 700 after filtering, around 800 after denoising to 
    find unique sequences, and around 300-500 sequences after merging sequences and removing chimeric sequences.

## Assigning Taxonomy

- To determine which taxon each  ASV belongs to DADA2 uses a naïve bayes classifier 
- This classifier uses a set of reference sequences with known taxonomy, here we use the SILVA database, as the training set and and outputs taxonomic assignments with bootstrapped confidence

**Code Chunk 11**

![](images/r-markdown-header.png)

```R
# Assigning Taxonomy

## dada2 uses a naive Bayes classifier when
## assigning taxonomy. This means we need a training
## set of sequences with known taxonomy information.
## here we use the silva database

taxa <- assignTaxonomy(seqtab.nochim, "/cluster/tufts/bio/data/metagenomes/silva/silva_nr99_v138.1_train_set.fa.gz")
```

## Databases

While we use the SILVA database here, there are other options databases:

- [EzBioCloud](https://help.ezbiocloud.net/ezbiocloud-16s-database/)
- [Greengenes](https://greengenes.secondgenome.com)

[Park SC, Won S. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6440677/) performed a benchmarking study assessing these three databases and found that Greengenes predicted fewer correct genera and that SILVA and EzBioCloud predicted roughly similar correct genera. However, SILVA had more false postitives. The authors note that this is probably due to the size of SILVA since it was roughly 4 times the size of EzBioCloud in 2018. 
