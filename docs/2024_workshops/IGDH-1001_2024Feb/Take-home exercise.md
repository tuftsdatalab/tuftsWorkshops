Date: Feb 29, 2024   
Author: Shirley Li, xue.li37@tufts.edu     
Class ID: Sp24-IDGH-1001-1-Bioinformatics    
Canvas link: https://canvas.tufts.edu/courses/55751/assignments

# Exercise 1 
## 1.1
Visit the NCBI to search for the BioProject [PRJNA891065](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA891065) and anwser the following questions:           
>1. What is this project about? Summarize in a few sentences.     
>2. How many biosamples are included in the project?
>3. Which gene is being sequenced, the 16S or 18S rRNA gene?

## 1.2
The paper for the bioproject is available here. [Combining omics tools for the characterization of the microbiota of diverse vinegars obtained by submerged culture: 16S rRNA amplicon sequencing and MALDI-TOF MS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9767973/)

>Summarize the abstract of the paper in no more than three sentences.

## 1.3
For this SRA run [SRR21926282](https://www.ncbi.nlm.nih.gov/sra/?term=SRR21926282), answer the following questions:     
>1. What sequencing instrument was used?
>2. What sequencing strategy was applied?
>3. Are the reads paired-end or single-end?
>4. How many reads were generated in this sequencing run?
>5. What is the biosample for this run? Which sample is it? Check the "environmental medium" on the biosample page.

# Exercise 2 
Review the materials on "[Hands-on sesssion.md](https://github.com/shirleyxueli41/Tufts_workshops/blob/main/IGDH-1001_2024Feb/Hands-on%20session.md)", especially the session "Exercise 3 Taxonomy assignment and interpretation." and "Exercise 4 Taxonomy visualization.". 

For the samples provided, please follow these steps: First, identify the Sequence Read Archive (SRA) run numbers associated with each sample. Next, utilize the Galaxy platform to process these SRA runs. Your goal is to ascertain the top three genera present in each sample. Additionally, for each sample, create a Krona chart to visually represent the taxonomic classification of the organisms found.

>SAMN31308859      
>SAMN31308863     
>SAMN31308856     
>SAMN31308853      
     
> [!WARNING]       
> *Identify the top three genera instead of species.*               
> *To do so, for the filtering step in Hands-on exercise.md , replace "c4==S" to "c4==G"*
> S stands for species, G stands for genera. 
           
![Screenshot 2024-02-29 at 12 42 18](https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/83959223-4da5-4ea2-b5fb-f7c7c50f2e1b)


Analyze and compare the results from the four SRA runs, particularly between samples from wineries and breweries.
Write a report including the answers to the questions above, results from Galaxy, and screenshots, in a Microsoft Word document.

# Assignment - Rubric (25 points total)               
1 Exercise 1.1 (1 point for each question, 3 points in total)      
2 Exercise 1.2 (2 points)     
3 Exercise 1.3 (1 point for each question, 5 points in total)      
4 Exercise 2:        
* Identify the SRA ID (1 point for each sample, 4 points in total)     
* Process the SRA in Galaxy and identify the top three genera for each sample. (1 point for each sample, 4 points in total)
* Visualize the results with Krona plots (1 point for each sample, 4 points in total)
* Summarize the results, with a comparison between winery and brewery samples. Remember to compare between winery samples and brewery samples. (3 points)

     
> [!NOTE]       
> *Please don't hesitate to reach out if you have any questions.*        


   
