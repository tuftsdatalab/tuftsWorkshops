
# Tufts Bioinformatics Resources

## Introduction
Welcome to the Tufts Bioinformatics Resources page! This resource is designed to help researchers, students, and faculty members navigate the bioinformatics tools, workshops, and resources available at Tufts. Whether you're just getting started or you're a seasoned bioinformatician, we hope this guide will be useful for your work.

## Mailing List (E-list)
To stay updated on bioinformatics education, software tools, and workshop notifications, subscribe to our e-list: [best@elist.tufts.edu](https://elist.tufts.edu/sympa/subscribe/best?previous_action=info).

In the future, we will post workshop notifications and resources through this email list.

## Workshops and Resources

### Current Workshops
- [Bioinformatics workshops by TTS Research Technology (2022-2024)](https://tuftsdatalab.github.io/tuftsWorkshops/)
  
### Archived Workshops
- [TTS Datalab Archived Bioinformatics Workshops](https://tuftsdatalab.github.io/Research_Technology_Bioinformatics/)
  
### External Resources
- [Bioinformatics Education and Services at Tufts](https://best-tufts.github.io/bioinformatics_workshops/)
- [Tufts Bioinformatics Website (Last updated July 2023)](https://it.tufts.edu/bioinformatics)

*Note: Some materials may be outdated, so always verify the relevance of content before applying it to your projects.*

## Tools on the Cluster

Use `module avail` to check the full list of tools available on the cluster. Below is a categorized list of some commonly used tools:

### Genomic Tools
- `kallisto/0.48.0`
- `kraken2/2.1.3`
- `bcftools/1.17`

### RNA-Seq Tools
- `salmon/1.5.2`
- `STAR/2.7.9a`
- `htseq/0.13.5`

### Phylogenetics
- `beast2/2.6.6`
- `raxml/8.2.12`

### Statistical Tools
- `R/4.4.0`
- `qiime2/2024.2`
- `mothur/1.46.0`

### A Few Tips:
1. Before installing your own tools, check if they are already available on the cluster using the `module avail` command.
2. Always be aware of the software versions, especially when using scripts from colleagues.
3. For less common tools, consider installing them yourself to ensure you have full control over the version and availability.

## How to Install Your Own Tools
We will cover how to install tools from source code in our upcoming workshop. If you need to install a tool not commonly used, it's best to do it yourself to avoid issues with maintenance. Stay tuned for a detailed guide!

## Cluster Usage Tips
- **Version Control:** Keep an eye on software versions, especially for tools like R. For example, R/4.4.0 is available, but many researchers may still use older versions like R/3.X.X.
- **Environment Management:** Use `conda` or `virtualenv` to manage your tools and environments effectively.
- **Job Submission Examples:** Check out our upcoming guide for submitting simple and complex job scripts using Slurm.

## Contact Us
We are part of the TTS Research Technology team. Our services are not limited to bioinformatics; we also support data science, statistics, and research data management.

For consultations, please submit a ticket to [tts-research@tufts.edu](mailto:tts-research@tufts.edu).

For more details, visit our page: [Research Technology Consulting Services](https://it.tufts.edu/research-technology/research-technology-consulting-services).
