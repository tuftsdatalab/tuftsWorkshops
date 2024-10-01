# Software Installation and Environment Management in HPC
For bioinformatics users, learning how to install software on HPC is very useful. 

Bioinformatics workflows often require specialized software tools, many of which have complex dependencies or need to be optimized for the HPC environment. Installing and managing these tools on an HPC system allows bioinformatics users to access the latest versions of key software. 

Understanding software installation on HPC ensures that users can customize their environments, streamline their analyses, and maintain reproducibility in their research, which is critical in bioinformatics given the rapid pace of data generation and evolving computational tools.

In this workshop, we will introduce how to install bioinformatics softwares on HPC using different ways. 



## Agenda

- [Application Installation from Source Codes](01_source.md)
- [R package Installation](02_R.md)
- [Conda Environment Management, Jupyer Kernel and Modules](03_conda.md)
- [Simplying Conda Environment Management with conda-env-mod](04_conda-env-mod.md)

!!! note "Containerization is an easier and recommended way" 
In this workshop, we will focus on R package installation and package installation from source codes or using a Conda environment. Actually, there is an easier way to run applications on clusters without any installation steps; this method involves containerization using Singularity/Apptainer. We will not introduce containerization in this workshop but plan to conduct a container workshop in Spring 2025. If you are interested in containers, you can check out our [Spring, 2024 container workshop](https://zhan4429.github.io/TuftsContainers.github.io/).


<div style="display: flex; align-items: center;">
  <img src="https://docs.sylabs.io/guides/3.8/user-guide/_static/logo.png" alt="singularityCE" style="height:200px; margin-right: 10px;"/>
  <img src="https://apptainer.org/docs/user/main/_static/logo.png" alt="apptainer" style="height:150px;" />
  <img src="https://www.vikingsoftware.com/wp-content/uploads/2024/02/Docker.png" alt="docker" style="height:200px;" />
</div>


## Presenters

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/shirleyxueli41"><img src="https://avatars.githubusercontent.com/u/88347911?v=4" width="100px;" alt=""/><br /><sub><b>Shirley Li</b></sub></a><br /></
    td>
    <td align="center"><a href="https://github.com/zhan4429"><img src="https://avatars.githubusercontent.com/u/90942318" width="100px;" alt=""/><br /><sub><b>Yucheng Zhang</b></sub></a><br /></td>    
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
[Next: fetchngs](01_fetchngs.md)
