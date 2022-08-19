# Project Management

Before we dive into R it is worth taking a moment to talk about project management. Often times data analysis is incremental and files build up over time resulting in messy directories:

![](images/messy.png)

Sifting through a non-organized file system can make it difficult to find files, share data/scripts, and identify different versions of scripts. To remedy this, It is reccomended to work within an R Project.

## R Project

To Create a new R project:

1. Go to `File` > `New Project`
2. `New Directory`
3. `New Project`
4. Create a name for your project (e.g. `new-project`)
5. `Create Project`

You will notice that your RStudio console switches to this project directory. When you log out of RStudio you can open this project again by clicking the `.Rproj` file in the project directory. 

!!! note
    The paths will be relative to this project directory as a safe guard against referencing data from outside sources. 

## File Organization

When working on a scientific project it is recommended that you put each project in its own directory and give it a name that is descriptive. Similarly, when naming scripts it is recommended that you also name these scripts after the function they are performing. When it comes to file structure within your project try following this folder structure:

- `doc` : folder for text documents associated with the project
- `data` : folder for your raw data/metadata
- `results` : folder for files generated from data/metadata
- `src` : folder for project's custom scripts/programs
- `bin` : folder for outside programs used in project

## Data Principles

When you deal with data treat it as read-only. Working with data files in something like excel can modify your original data source without any record of what was done to it. That being said, often times you will need to do some data cleaning. When you need to significantly modify your data source make a separate folder withing `data` for the `raw_data` and the `cleaned_data`. Also ensure that the scripts you used to clean the data are placed in a separate folder (e.g. `src/data_cleaning_scripts/`). Data that is generated from this raw data should be deposited in your `results` folder and should be treated as disposable. These files should be reproducible from your raw data using your scripts and are good candidate files to cut if you are getting low on storage.

## Script Management

When performing analyses you'll note that some code blocks are useful in multiple scenarios. It is a good idea to store these reusable chunks in a separate folder to use in other analysis scripts. 
