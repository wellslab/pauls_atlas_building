# Basic repository template

Primarily designed for new students and those new to version control, this repository can work as a basic template to start working on a bioinformatics project within the Wells Lab. The idea is to clone this project, change its name, and use it to store analyses and data, while pushing changes to the repository. If you're already familiar with repositories or are using one already, you can skip the set up instructions and have a look at extra set up steps or recommended folder structures.

## Setup instructions

1. Set up an environment on your computer

    This step is somewhat independent of this repository, so can be skipped, but we generally use a conda environment for our work. An environment is where all your programs are installed, and conda makes it easy to manage this. Having separate environments for separate projects mean you can reproduce your work more easily later on, and have less chance of program versions conflicting across different projects.

    Download and install [Anaconda](https://www.anaconda.com/products/individual) (click on Download button and select your computer). More details can be found under howtos/ManagingEnvironment.md.

2. Clone this repo and make it yours

    You can either download the .zip file of this repo and unzip it or use the git clone command. Note that the location of this project folder is independent of where your environments are. Now change the name of the folder to something meaningful (in the example below we're using 'dc_atlas'), remove the .git subfolder inside the folder, then run git init to re-create a fresh one. On the command line you would run these for all these steps:
    ```
    > git clone git://github.com/wellslab/basicprojecttemplate
    > mv basicprojecttemplate dc_atlas
    > cd dc_atlas
    > rm -r .git
    > git init
    ```

3. Create a new repo on github and sync with this repo
    
    Login to github and go to github.com/wellslab and create a new repo, calling it the same as in step 2. Then you tell your local repo what your origin will be:
    ```
    > git remote add origin git@github.com:wellslab/dc_atlas.git
    ```
    
That's it! You can now start working in this repo: make changes, commit them, then push them to github. Don't forget to replace this readme with one more appropriate for your own project too.

## Extra setup steps

- For managing large data files in this repo, read howtos/UsefulGitTips.md. You can use git lfs to tell git to track specific files or folders differently to avoid the repository getting too large.

- For handling jupyter notebooks which can be tricky to version control, read howtos/ManagingNotebooks.md. We propose that you track .html and .py files as well as .ipynb files, and jupyter_notebook_config.py has been provided to make this easier.


## Folder structure

You can change the subfolders in this project in any way you want, but for beginners, the current structure is recommended to deal with most situations in various projects. This template comes with example folders and files to get you started.

* howtos/

    Contains notes you can keep as you discover tricks and work out how to do tasks (markdown format makes it easier to read on the website directly).

* notebooks/

    Jupyter notebooks are placed here. If there are many notebooks in a big project, it's useful to create subfolders here too (eg. 0_Initial_analysis/, 1_Clustering/). For tips on how to version control jupyter notebooks, see howtos/ManagingNotebooks.md.

    Output files from notebooks can be placed in subfolders within notebooks/. This actually makes it easy to identify which output files were produced by which notebooks. These output files can then be shared with collaborators.

* received/

    Contains data files which often act as inputs to programs you write. These may be downloaded from websites or received by email. Subfolders here can be used to distinguish the different sources. For example, we may want to have a subfolder called "GEO" here for files downloaded from GEO website, and another called "Paul" for files received from Paul. It's good practise to record where each file came from, as this can easily be forgotten later and it can put much work downstream in doubt when you're not sure where the input files originate.
