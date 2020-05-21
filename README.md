# **What is ESD_thermotrace?**

A [jupyter notebook](https://jupyter.org/) that helps interpreting detrital thermochronometric ages.
Briefly, throughout the notebook one can do the following:

* Interpolate a surface bedrock age map over your study area.
* Forward model the expected detrital age distribution for the desired catchment.
* Load mineral fertility and erosion maps to test their effect on the age distribution.
* Evaluate the confidence of each scenario to your observed detrital age distribution

# **Ok, great, but what *exactly* does it do?**

The *jupyter notebook* hosts all the code, which is well-commented and allows anyone to understand all the operations.

# **What do I need to install?**
Perhaps you already have a Python 3 installation with jupyter,
but ESD_thermotrace might use some libraries you still have not installed.
To ensure a painless installation, this is the way to go:

### Make sure you have an Anaconda installation with python 3.
### If not, visit https://docs.anaconda.com/anaconda/install/ and follow the instructions
### Clone or download this repository to your preferred directory
### Open a terminal (MacOS, Linux) or Anaconda-Prompt window (Windows).
### Navigate to the *ESD_thermotrace* directory.
### Create a new environment from the provided .yml file by entering the following command:
```
conda env create -f ESD_thermotrace_environment.yml
```
### Activate the environment like so
```
source activate ESD_thermotrace
```
or (depending on Anaconda version and operating system)
```
conda activate ESD_thermotrace
```
### launch jupyter
```
jupyter notebook
```
### Open the notebook *ESD_thermotrace.ipynb* in the browser Windows
Press **SHIFT+ENTER** to run the current cell
To close the program, enter CTRL+C in the terminal, or press the Quit button in the User Interface.
Please refer to the [jupyter documentation](https://jupyter-notebook.readthedocs.io/en/stable/) for all other questions on how to use notebooks

### To return to the base Anaconda environment, enter the activation command without specifying the env. name:
```
source activate
```
or
```
conda activate
```
