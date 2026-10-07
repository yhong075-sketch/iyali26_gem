# Comparative Flux Sampling Analysis (CFSA)

Extension of the Cobrapy flux sampling functionality for constraint-based models.

Includes extras such as:
* Comparison of different scenarios (growth vs production)
* Filtering reactions based on different criteria
* Parallel execution of multiple scenarios
* Caching
* Convergence testing

### Installation:
* The easiest way to install is via 'pip':
    `pip install git+https://gitlab.com/wurssb/Modelling/sampling-tools.git#egg=sampling-tools`
* Alternatively:
    1. Install the dependencies below via your preferred method.
    2. Clone this repository: `git clone https://gitlab.com/wurssb/Modelling/sampling-tools`
    3. Run `python setup.py install` in the cloned repository folder. 

* Reaction filtering methods are provided as Jupyter notebooks. These can be installed using conda or Anaconda or using 'pip' (https://jupyter.org/install)

### Usage:
Examples are included demonstrating how to:
1. **Generate the samples.** Python scripts _example_sampling_coleaginosus.py_, _example_sampling_yeast.py_
    - User can change the general settings including:
        - Working directories
        - Model settings (e.g. target_id, exchange reactions...)
        - Sampling settings (e.g. number of samples, optimality, flux fraction...)
    - These scripts results in the generation of sampling_results.xlsx files that 
    contain un-filtered sampling results.
    - It also generate sample files that can be used for plotting distribution graphs
2. **Reaction filtering**: Jupyter notebooks _example_filtering_coleaginosus.ipynb_, _example_filtering_yeast.ipynb_
    - Takes as input the sampling_results.xlsx files and filter reaction targets
    based un multiple criteria
    - Generates filtered_results.xlsx
3.** Generation of distribution graphs: ** The python scripts_ example_distributions_coleaginosus.py_, _example_distributions_yeast.py_
    - Use as input the filtered_results.xlsx file and the generated samples
    - Generate KD.pdf and OV.pdf wtith distribution graphs for down-regulation
    and over-expression targets.


### Dependencies:
* python (3.7+)
* pandas (1+)
* cobra (available on the conda forge channel when using conda)
* numpy
* scipy
* seaborn
* h5py
* openpyxl
* pytables

### License:
* MultiLevelPool (`sampling_tools/multilevelpool.py`) is licensed under CC BY-SA 4.0, see the file header for details.
* Everything else is licensed under the MIT license (see `LICENSE`)
