# QuaDS

Qualitative/quantitative Descriptive Statistics is a Python implementation of the [catdes()](http://factominer.free.fr/factomethods/description-des-modalites.html) function from the [FactoMiner R package](http://factominer.free.fr) with extras.

## Installation
In order to use the pipeline you first of all have to `clone` the git
repository or download it.

The pipeline has been written for Python 3.9. 
You have to create a Python environment to run the pipeline.

    conda create -n environment_quads
    conda activate environment_quads
    conda install python=3.9

It relies on several libraries which are listed in the `requirement.txt` file.
Alternatively the dependencies can by installed using pip:

    python -m pip install -r requirements.txt

## Usage
You have to put your data files in the case : /data
You have to say in the config_file.yml the different parameters :
  - directory of data and results
  - names of datafile, and output files
  - list of qualitative variable names
  - list of quantitative variable names
  - cluster variable
  - different thresholds of the tests
  - colors of the visual output

descriptives statistics :

    python -m scripts.launch_quads
    
## Visuals
When you have your tables support and you want to see the visualisation

    python -m scripts.visualisation


## Deactivation of conda
You have finish to use the pipeline.

    conda deactivate
