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

    pip install -r requirements.txt

## Usage
In the example test, the datafile is in the repository: /data
You have to say in the config_file.yml the different parameters :
  - directory of data and results
  - names of datafile, and output files
  - separator of your datafile
  - presence of an index in your datafile
  - list of qualitative variables names
  - list of quantitative variables names
  - cluster variable
  - different thresholds of the tests
  - colors of the visual output

descriptives statistics :

    python3 scripts/launch_quads.py

## Outputs 
For your qualitative analysis, you will obtain 2 outputs (datatest explication):
  - x2.csv : informs which variable is implicated in the groups
  - qualitative_hypergeometric.csv informs if the variable is implicated in the groups, 
    this file informs in each groups if the modality is:
    - over-represented
    - under-represented
    - not significant
    - not present

For your quantitative analysis, you will obtain 2 outputs (datatest explication): 
  - anova.csv : informs if a variable have a significant higer or lower average than the average of all the groups.
  - quantitative_gaussian.csv informs if the average of the variable is:
    - above from the average for all individuals
    - below from the average for all individuals
    - Not significantly different from the average for all individuals


## Visuals
When you have your tables support and you want to see the visualisation

    python3 scripts/visualisation.py


## Deactivation of conda
You have finish to use the pipeline.

    conda deactivate
