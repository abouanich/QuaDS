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
  - factor variable
  - different thresholds of the tests
  - colors of the visual output

descriptives statistics :

    python3 scripts/launch_quads.py

## Outputs 
For your qualitative analysis, you will obtain maximum 4 outputs:
  - Chi2.csv: informs which variable is implicated in the factor's modalitities
  - fisher_exact.csv: informs which variable (with low frequencies) is implacated in the factor's modalities  
  - qualitative_hypergeometric.csv: informs if the variable is implicated in the factor's     modalities, this file informs in each factor's modality if the variable modality is:
    - over-represented
    - under-represented
    - not significant
    - not present
  - weight.csv: informs the qualitative variables contribution to the factor's modalities
  
For your quantitative analysis, you will obtain maximum 5 outputs: 
  - normality.csv: informs if the quantitative variables have a normal distribution in the different factor's modalitities 
  - homoscedasticity.csv: informs if the quantitative variables' standard deviation are the same in the different factor's modalitities 
  - anova.csv : informs if a variable have a significant higer or lower average than the average of all the groups.
  - kruskal_wallis.csv: informs if a variable have a significant higer or lower average than the average of all the groups for variables that are not normal distributed.
  - quantitative_gaussian.csv informs for significative variable (to ANOVA or kruskal wallis) if the average of the variable is:
    - above from the average for all individuals
    - below from the average for all individuals
    - Not significantly different from the average for all individuals


## Visuals
When you have your tables support and you want to see the visualisation

    python3 scripts/visualisation.py


## Deactivation of conda
You have finish to use the pipeline.

    conda deactivate
