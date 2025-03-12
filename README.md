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
    conda install r-base

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

For your qualitative analysis, you may obtain up to 4 output files:
  - Chi2.csv: contains the results of the Chi-squared test of independence which assess whether a variable is dependent on the factor.
  - fisher_exact.csv: an alternative to the Chi-squared when the expected frequencies in the contingency table between the variable and the factor are very low (<5). It also tests the dependency between the variable and the factor.
  - qualitative_results.csv: for dependent variables, this file describes the dependency of the factor levels. For each factor level, it indicates whether the variable modality is:
    - over-represented
    - under-represented
    - not significant
    - not present
  - weight.csv: indicates the contribution of the qualitative variables to the factor levels
  
For your quantitative analysis, you may obtain up to 5 outputs files:
  - normality.csv: contains the results of the Shapiro-Wilk test, which assesses whether each quantitative variable meets the normality assumption within the compared factor levels.
  - homoscedasticity.csv: contains the results of Bartlett's test, which verifies the equality of variances between the factor levels.
  - anova.csv : when assumptions are met, this test determines whether there is a significance difference in at least one factor level for each quantitative variable.
  - kruskal_wallis.csv: the non-parametric alternative of anova, used when at least one assumption is not met.
  - quantitative_results.csv: provides information on variables with significant differences (based on anova and Kruskal-Wallis results) and only when the homoscedasticity assumption is verified. If homoscedasticity is not verified, Kruskal-Wallis is applied but, no further description of the factor levels is performed. This file indicates whether the variable mean is:
    - above the overall average
    - below the overall average
    - Not significantly different from the overall average


## Visuals
When you have your tables support and you want to see the visualisation

    python3 scripts/visualisation.py


## Deactivation of conda
You have finish to use the pipeline.

    conda deactivate
