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

descriptives statistics :

for gower distance:

    python -m scripts.gower5
    python -m scripts.gower6
    python -m scripts.gower7
    
for semantic distance :

    python -m scripts.semantic5
    python -m scripts.semantic6
    python -m scripts.semantic7
    
for genetic groups :

    python -m scripts.genetic_groups
    
## Visuals
When you have your tables support and you want to see the visualisation

for gower distance :

    python -m scripts.visu_gower5
    python -m scripts.visu_gower6
    python -m scripts.visu_gower7

for semantic distance :

    python -m scripts.visu_semantic5
    python -m scripts.visu_semantic6
    python -m scripts.visu_semantic7
    
for genetic groups :

    python -m scripts.visu_genetic_groups


## Deactivation of conda
You have finish to use the pipeline.

    conda deactivate
