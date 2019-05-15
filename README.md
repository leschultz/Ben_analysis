# Bash Executables for Analysing Ben's Data

Bash scripts that make use of python executables to gather descriptor data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The following are need to run the python scripts.

```
Python 3.6.7
Ovito 3 (python packages)
```

### Installing

Make sure to include this folder and the python folder in the environment variable PATH.
If on a linux machine, include the following to either a .bashrc or .zshrc file in the home directory:

```
# Path to analysis_scripts
export PATH="/home/nerve/Tools/analysis_scripts":$PATH
export PATH="/home/nerve/Tools/analysis_scripts/python":$PATH
```

A sample use of the bash scripts is the following:

```
cd data
calculate_tg
calculate_apd
calculate_variety_and_variance
gather_tg
gather_apd
calculate_mean_tg
calculate_mean_apd
plot_composition_tg
```
Each of the commands calculates a different set of data. All data is stored in a directory called export.

## Authors

* **Lane Schultz** - *Initial work* - [leschultz](https://github.com/leschultz)

## Acknowledgments

* Benjamin Afflerbach (data)
