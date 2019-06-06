# Bash Executables for Analysing Ben's Data

This documentation may be out of date

Bash scripts that make use of python executables to gather descriptor data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The following are need to run the python scripts.

```
Python 3.6.7
OVITO 3.0.0 (python packages)
```

The following packages are needed to ensure code functionality:

```
absl-py==0.7.1
apt-xapian-index==0.47
asn1crypto==0.24.0
astor==0.7.1
attrs==19.1.0
backcall==0.1.0
bleach==3.1.0
certifi==2018.1.18
chardet==3.0.4
command-not-found==0.3
cryptography==2.1.4
cupshelpers==1.0
cycler==0.10.0
decorator==4.4.0
defusedxml==0.5.0
distro-info===0.18ubuntu0.18.04.1
entrypoints==0.3
gast==0.2.2
grpcio==1.20.0
h5py==2.9.0
httplib2==0.9.2
idna==2.6
ipykernel==5.1.0
ipython==7.4.0
ipython-genutils==0.2.0
ipywidgets==7.4.2
jedi==0.13.3
Jinja2==2.10.1
jsonschema==3.0.1
jupyter==1.0.0
jupyter-client==5.2.4
jupyter-console==6.0.0
jupyter-core==4.4.0
Keras==2.2.4
Keras-Applications==1.0.7
Keras-Preprocessing==1.0.9
keyring==10.6.0
keyrings.alt==3.0
kiwisolver==1.0.1
language-selector==0.1
m64py==0.2.5
Markdown==3.1
MarkupSafe==1.1.1
matplotlib==3.0.2
mistune==0.8.4
mock==2.0.0
monty==2.0.4
mpmath==1.1.0
nbconvert==5.4.1
nbformat==4.4.0
netifaces==0.10.4
networkx==2.3
notebook==5.7.8
numpy==1.16.1
olefile==0.45.1
palettable==3.1.1
pandas==0.24.2
pandocfilters==1.4.2
parso==0.4.0
pbr==5.1.3
pexpect==4.2.1
pickleshare==0.7.5
Pillow==5.1.0
prometheus-client==0.6.0
prompt-toolkit==2.0.9
protobuf==3.7.1
ptyprocess==0.6.0
pycairo==1.16.2
pycodestyle==2.3.1
pycrypto==2.6.1
pycups==1.9.73
PyDispatcher==2.0.5
Pygments==2.3.1
pygobject==3.26.1
pymatgen==2019.4.11
pyparsing==2.3.1
PyQt5==5.12
PyQt5-sip==4.19.14
pyrsistent==0.14.11
PySDL2==0.9.6
python-apt==1.6.3+ubuntu1
python-dateutil==2.8.0
python-debian==0.1.32
pytz==2019.1
pyxdg==0.25
PyYAML==3.12
pyzmq==18.0.1
qtconsole==4.4.3
reportlab==3.4.0
requests==2.18.4
requests-unixsocket==0.1.5
ruamel.yaml==0.15.94
scikit-learn==0.20.3
scipy==1.2.1
SecretStorage==2.3.1
Send2Trash==1.5.0
six==1.11.0
sklearn==0.0
spglib==1.12.2.post0
ssh-import-id==5.7
sympy==1.4
systemd-python==234
tabulate==0.8.3
tensorboard==1.13.1
tensorflow==1.13.1
tensorflow-estimator==1.13.0
termcolor==1.1.0
terminado==0.8.2
testpath==0.4.2
tornado==6.0.2
traitlets==4.3.2
ubuntu-drivers-common==0.0.0
ufw==0.36
unattended-upgrades==0.1
urllib3==1.22
virtualenv==15.1.0
wcwidth==0.1.7
webencodings==0.5.1
Werkzeug==0.15.2
widgetsnbextension==3.4.2
xkit==0.0.0
xlrd==1.2.0
```

### Installing

Make sure to include this folder and the python folder in the environment variable PATH.
If on a linux machine, include the following to either a .bashrc or .zshrc file in the home directory:

# Path to analysis_scripts

```
export PATH="/home/nerve/Tools/analysis_scripts":$PATH
export PATH="/home/nerve/Tools/analysis_scripts/python_rc":$PATH
export PATH="/home/nerve/Tools/analysis_scripts/python_diffusion":$PATH
export PATH="/home/nerve/Tools/analysis_scripts/python_crystal":$PATH
export PATH="/home/nerve/Tools/analysis_scripts/python_misc":$PATH
```

The wrapper for OVITO is needed aswell. Make sure to include OVITO's python scripts in $PYTHONPATH.

# Path to OVITO python wrapper

```
export PYTHONPATH="/home/nerve/Tools/ovito/lib/ovito/plugins/python":$PYTHONPATH
```

A sample use of the bash scripts is the following for Ben's Rc runs:

```
cd data
calculate_tg
calculate_apd
calculate_variance
gather_tg
gather_apd
calculate_mean_tg
calculate_mean_apd
gather_enthalpy_glass
calculate_mean_enthalpy_glass
plot_composition_tg
```

A sample use of the bash scripts is the following for Ben's diffusion runs:

```
calculate_diffusion
gather_diffusion
calculate_ico_at_tg
gather_ico_at_tg
calculate_ico_at_tlow
gather_ico_at_tlow
calculate_msd
```

A sample use of the bash scripts is the following for Ben's crystal enthalpy runs:

```
gather_enthalpy_crystal
```

Each of the commands calculates a different set of data. All data is stored in a directory called export. Each of the bash scripts in commented for python tool use. The order listed ensure functionality because some scripts need data produced from other scripts to function. Furthermore, some of these scripts are only applicable for a specific set of runs.

## Coding Style

Python scripts follow PEP 8 guidelines. A usefull tool to use to check a coding style is pycodestyle.

```
pycodestyle script.py
```

## Authors

* **Lane Schultz** - *Initial work* - [leschultz](https://github.com/leschultz)

## Acknowledgments

* The Computational Materials Group (CMG) at the University of Wisconsin - Madison
* Dr. Dane Morgan for computational material science guidence
* Dr. Izabela Szlufarska for computational material science guidence
* Benjamin Afflerbach (data)
