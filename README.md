# PyFolding


PyFolding is a simple Python based framework for fitting kinetic and thermodynamic models to 
protein folding data.  It includes several basic models and is extensible to enable fitting
of more interesting models.



## Requirements

## Installation

You can install Pyfolding by cloning the repo and running the setup script:
```sh
$ git clone https://github.com/quantumjot/PyFolding.git
$ python setup.py -install
```

## Testing the installation

Once installed, a simple script will execute the test function. This tests
whether the package has installed correctly.

```python
# import the pyfolding package
import pyfolding
# run a self test to check everything is installed
pyfolding.test()
```
