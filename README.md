# PyFolding


PyFolding is a simple Python based framework for fitting kinetic and thermodynamic models to 
protein folding data.  It includes several basic models and is extensible to enable fitting
of more interesting models.

---

### Requirements

Pyfolding has been tested on Python 2.7, and requires the following additional packages:
+ Numpy
+ Scipy
+ Matplotlib
+ (Optional) Jupyter

### Installation

You can install Pyfolding by cloning the repo and running the setup script:
```sh
$ git clone https://github.com/quantumjot/PyFolding.git
$ cd PyFolding
$ python setup.py install
```

### Testing the installation

Once installed, a simple script will execute the test function. This tests
whether the package has installed correctly.  The following script is found in /tests/test_installation.py

```python
# import the pyfolding package
import pyfolding
# run a self test to check everything is installed
pyfolding.test()
```

Execute the script as follows:
```sh
$ cd tests
$ python test_installation.py
```

Upon executing the script, the following output should be generated:

```sh
--------------------
 Fitting results
--------------------
Model: TwoStateEquilibrium
alpha_f: 0.00064 ± 0.00036
beta_f: -0.00065 ± 0.00030
alpha_u: 1.03272 ± 0.00184
beta_u: -0.00333 ± 0.00028
m: 1.48006 ± 0.00300
d50: 5.01164 ± 0.00177
--------------------
--------------------
 Fitting results
--------------------
Model: TwoStateChevron
kf: 99.08423 ± 0.31144
mf: 0.99802 ± 0.00127
ku: 0.00544 ± 0.00005
mu: 0.98653 ± 0.00123
--------------------

Test Completed!
```

### Data format

Raw data for PyFolding should be provided in .csv files.

---

### Current models supported

+ Two-state equilibrium unfolding
+ Three-state equilibrium unfolding
+ Homopoylmer Ising model for equilibrium unfolding
+ Heteropolymer Ising model for equilibrium unfolding
+ Two-state chevron
+ Two-state with moving transition state chevron
+ Three-state chevron with fast pre-equilibrium chevron 
+ Three-state chevron with fast phase chevron
+ Three-state sequential barriers chevron
+ Parallel two-state chevrons

---
### Example scripts

There are numerous example scripts in the examples folder. A simple script to
fit the equilbrium denaturation and kinetic chevron to two-state models is 
shown below:

```python
import pyfolding
import models

# load the kinetic and equilibrium folding data
chevron = pyfolding.read_kinetic_data("WT_chevron.csv")
equilibrium = folding.read_equilibrium_data("WT_equilibrium.csv")

# fit the equilibrium data to a two-state model
equilibrium.fit_func = models.TwoStateEquilibrium
equilibrium.fit()

# print the m-value, transition midpoint and stability
print equilibrium.m_value
print equilibrium.midpoint
print equilibrium.deltaG

# use the midpoint (D_50) of the equilibrium curve as the kinetic midpoint
chevron.midpoint = equilibrium.midpoint

# now fit the chevron to a two-state model
chevron.fit_func = models.TwoStateChevron
chevron.fit()

# plot the output
folding.plot_figure(equilibrium, chevron)
```

---
### IPython Notebooks

PyFolding can be used within IPython Notebooks. Several example notebooks are provided in the examples folder.