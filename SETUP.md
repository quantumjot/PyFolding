
# Installation on a Mac or PC

We have included links to tutorial videos for setting up on both PC or Mac on the PyFolding
Wiki Page (https://github.com/quantumjot/PyFolding/wiki).

The simplest way to automatically install all the packages that PyFolding requires
and Jupyter notebooks is to download and install "Anaconda":

https://www.anaconda.com/download/

+ You need to install the version with Python 2.7
+ *When installing on a PC, you will need to set up your environmental variables*
  - please see video tutorial linked from the Pyfolding WIKI page.


Once you have installed Anaconda you can then install PyFolding as follows:

#### Setup PyFolding on a Mac
+ Open the Terminal.app and use the following commands:

```
  git clone https://github.com/quantumjot/PyFolding.git
  cd PyFolding
  python setup.py install
```

#### Setup PyFolding on a PC
+ Go to https://github.com/quantumjot/PyFolding
+ Download PyFolding as a zip file
+ unzip where you want the PyFolding to reside
+ Open the command line prompt (CMD.exe) and navigate to where you put PyFolding
+ Then Use the following commands:

```
  cd PyFolding
  python setup.py install
```

---


## Using Jupyter Notebooks with PyFolding

This is the easiest way to use PyFolding. This will enable the less script/computer orientated user,
to use the various scripts we have prepared.  The user can then simply change the data paths / initial variables, etc
for their proteins and re-run the jupyter notebook to obtain their analyse.


#### Opening Up Jupyter Notebooks with Anaconda

+ Once Anaconda and PyFolding are installed, there are a number of ways to start and
use PyFolding as part of the jupyter notebooks.

+ The easiest are:
    - On a PC either open jupyter notebook desktop app or open the command line prompt (CMD.exe)
      then type "jupyter notebook".
    - On a Mac open a Terminal (Terminal.app) then type "jupyter notebook".
    - All will open your default web browser with the jupyter notebook interface.

+ You can then navigate from the jupyter interface on your web browser to where you installed PyFolding
(the notebooks folder) and open the jupyter notebooks we have provided by simply double clicking on them.

+ You can also open a new notebook and start from scratch too!

+ All scripts were written on a Mac, therefore on a PC you will need to change the path of your data
to include "C:\" at the beginning.

+ To run the notebook select "kernel" from the menu and use the "Restart and run all" command.
Once you have done this you can re-run individual cells by selecting the cell and
pressing SHIFT + RETURN keys together at the bottom of the cell.


### Testing Jupyter Notebook & Pyfolding installation

+ Open the jupyter notebook script named "PyFolding-Test.ipynb" from the notebook folder.
+ Select kernel from the menu and use the "Restart and run all" command

Upon executing the script, the output generated should be the same as in the README.md file

i.e

```
==================================================
 Fitting results
==================================================
ID: Simulated protein
Model: TwoStateEquilibrium
Method: scipy.optimize.curve_fit

m: 1.49978 ± 0.00000
d50: 5.00013 ± 0.00000
--------------------------------------------------
==================================================
==================================================
 Fitting results
==================================================
ID: Simulated protein
Model: TwoStateChevron
Method: scipy.optimize.curve_fit

kf: 100.00017 ± 0.00008
mf: 1.00002 ± 0.00000
ku: 0.00500 ± 0.00000
mu: 0.99992 ± 0.00000
--------------------------------------------------
==================================================
Test completed!
```

---
