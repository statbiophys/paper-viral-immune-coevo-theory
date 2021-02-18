# Antigenic waves of virus-immune co-evolution


This repository contains the source code associated with the manuscript

Marchi, Laessig, Walczak, Mora: Antigenic waves of virus-immune co-evolution 2021

It allows reproduction of all figures of the manuscript from pre-processed numerical data and also provides the simulation code for transparency and review.


## Figures reproduction code

### Installation requirements

The visualization code uses Python 2.7+.

A number of standard scientific python packages are needed for the numerical simulations and visualizations. An easy way to install all of these is to install a Python distribution such as [Anaconda](https://www.anaconda.com/):

- [numpy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/stable/index.html)


### Files organization/running the code

All files necessary to reproduce the figures is in [figs_paper_coarse_gr_clean](./figs_paper_coarse_gr_clean)
The python scripts are in [python_code](./figs_paper_coarse_gr_clean/python_code), where [python_code](./figs_paper_coarse_gr_clean/python_code)

Every folder contains a file plot.py which needs to be run to produce the figures. For some figures cosmetic changes were done in inkscape as a postprocessing step. In these cases the figures will not be reproduced precisely.

## Simulation and pre-processing code


The simulation and analysis code uses a combination of C, C++, Bash and Python 2.7+. 

The code is posted here for transparency, it is not meant to be run on local machines as it is. It may not run across platforms as some synthax is Linux specific and is not necessarily POSIX compliant. 
Some snippets also rely on the specific folder structure on the machine. It should also be noted that for some parameters choices simulations are computationally heavy and need to be run with the aid of a computing cluster.

The C++ code depends on the packages [FFTW](http://www.fftw.org/) and [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) that need to be installed on the local machine and properly compiled (see Makefile).



### Files organization

