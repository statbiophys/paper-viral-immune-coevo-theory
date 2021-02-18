# Antigenic waves of virus-immune co-evolution


This repository contains the source code associated with the manuscript

Marchi, Laessig, Walczak, Mora: Antigenic waves of virus-immune co-evolution 2021

It allows reproduction of all figures of the manuscript from pre-processed numerical data and also provides the simulation code for transparency and review.


## Figures reproduction code

#### Installation requirements

The visualization code uses Python 2.7+.

A number of standard scientific python packages are needed for the numerical simulations and visualizations. An easy way to install all of these is to install a Python distribution such as [Anaconda](https://www.anaconda.com/):

- [numpy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/stable/index.html)


#### Files organization/running the code

All files necessary to reproduce the figures are in [figs_paper_coarse_gr_clean](./figs_paper_coarse_gr_clean) in the required relative positions.
The python scripts are in [python_code](./figs_paper_coarse_gr_clean/python_code), where [lib](./figs_paper_coarse_gr_clean/python_code/lib) contains some plotting cosmetics definitions, and [plots](./figs_paper_coarse_gr_clean/python_code/plots) contains the scripts to run to produce the figures in the corresponding folder. Pre-processed data for Figs 2,3 and 4 can be found in the subdirectories of [fig2](./figs_paper_coarse_gr_clean/fig2).

For some figures cosmetic changes were done in inkscape as a postprocessing step. In these cases the figures will not be reproduced precisely.

## Simulation and pre-processing code


The simulation and analysis code uses a combination of C, C++, Bash and Python 2.7+. 

The C++ code depends on the packages [FFTW](http://www.fftw.org/) and [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) that need to be installed on the local machine and properly compiled (see Makefile).


##### Don't try this at home:
The code is posted here for transparency, it is not meant to be run on local machines as it is. It may not run across platforms as some synthax is Linux specific and is not necessarily POSIX compliant. 
Some snippets also rely on the specific folder structure on the machine. So if the reader wants to try and use this code they should be careful to change the code to reflect their OS and folder structure. 
It should also be noted that for some parameters choices simulations are computationally heavy both in time and memory requirements, and need to be run with the aid of a computing cluster.

#### Files organization

- In  [bash_scripts_handlers_github](./bash_scripts_handlers_github) there are three bash handlers as examples on how we: 1) sweeped through parameters, 2) ran the model C++ code handling extinctions and explosions and output directory structure, 3) performed analysis and plots on the model outputs. These scripts are meant to be run on a  computing cluster mounting a SLURM queuing engine.

