Simulation and Figure rendering files for "Accelerated simulation of evolutionary trajectories in origin--fixation models" by Ashley I. Teufel and Claus O. Wilke
====================================================

This repository contains the simulation scripts written in Python as well as the figure rendering scripts written in R for the above-mentioned manuscript. This manuscript introduces a new method to simulate evolution under the origin--fixation model. We applied this method to simulations of protein evolution using an all-atom description of protein structure. In order to compare our accelerated method to the traditional one, the procedure by which Rosetta makes a mutation to a protein was chosen for its speed rather than its relationship to biological reality. Hence, if you are interested in accurately simulating protein evolution under realistic biophysical constraints modification to the simulation scripts is necessary.   

# Dependencies

To run the simulation scripts you must have python 2.7, Rosetta, and Pyrosetta installed. Pyrosetta must be included in the path that you are running the scripts from.  Most of the figure scripts depend on the R libraries cowplot, ggplot2, and grid. Figures of directed graphs also need qgraph and igraph.  
