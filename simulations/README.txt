This folder contains R code and data to reproduce the simulations of Degras (2021) "Scalable feature matching across large data collections". 

Gurobi 
Some the methods used in the simulations utilize integer linear or quadratic program. To be able to run these methods, please install the software Gurobi and its R interface. Free academic licenses can be obtained at https://www.gurobi.com
The R code will run if Gurobi is not installed, but the methods that need it will not be calculated. 

To reproduce the simulations of the paper:
1) Execute the file simulations.R. To avoid possible forking during parallel execution, we recommend executing the file in a Terminal session: open a Terminal, move the simulation folder and type
Rscript simulations.R
2) Execute the file simulations-summary.R to reproduce the tables and figures in the simulation section of the paper. 

The file relaxation.R contains additional experiments on relaxation and lower bounds for the matching problem.

  