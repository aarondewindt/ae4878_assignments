# AE4-878 Mission Geometry and Orbit Design assignments
This repository contains my assignments for the TU Delft AE4-878 course.

## Repository structure
The `ae4878_assignments` folder is a python package containing the code 
and reports for the assignments. It is split into modules, one for each
individual assignment.

Each assignment contains one or more python scripts with the functions 
and/or classes used. Unittests can also be found, these are the files 
prefixed with `test_`. 
The main entry point for each assignment is the Jupyter notebook 
(`.ipynb` file) named after the assignment title. For convenience a PDF 
printout of the notebook is provided, but it's necessary to use Jupyter
in order to run the code in the notebook.
A pdf of the report may also be found in the assignment module.

## Requirements
 - python 3.8 or higher
 - numpy
 - scipy
 - tabulate

## Running using docker
To ensure a consistent development environment a docker container is used 
to host the python interpreter and JupyterLab server used to develop and 
run the code. To start the container first ensure that both docker and 
docker-compose are installed on the host machine. Starting the JupyterLab 
instance can be done by running the following command in a terminal from the 
repository root directory (the directory containing the `docker-compose.yml` file).

`docker-compose up`

Instructions on how to open JupyterLab will appear in the terminal after the 
container has successfully started. This JupyterLab instance will have all 
requirements preinstalled and the assignment code should run without issues.
