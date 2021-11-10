# Subduction-Dynamics-Tools
This is a repository for codes used in the support of subduction dynamics research. It includes python scripts, jupyter notebooks and matlab scripts. The topics covered range from reproducing analytic solutions, plotting and comparing rheological flow laws, documenting how to define input for model simulations, and more. 

Run Jupyter Notebooks using Binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mibillen/Subduction-Dynamics-Tools/master)

## Contents in the python_scripts folder

### The Affinity-test: run affinity test
This folder contains the scripts and prm file to run affinity test on servers.
The test (i.e. the prm file) is created by Rene Gassmoeller. I got the original bash and python script from Max Rudolph

Author: Haoyuan Li

Updated: 11/09/2021

* AffinityTest.py: the script used to run and analyze the affinity test.
* spherical_shell_expensive_solver.prm: the prm file for the test.
* parse_block_output.sh: a bash script the py script depends on that parses the block outputs from aspect's log file.

#### Instructions on how to run:
The scripts is meant to operate on both server and desktop. Some of the plotting packages might be missing on the server, that why I choose to do analysis on laptop. While on server, I just request a few basic packages to run the test.

A workflow includes:
* run on server and get test results
* download test results to a local laptop
* run on the laptop to analyze the results and plot figures

First run with the '-h' option to see the usage.

    python AffinityTest.py -h

This will give explicit instructions as well as examples of using the script.