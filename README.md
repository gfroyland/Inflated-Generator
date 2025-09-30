# Inflated-Generator
This repository contains Julia code used for the numerical implementation of the inflated generator, a time-expanded version of the generator of the transfer operator, as described in the papers:

Aleksandar Badza and Gary Froyland. "Identifying the onset and decay of quasi-stationary families of almost-invariant sets
with an application to atmospheric blocking events", <i>Chaos</i> 34(12):123153, 2024.  https://arxiv.org/abs/2407.07278

Aleksandar Badza, Gary Froyland, Roshan J. Samuel and Joerg Schumacher. "Transient almost-invariant sets reveal convective heat transfer patterns in plane-layer Rayleigh-Benard convection", arXiv details TBA.

The code provided here is for three systems:

1. The Switching Double Gyre, an idealised flow system involving two counter-rotating gyres of unequal size within a rectangular flow domain as defined in [1];
2. ECMWF wind velocity data [2] used to identify two atmospheric blocking events occurring during the European Summer of 2003. The two events in question are a "West block", centred around France (further West in Europe) and taking place between 2-12 August 2003; and an "East block", centred around Western Russia (further East in Europe) and taking place between 29 July-3 August 2003; and
3. Three-dimensional Rayleigh-Benard Convection velocity data, used to identify quasi-stationary, almost-invariant plumes of fluid within a convection cell, and also to identify regions within this cell which pertain to less metastable fluid movement.

Included within this repository are function files in Julia used for the numerical construction of the inflated generator, and separate script files used to implement the method on each dynamic system of interest (in the atmospheric case, the code for the West and East blocks is split into separate scripts as two separate spatial domains are used to identify and track each of these blocking events).

# Downloading and Setting Up The Repository

Instructions for setting up the Julia environment and running the inflated generator code for each system of interest are included within the README.md files within each of the three subdirectories present within this repository. To run the Julia code in each case, you'll need to download the relevant subfolder of code/data from this repository and set up the required Julia environment featuring the packages necessary for the successful implementation of this code. The packages required to set up this environment feature within the Project.toml file, while information regarding the versions of these packages used and their subdependencies is contained within Manifest.toml. Each subfolder has its own Julia environment files with packages and subdependencies specific to each system.

Note that if you wish to alternate between each of these three systems you will have to restart Julia. You do not have to restart Julia if you wish to run the atmospheric blocking version of the script for both the East and West blocks.

# References

[1] Jason Atnip, Gary Froyland, and Peter Koltai. An inflated dynamic Laplacian to track the emergence and disappearance of semi-material coherent sets, 2024

[2] H. Hersbach, B. Bell, P. Berrisford, G. Biavati, A. Horanyi, J. Munoz-Sabater, J. Nicolas, C. Peubey, R. Radu, I. Rozum, D. Schepers, A. Simmons, C. Soci, D. Dee, and J-N Thepaut. ERA5 hourly data on pressure levels from 1940 to present, Copernicus Climate Change Service (C3S) Climate Data Store (CDS). https://doi.org/10.24381/cds.bd0915c6, 2023.

# Acknowledgements

The authors of this repository wish to acknowledge Prof. Joerg Schumacher and Roshan John Samuel of Technische Universitat Ilmenau for their curation of the Rayleigh-Benard Convection velocity data, their analysis of the results and their helpful assistance and suggestions pertinent to this research.