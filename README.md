# Inflated-Generator
This repository features Julia and MATLAB code used for the numerical implementation of the inflated generator, a time-extended version of the transfer operator, to two dynamic velocity systems in an effort to detect quasi-stationary families of almost-invariant sets within each of these systems. The two systems in question are:

1. The Switching Double Gyre, an idealised flow system involving two counter-rotating gyres of unequal size within a rectangular flow domain as defined in [1]; and
2. ECMWF wind velocity data [2] used to identify two atmospheric blocking events occurring during the European Summer of 2003. The two events in question are a "West block", centred around France (further West in Europe) and taking place between August 2-12 2003; and an "East block", centred around Western Russia (further East in Europe) and taking place between July 29-August 3 2003.

Included within this repository are function files in Julia used for the numerical construction of the inflated generator, and separate script files used to implement the method on each dynamic system of interest (in the atmospheric case, the West and East blocks are split into separate scripts as two separate spatial domains are used to identify and track each of these blocking events).

# Downloading and Setting Up The Repository

Before using the code featured in this repository on your system, you'll need to set up the required Julia environment featuring the packages necessary for the successful implementation of this code. The packages required to set up this environment feature within the Project.toml file, while information regarding the versions of these packages used and their subdependencies is contained within Manifest.toml.

**Using a Stand Alone Julia REPL Window**

1. Download and save the "Inflated-Generator" repository to your system
2. Open a new Julia REPL window and move to the "Inflated-Generator" directory
3. Type "]", followed by the commands "activate ." and "instantiate" to set up the Julia environment for this repository
4. Change directory to "atmospheric blocking" or "double gyre", depending on which generator script you'd prefer to run
5. Run the generator script of your choice using the include() function.

**Using VS Code**

1. Download and save the "Inflated-Generator" repository to your system
2. Open VS Code, and open the "Inflated-Generator" folder in your workspace
3. Click on "Julia env" at the bottom of your VS Code window, select "(pick a folder)" from the drop down menu appearing at the top of the window, and find and select the Inflated-Generator folder in the system dialog
4. Start a new Julia REPL in VS Code, then type "]" followed by the commands "activate ." and "instantiate" to complete set up of the Julia environment for this repository
5. Change directory to "atmospheric blocking" or "double gyre", depending on which generator script you'd prefer to run
6. Run the generator script of your choice using the include() function.

# Double Gyre Sample Figure
Here is a three-dimensional visualisation of quasi-stationary families of almost-invariant sets detected for the switching Double Gyre system. This image can be reproduced using the MATLAB script "plot_slices_3D.m" and output data obtained from running "generator_script_DG.jl". All necessary code for the inflated generator calculations on the Double Gyre system can be found in the "double gyre" directory.

<img src = "https://github.com/gfroyland/Inflated-Generator/assets/168791783/62afbbad-f46a-4e7f-baf4-9c85405ea945" width=600 >

# Atmospheric Blocking Sample Figure
Here is an image featuring time slices (with a step of 24 hours between each one) of the 10th eigenvector of the inflated generator for the ECMWF velocity data under the East block spatial domain. This image can be reproduced upon running "generator_script_east_block.jl". All necessary code and data for the inflated generator calculations on the ECMWF atmospheric velocity system can be found in the "atmospheric blocking" directory.

<img src = "https://github.com/gfroyland/Inflated-Generator/blob/aleks/Sample_Block_Eigvec_Fig.svg?raw=true" width=600 >

# References

[1] Jason Atnip, Gary Froyland, and Peter Koltai. An inflated dynamic Laplacian to track the emergence and disappearance of semi-material coherent sets, 2024

[2] H. Hersbach, B. Bell, P. Berrisford, G. Biavati, A. Horanyi, J. Munoz-Sabater, J. Nicolas, C. Peubey, R. Radu, I. Rozum, D. Schepers, A. Simmons, C. Soci, D. Dee, and J-N Thepaut. ERA5 hourly data on pressure levels from 1940 to present, Copernicus Climate Change Service (C3S) Climate Data Store (CDS). https://doi.org/10.24381/cds.bd0915c6, 2023.
