# Inflated-Generator
This repository features Julia and MATLAB code used for the identification of quasi-stationary 
families of almost-invariant sets within velocity systems through the inflated generator. Code is available
for applying the inflated generator technique to two dynamic velocity systems:

1. The Switching Double Gyre, an idealised flow system as defined in Atnip, Froyland and Koltai (2024); and
2. ECMWF atmospheric data, observational wind velocity data used to identify two atmospheric blocking events occurring during the European Summer of 2003. The two events in question are a "West block", centred around France and taking place between August 2-12 2003, and an "East block", centred around Western Russia and taking place between July 29-August 3 2003.

# Downloading and Setting Up The Repository

Before using the code featured in this repository, you'll need to set up the required Julia environment on your system. The packages required to set up this environment feature within the Project.toml file, while the versions of these packages used and their subdependencies feature within Manifest.toml.

**Using a Julia REPL Stand Alone Window**

1. Download and save the "Inflated-Generator" repository to your system
2. Open a new Julia window and move to the "Inflated-Generator" directory
3. Type "]", followed by the commands "activate ." and "instantiate" to set up the Julia environment
for this repository.

**Using VS Code**

1. Download and save the "Inflated-Generator" repository to your system
2. Open VS Code, and open the "Inflated-Generator" folder in your workspace
3. Click on "Julia env" at the bottom of your VS Code window, select "(pick a folder)" from the drop down menu
appearing at the top of the window, and find and select the Inflated-Generator folder in the system dialog
4. Start a new REPL in VS Code, then type "]" followed by the commands "activate ." and "instantiate" to complete
set up of the Julia environment for this repository.

# Double Gyre Sample Figure
Here is a three-dimensional visualisation of time slices of quasi-stationary families of almost-invariant sets for the switching Double Gyre system. This image can be replicated using the MATLAB script "plot_slices_3D.m" using output data obtained from running "generator_script_DG.jl". All necessary code for the Double Gyre inflated generator calculations can be found in the "double gyre" directory.
<img src = "https://github.com/gfroyland/Inflated-Generator/assets/168791783/62afbbad-f46a-4e7f-baf4-9c85405ea945" width=600 >

# Atmospheric Blocking Sample Figure
Here is an image featuring 24-hour time slices of the 10th eigenvector of the inflated generator applied to the ECMWF velocity data Eastern block spatial domain. This image will be produced upon running "generator_script_east_block.jl". All necessary code and data for the Atmospheric Blocking inflated generator calculations can be found in the "atmospheric blocking" directory.
<img src = "https://github.com/gfroyland/Inflated-Generator/blob/aleks/Sample_Block_Eigvec_Fig.svg?raw=true" width=600 >

# References
