# Readme File for Atmospheric Blocking

This folder contains Julia code designed for the numerical implementation of the inflated generator method to identify two atmospheric blocking events which occurred during the European Summer of 2003.

This folder contains the following files/folders:

1. "Project.toml" and "Manifest.toml" files which detail the packages and version/subdependency information for these packages which are necessary for the replication of the Julia environment generated for the successful execution of this code.
2. "generator_functions.jl", which contains Julia functions used to (among other things) construct the inflated generator, classify the inflated generator eigenvalues/eigenvectors produced (as spatial or temporal), generate sample Figures and movies for the results and save relevant output data once the code has finished running.
3. "generator_script_east_block.jl", used to execute the full inflated generator method on our atmospheric velocity data by way of constructing the inflated generator, calculating the leading 10 eigenvalues and eigenvectors of this inflated generator (i.e. the 10 eigenvalues with largest real part and their corresponding eigenvectors); and generating a SEBA basis for the collection of real-valued spatial eigenvectors computed for the inflated generator. In this script, we execute the inflated generator method to identify the "East Block", which was centred around Western Russia (further East in Europe) and took place between 29 July-3 August 2003.
4. "generator_script_west_block.jl", a similar script to the one described above, however in this case we are trying to identify the "West Block", centred around France (further West in Europe) and taking place between 2-12 August 2003.
5. "SEBA.jl", which contains code used to run the SEBA algorithm. This algorithm is used to generate a sparse eigenbasis approximation from the spatial eigenvectors calculated for the inflated generator, in an effort to gain a clearer visualisation of the quasi-stationary families of almost-invariant sets present within our atmospheric flow system.
6. "data", a folder containing the ECMWF wind velocity data [1] files required to run the inflated generator method for this system.

# Downloading the Repository and Running the Code

Here are the instructions to follow for the successful execution of this code on your system: 

**Using a Stand Alone Julia REPL Window**

1. Download and save the "atmospheric blocking" repository to your system.
2. Open a new Julia REPL window and move to the "atmospheric blocking" directory.
3. Type "]", followed by the commands "activate ." and "instantiate" to set up the Julia environment for this repository.
4. Run the atmospheric blocking script with: include("generator_script_east_block.jl") for the Eastern block. Substitute “east” with “west” to run the script for the Western block.
5. Data, images, and movies will be stored within the atmospheric blocking directory.

**Using VS Code**

1. Download and save the "atmospheric blocking" repository to your system
2. Open VS Code, and open the "atmospheric blocking" folder in your workspace
3. Start a new Julia REPL in VS Code. Click on "Julia env" at the bottom of your VS Code window, select "(pick a folder)" from the drop down menu appearing at the top of the window, and find and select the "atmospheric blocking" folder in the system dialog
4. Type "]" followed by the commands "activate ." and "instantiate" to complete set up of the Julia environment for this repository
5. In the VS Code explorer sidebar, left-click on either one of the .jl files containing the word “script” (depending on which of the two blocking events you are most interested in). The code should open at the right. Click on the right-pointing triangle icon near the top right to run.
6. Data, images, and movies will be stored within the atmospheric blocking directory.

# Details on Code Output

Once either one of the inflated generator script files has been executed, a collection of images, movies and data files will be saved to the "atmospheric blocking" directory. A brief rundown of each of the files saved is as follows:

1. "InfGen_Results_EuroBlock_East.jld2" and "InfGen_Results_EuroBlock_East.h5" (or "InfGen_Results_EuroBlock_West.jld2" and "InfGen_Results_EuroBlock_West.h5"), a JLD2 file and an HDF5 file each containing the inflated generator eigenvalues, and data pertaining to both the eigenvectors of the inflated generator and the basis of SEBA vectors computed from spatial eigenvectors of the inflated generator for either the East or West block. Also included is the grid data for the spatial domain and the date range taken for the inflated generator calculations (each of which will be different between the two blocking events).
2. "Inflated Generator Eigenvalue Spectrum for the East Block.png" (or "Inflated Generator Eigenvalue Spectrum for the West Block.png"), a plot of the eigenvalue spectrum for the inflated generator, showing the leading 10 eigenvalues distinguished by their type (spatial or temporal).
3. "The East Block illustrated through SEBA vector $index_to_plot.png" (or "The West Block illustrated through SEBA vector $index_to_plot.png"), a Figure showing time slices of the SEBA vector which best illustrates the East (or West) block over its entire lifespan.
4. "Movie of the East Block illustrated through SEBA vector $index_to_plot.mp4" (or "Movie of the West Block illustrated through SEBA vector $index_to_plot.mp4"), a movie showing the temporal evolution of the SEBA vector which best illustrates the East (or West) block over its entire lifespan.

# Atmospheric Blocking Sample Figure

Below are several time slices of the SEBA vector highlighting the East block, spaced 24 hours apart.  One sees the block (highlighted in red) begin to appear on 29 July, strengthen from 30 July through to 2 July, and fade on 3 July. This image can be reproduced upon running "generator_script_east_block.jl". All necessary code and data for the inflated generator calculations on the ECMWF atmospheric velocity system can be found in the "atmospheric blocking" directory.

<img src = "https://github.com/gfroyland/Inflated-Generator/assets/168791783/95bcbd8b-103d-45cf-bbff-024be94c851e" width=600 >

# References

[1] H. Hersbach, B. Bell, P. Berrisford, G. Biavati, A. Horanyi, J. Munoz-Sabater, J. Nicolas, C. Peubey, R. Radu, I. Rozum, D. Schepers, A. Simmons, C. Soci, D. Dee, and J-N Thepaut. ERA5 hourly data on pressure levels from 1940 to present, Copernicus Climate Change Service (C3S) Climate Data Store (CDS). https://doi.org/10.24381/cds.bd0915c6, 2023.