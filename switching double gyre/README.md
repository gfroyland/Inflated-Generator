# Readme File for the Switching Double Gyre

This folder contains Julia code designed for the numerical implementation of the inflated generator method to identify quasi-stationary families of almost-invariant sets present within a switching Double Gyre flow system, along with MATLAB code used to produce a sample Figure for the generated results.

This folder contains the following files:

1. "Project.toml" and "Manifest.toml" files which detail the packages and version/subdependency information for these packages which are necessary for the replication of the Julia environment generated for the successful execution of this code.
2. "generator_functions.jl", which contains Julia functions used to (among other things) construct the inflated generator, calculate its eigenbasis using the eigs function from Arpack.jl, classify the eigenvalues/eigenvectors produced (as spatial or temporal), generate sample Figures and movies for the results and save relevant output data once the code has finished running.
3. "generator_script.jl", used to execute the full inflated generator method on the switching Double Gyre system by way of constructing the inflated generator, calculating the leading 10 eigenvalues and eigenvectors of this inflated generator (i.e. the 300 eigenvalues with largest real part and their corresponding eigenvectors); and generating a SEBA basis for the collection of real-valued spatial eigenvectors computed for the inflated generator.
4. "SEBA.jl", which contains code used to run the SEBA algorithm. This algorithm is used to generate a sparse eigenbasis approximation from the real-valued spatial eigenvectors calculated for the RBC inflated generator, in an effort to gain a clearer visualisation of the quasi-stationary families of almost-invariant sets present within the 3D RBC flow system.
5. "plot_slices_3D.m", a MATLAB file used to produce the below sample Figure showing two quasi-stationary families of almost-invariant sets identified for this system.
6. "InfGen_Results_SwitchingDoubleGyre.h5", an HDF5 file containing output from the inflated generator method applied to this system which can be used immediately to prepare the below sample Figure.
7. "bluewhitered.m", an additional MATLAB file used to define the blue-to-red colourmap featured in the below sample Figure.

# Downloading the Repository and Running the Code

Here are the instructions to follow for the successful execution of this code on your system: 

**Using a Stand Alone Julia REPL Window**

1. Download and save the "switching double gyre" repository to your system.
2. Open a new Julia REPL window and move to the "switching double gyre" directory.
3. Type "]", followed by the commands "activate ." and "instantiate" to set up the Julia environment for this repository.
4. Run the double gyre script with: include("generator_script.jl").
5. Data, images, and movies will be stored within the double gyre directory.

**Using VS Code**

1. Download and save the "switching double gyre" repository to your system
2. Open VS Code, and open the "switching double gyre" folder in your workspace
3. Start a new Julia REPL in VS Code. Click on "Julia env" at the bottom of your VS Code window, select "(pick a folder)" from the drop down menu appearing at the top of the window, and find and select the "switching double gyre" folder in the system dialog
4. Type "]" followed by the commands "activate ." and "instantiate" to complete set up of the Julia environment for this repository
5. In the VS Code explorer sidebar, left-click on the "generator_script.jl" file. The code should open at the right. Click on the right-pointing triangle icon near the top right to run.
6. Data, images, and movies will be stored within the double gyre directory.

# Details on Code Output

Once the inflated generator script file has been executed, a collection of images, movies and data files will be saved to the "switching double gyre" directory. A brief rundown of each of the files saved is as follows:

1. "InfGen_Results_SwitchingDoubleGyre.jld2" and "InfGen_Results_SwitchingDoubleGyre.h5", a JLD2 file and an HDF5 file each containing the inflated generator eigenvalues, and data pertaining to both the eigenvectors of the inflated generator and the basis of SEBA vectors computed from spatial eigenvectors of the inflated generator. Also included is the grid data for the spatial domain and the time range taken for the inflated generator calculations.
2. "Inflated Generator Eigenvalue Spectrum for the Double Gyre.png", a plot of the eigenvalue spectrum for the inflated generator, showing the leading 10 eigenvalues distinguished by their type (spatial or temporal).
3. "Leading Spatial Eigenvector for the Double Gyre.png", a Figure showing time slices of the leading real-valued spatial eigenvector of the inflated generator for the switching Double Gyre system spaced 0.1 units of time apart.
4. "Movie of leading spatial inflated generator eigenvector for the Double Gyre.gif", a movie showing the temporal evolution of the leading real-valued spatial eigenvector of the inflated generator for the switching Double Gyre system.
5. "SEBA vector *N* for the Double Gyre.png", a Figure showing time slices of the *N*-th SEBA vector computed for the switching Double Gyre system spaced 0.1 units of time apart. *N* will range from 1 to 2, as it is known a priori that two quasi-stationary families of almost-invariant sets exist within this system.
6. "Movie of SEBA vector *N* for the Double Gyre.gif", a movie showing the temporal evolution of the *N*-th SEBA vector computed for the switching Double Gyre system. Again, *N* will range from 1 to 2, as it is known a priori that two quasi-stationary families of almost-invariant sets exist within this system.
7. "Maxima of SEBA vectors for the Double Gyre.png", a Figure showing time slices of the maxima of all SEBA vectors computed for the switching Double Gyre system spaced 0.1 units of time apart.
8. "Movie of SEBA vector maxima for the Double Gyre.gif", a movie showing the temporal evolution of the maxima of all SEBA vectors computed for the switching Double Gyre system.

# Double Gyre Sample Figure

Below is a three-dimensional visualisation of time slices of two almost-invariant sets detected for the switching Double Gyre system using the inflated generator method. The left axis of the figure displays the leading spatial eigenvector of the inflated generator, while the right axis of the figure displays two SEBA vectors generated from the first two eigenvectors of the inflated generator. This image can be reproduced using the MATLAB script "plot_slices_3D.m" and output data obtained from running "generator_script.jl".

<img src = "https://github.com/gfroyland/Inflated-Generator/assets/168791783/9c79fbd8-ee85-4250-be97-03af57e6221e" width=600 >

# Changing Parameters Within The Scripts

After running "generator_script.jl" in its current form, you can rerun this script (and hence the inflated generator method for the switching Double Gyre system) after making some changes to some of the method's most crucial computational parameters. To make these changes, open "generator_script.jl" in VS Code or using a Nano Text Editor (preferably the former so that you can more easily navigate through the script file), and make the following alterations to the script to satisfy your needs:

1. On Line 6 of "generator_script.jl", change "Δt" to set the desired time step for the inflated generator method.
2. On Line 11 of "generator_script.jl", change "ℓ" to set the desired side length for the boxes we subdivide the spatial domain for the switching Double Gyre system into (and hence, set the resolution level of the results).
3. On Line 23 of "generator_script.jl", choose a new (positive) value for the temporal diffusion parameter a. This parameter is currently set to 0.115, as this gives us a reasonably good match between the leading spatial and temporal eigenvalues of the inflated generator.
4. On Line 37 of "generator_script.jl", change "num_of_Λ" to set the number of eigenvalues/eigenvectors of the inflated generator you wish to compute using the Arpack eigs function. This is currently set to 10, which should be sufficient for this system, but increase this if desired.
5. On Line 48 of "generator_script.jl", change "time_slice_spacing" to some integer greater than or equal to 1, if you only wish to include every time_slice_spacing-th slice after t_0 = 0 in the plots/movies of the eigenvectors/SEBA vectors produced.