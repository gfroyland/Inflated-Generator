# Readme File for Rayleigh-Benard Convection

This folder contains Julia code designed for the numerical implementation of the inflated generator method to identify quasi-stationary families of almost-invariant sets present within a three-dimensional Rayleigh-Benard Convection (RBC) flow system. 

There are two versions of the Julia code required to run this method on the RBC system, contained within the two subfolders (named "InfGen_RBC_CPU" and "InfGen_RBC_GPU") which each contain the following files:

1. "Project.toml" and "Manifest.toml" files which detail the packages and version/subdependency information for these packages which are necessary for the replication of the Julia environment generated for the successful execution of this code.
2. "generator_functions.jl", which contains Julia functions used to (among other things) construct the inflated generator, calculate its eigenbasis using the Arnoldi method, classify the eigenvalues/eigenvectors produced (as spatial, temporal or complex), generate sample Figures for the results and save relevant output data once the code has finished running.
3. "generator_script.jl", used to execute the full inflated generator method on the 3D RBC flow system by way of constructing the inflated generator, calculating the leading 300 eigenvalues and eigenvectors of this inflated generator (i.e. the 300 eigenvalues with largest real part and their corresponding eigenvectors); and generating a SEBA basis for the collection of real-valued spatial eigenvectors computed for the inflated generator.
4. "SEBA.jl", which contains code used to run the SEBA algorithm. This algorithm is used to generate a sparse eigenbasis approximation from the real-valued spatial eigenvectors calculated for the RBC inflated generator, in an effort to gain a clearer visualisation of the quasi-stationary families of almost-invariant sets present within the 3D RBC flow system.

The difference between the code files in each folder is that in the version of the code contained within the "InfGen_RBC_GPU" directory, computationally expensive mathematical processes, in particular the multiplication of large matrices as part of the Arnoldi method or the SEBA algorithm, are performed using arrays defined on GPU cores rather than a standard CPU processor in order to improve the efficiency of the inflated generator method. If you do not have access to a PC or a HPC cluster containing GPU cores, you can still execute the inflated generator method by downloading the code contained within the "InfGen_RBC_CPU" directory. Users should download only one of these subfolders depending on the capabilities of the system on which you are running the code.

When running the code featured in this repository, please ensure that you run it with a (as of October 2025) new version of Julia such as v1.10.9 or v1.11.6. Do not run this code with Julia versions v1.10.7 or v1.11.1, or the code featured within this repository will fail due to circular dependency issues present with the CUDA.jl package on these versions of Julia. For more information on these circular dependency issues, we refer you to this link: https://discourse.julialang.org/t/circular-dependency-warning/123388.

For the CPU version of this code, given a value of the parameter ℓ = 1/4 (the box side length used to construct the spatial generators for the overall inflated generator), on a laptop we expect this code to take roughly eight minutes. Of all of the processes which form this method, eigensolving the inflated generator using the Arnoldi method is the longest process, coming in at roughly six minutes.

# Downloading the RBC Velocity Data

Unfortunately, owing to their large file sizes, the RBC velocity data required to run this code cannot be downloaded from this repository. This velocity data, along with a subset of results data which can be used to replicate the Figures seen in 

Aleksandar Badza, Gary Froyland, Roshan J. Samuel and Joerg Schumacher. "Transient almost-invariant sets reveal convective heat transfer patterns in plane-layer Rayleigh-Benard convection", arXiv details TBA,

can be downloaded from the Zenodo repository available at this link: *Details to follow...*

# Downloading the Repository and Running the Code

Here are the instructions to follow for the successful execution of this code on your system: 

**Using a Stand Alone Julia REPL Window**

1. Download and save either the "InfGen_RBC_CPU" or "InfGen_RBC_GPU" subrepository to your system. If your PC has GPU capabilities or you are running this code on a HPC cluster, download and save the latter, otherwise if you do not have these GPU capabilities at hand, download and save the former.
2. Download the "RBC Velocity Data" folder from the aforementioned Zenodo repository and save it to the "InfGen_RBC_CPU" (or "InfGen_RBC_GPU") folder on your system.
3. Open a new Julia REPL window and move to the "InfGen_RBC_CPU" (or "InfGen_RBC_GPU") directory.
4. Set up the Julia environment for this repository by hitting the "]" key on your keyboard, followed by typing the commands "activate ." and "instantiate", and hitting the backspace key once these processes have been completed.
5. Run the RBC inflated generator method script with: include("generator_script.jl")
6. Data, images, movies and text files will be stored in the "InfGen_RBC_CPU" (or "InfGen_RBC_GPU") directory.

**Using VS Code**

1. Download and save either the "InfGen_RBC_CPU" or "InfGen_RBC_GPU" subrepository to your system. If your PC has GPU capabilities or you are running this code on a HPC cluster, download and save the latter, otherwise if you do not have these GPU capabilities at hand, download and save the former.
2. Download the "RBC Velocity Data" folder from the aforementioned Zenodo repository and save it to the "InfGen_RBC_CPU" (or "InfGen_RBC_GPU") folder on your system.
3. Open VS Code, and open the "InfGen_RBC_CPU" (or "InfGen_RBC_GPU") folder in your workspace.
4. Start a new Julia REPL in VS Code. Click on "Julia env" at the bottom of your VS Code window, select "(pick a folder)" from the drop down menu appearing at the top of the window, and find and select the "InfGen_RBC_CPU" (or "InfGen_RBC_GPU") folder in the system dialog.
5. Complete the set up of the Julia environment for this repository by hitting the "]" key on your keyboard, followed by typing the commands "activate ." and "instantiate", and hitting the backspace key once these processes have been completed.
6. In the VS Code explorer sidebar, left-click on the "generator_script.jl" file within that folder. The code should open at the right. Click on the right-pointing triangle icon near the top right to run the script file.
7. Data, images, movies and text files will be stored in the "InfGen_RBC_CPU" (or "InfGen_RBC_GPU") directory.

# Details on Code Output

Once an inflated generator script file has been executed, a text file along with a collection of images, movies and data files will be saved to either the "src/cpu_code" or "src/gpu_code" directories; depending on which version of the script has been run. A brief rundown of each of the files saved is as follows:

1. "InfGen Parameters.txt", a text file detailing key output parameters pertinent to the inflated generator method including, but not limited to: the spatial domain and time interval for the inflated generator calculations; the side length of each box that the spatial domain is subdivided into, the length of each time step along with the total number of time steps taken; values calculated for the spatial and temporal diffusion parameters ϵ and a; computation times for the construction of the spatial generators, calculation of the inflated generator eigenvalues/eigenvectors and generation of the SEBA vectors; indices of the spatial and temporal inflated generator eigenvectors; and the minima of the SEBA vectors (listed in order from closest to zero to furthest from zero; which is how the SEBA vectors are ranked and assigned their indices).
2. "RBC_InfGen_Spatial_Generator_Data.jld2", a JLD2 file containing the spatial generators computed for the RBC flow at each time step, along with the grid data for the spatial domain and the time range and time step taken for the inflated generator calculations.
3. "RBC_InfGen_Eigenbasis_Results.h5", an HDF5 file containing the inflated generator eigenvalues, vectors with indices used to identify the type of each eigenvalue (spatial, temporal or complex) and data pertaining to the real-valued spatial eigenvectors of the inflated generator (these are the only eigenvectors saved, in an effort to reduce the size of this data file and also because these eigenvectors are the only ones that are practical for use in finding quasi-stationary families of almost-invariant sets). Also included is the grid data for the spatial domain and the time range taken for the inflated generator calculations.
4. "RBC_InfGen_SEBA_Data.h5", an HDF5 file containing the SEBA vectors calculated from the real-valued spatial inflated generator eigenvectors, and interpolated versions of these vectors on a finer grid corresponding to half the box side length used in the original inflated generator calculations. Also saved is the spatial grid data for the original grid, the spatial grid data for the finer grid over which the SEBA data is interpolated, and the time range data.
5. "RBC Inflated Generator Spectrum.png", a plot of the eigenvalue spectrum for the inflated generator, showing the leading 300 eigenvalues distinguished by their type (spatial, temporal or complex).
6. "RBC SEBA Maxima xy Midplane.mp4", a movie showing the temporal evolution of the maxima of all SEBA vectors computed for the RBC flow system taken along the xy-plane with z = 0.5 (in other words, the xy "midplane").
7. "RBC SEBA Maxima xy Midplane.png", a gridded Figure showing time slices of the maxima of all SEBA vectors computed for the RBC flow system in the xy-plane with z = 0.5 spaced 3 T_f (units of time) apart.
8. "RBC SEBA Maxima Five Planes.png", a Figure showing the maxima of all SEBA vectors computed for the RBC flow system taken at the middlemost moment of our time interval restricted to five two-dimensional planes within our full three-dimensional domain. The five planes in question are: the xy-Bottom (or "xy-Floor", a horizontal xy-plane taken within the lower thermal boundary layer of our fluid cell), the xy-midplane (the horizontal xy-plane with z = 0.5), the xy-Top (or "xy-Ceiling", a horizontal xy-plane taken within the upper thermal boundary layer of our fluid cell), the xz-midplane (the vertical xz-plane with y = 0) and the yz-midplane (the vertical yz-plane with x = 0).

# Sample Figures

The Figure below shows time slices of the maxima of 30 SEBA vectors generated for the RBC flow system taken along the xy-plane with z = 0.5 (the xy-midplane) and spaced 3 T_f (units of time) apart. Quasi-stationary families of almost-invariant sets are identifiable through sub-regions of the domain coloured in deep red which roughly maintain their shape and do not fade in colour for at least 30-40 units of time. This Figure can be reproduced by running the GPU version of "generator_script.jl" in its current form. The CPU version of this script will also produce this Figure, but the code will take much longer to execute if we use the same resolution level (governed by ℓ) as that taken in the GPU version of the code.

<!--<img width="5000" height="4000" alt="Image" src="https://github.com/user-attachments/assets/adcbd6ba-6513-40f2-81f8-105e5cfc250e" />-->

<!--<img width="5000" height="4000" alt="Image" src="https://github.com/user-attachments/assets/98a10f6a-69b8-4c95-8f38-5193f615746f" />-->

<img width="5000" height="4000" alt="Image" src="https://github.com/user-attachments/assets/c919f73f-5213-4cbb-bb39-c6072b622899" />

The Figure below shows the maxima of 30 SEBA vectors generated for the RBC flow system taken along five two-dimensional restrictions of our three-dimensional spatial domain at the middlemost point of our time interval, 2052 T_f. This Figure can also be reproduced by running the GPU version of "generator_script.jl" in its current form. Again, the CPU version of this script will also produce this Figure, however it will take much longer for the code to execute at this fine resolution level.

<!--<img width="1000" height="800" alt="Image" src="https://github.com/user-attachments/assets/89d58721-f307-4d19-835c-e70f518b652b" />-->

<!--<img width="1000" height="800" alt="Image" src="https://github.com/user-attachments/assets/b99cb872-8ea5-4ede-a81c-d929a20f1bc1" />-->

<img width="4999" height="4999" alt="Image" src="https://github.com/user-attachments/assets/a6139d1c-cdad-45b8-a8c1-cd1b462f331e" />

# Changing Parameters Within The Scripts

After running either version of "generator_script.jl" in their current forms, you can rerun this script (and hence the inflated generator method for the RBC flow system) after making some changes to some of the method's most crucial computational parameters. To make these changes, open the preferred version of "generator_script.jl" (depending on your system's capabilities) in VS Code or using a Nano Text Editor (preferably the former so that you can more easily navigate through the script file), and make the following alterations to the script to satisfy your needs:

1. On Line 17 of "generator_script.jl", change "v_orientation" to 1 to run the inflated generator method for +v, or positive RBC velocity data (this is currently set to -1, to run the method on -v or negative RBC velocity data).
2. On Line 30 of "generator_script.jl", change "ℓ" to set the desired side length for the cubes we subdivide the 3D spatial domain for the RBC flow system into (and hence, set the resolution level of the results). Ideally, ℓ should take a value of 1/(2^(n)), where n is a positive integer greater than or equal to 1, and the smaller the value of ℓ is, the better the resolution of the results will be. However, it must be stressed that when using the CPU version of the script the runtime of the inflated generator method will start to increase significantly as the value of ℓ decreases. The difference in runtimes between the two versions of the script will not differ too significantly for the currently set value of ℓ (which is 1/4), but as ℓ gets smaller the size of the inflated generator matrix will get larger, and the large matrix multiplications that will not be performed on GPU cores when we reach the Arnoldi method and SEBA algorithm will cause the runtime for the CPU version of the script to increase at a much faster rate than the runtime for the GPU version of the script.
3. On Line 44 of "generator_script.jl", choose a new value for the temporal diffusion parameter a. Choose any positive value of a that you wish, though ideally it should be larger (but not too much larger) than the heuristic for a previously computed when the script was first run. Consult the "InfGen Parameters.txt" text file in the "src/cpu_code" or "src/gpu_code" folder to uncover what that heuristic of a was computed to be.
4. On Line 90 of "generator_script.jl", change "num_of_Λ" to set the number of eigenvalues/eigenvectors of the inflated generator you wish to compute using the Arnoldi method (this is currently set to 300, which should be sufficient, but increase this if desired, bearing in mind that the runtime of the code will increase for larger values of this parameter regardless of which version of the script is used).
5. On Line 91 of "generator_script.jl", change "tol" to set a different eigenvalue/eigenvector residual tolerance for the Arnoldi method if you wish. tol is currently set to √eps() or roughly 1.49e-08, which is the default tolerance for the Arnoldi method.

