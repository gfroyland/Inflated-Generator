# Readme File for the Switching Double Gyre

**Using a Stand Alone Julia REPL Window**

1. Download and save the "Inflated-Generator" repository to your system.
2. Open a new Julia REPL window and move to the "Inflated-Generator" directory.
3. Type "]", followed by the commands "activate ." and "instantiate" to set up the Julia environment for this repository.
4. Either:
    1. Run the double gyre script with: include("double gyre/generator_script_DG.jl")
    2. Run the atmospheric blocking script with: include("atmospheric blocking/generator_script_east_block.jl") for the Eastern block.  Substitute “east” with “west” to run the script for the Western block.
5. Data, images, and movies will be stored in either the double gyre or atmospheric blocking directory, depending on which script has been run.

**Using VS Code**

1. Download and save the "Inflated-Generator" repository to your system
2. Open VS Code, and open the "Inflated-Generator" folder in your workspace
3. Start a new Julia REPL in VS Code. Click on "Julia env" at the bottom of your VS Code window, select "(pick a folder)" from the drop down menu appearing at the top of the window, and find and select the Inflated-Generator folder in the system dialog
4. Type "]" followed by the commands "activate ." and "instantiate" to complete set up of the Julia environment for this repository
5. In the VS Code explorer sidebar, navigate to either the "atmospheric blocking" or "double gyre" folder, depending on which generator script you wish to run. 
6. In the VS Code explorer sidebar, left-click on a .jl file containing the word “script”. The code should open at the right. Click on the right-pointing triangle icon near the top right to run.
7. Data, images, and movies will be stored in either the double gyre or atmospheric blocking directory, depending on which script has been run.

# Double Gyre Sample Figure
Below is a three-dimensional visualisation of time slices of two almost-invariant sets detected for the switching Double Gyre system using the inflated generator method. The left axis of the figure displays the leading spatial eigenvector of the inflated generator, while the right axis of the figure displays two SEBA vectors generated from the first two eigenvectors of the inflated generator. This image can be reproduced using the MATLAB script "plot_slices_3D.m" and output data obtained from running "generator_script_DG.jl". All necessary code for the inflated generator calculations on the Double Gyre system can be found in the "double gyre" directory.

<img src = "https://github.com/gfroyland/Inflated-Generator/assets/168791783/9c79fbd8-ee85-4250-be97-03af57e6221e" width=600 >