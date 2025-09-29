# Readme File for Atmospheric Blocking

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

# Atmospheric Blocking Sample Figure
Below are several time slices of the SEBA vector highlighting the East block, spaced 24 hours apart.  One sees the block (highlighted in red) begin to appear on 29 July, strengthen from 30 July through to 2 July, and fade on 3 July. This image can be reproduced upon running "generator_script_east_block.jl". All necessary code and data for the inflated generator calculations on the ECMWF atmospheric velocity system can be found in the "atmospheric blocking" directory.

<img src = "https://github.com/gfroyland/Inflated-Generator/assets/168791783/95bcbd8b-103d-45cf-bbff-024be94c851e" width=600 >