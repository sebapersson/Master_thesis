#!/bin/bash

# -------------------------------------------------------------------------
# Start of functions 
# -------------------------------------------------------------------------
# Function that will check if a directory exists. If the directory doesn't
# exist the function will create it. Note that input order matters 
# Args:
#     $1 Directory name 
check_if_dir_exists ()
{
    # Create directory if it doesn't exist
    if [ ! -d $1 ]; then
	mkdir $1
    fi
}

# Function that will create the mesh for a file
# Args:
#    $1, the file name
create_mesh ()
{
    ~/gmsh/gmsh-4.4.1-Linux64/bin/gmsh -2 $1 > /dev/null 
}

# -------------------------------------------------------------------------
# End of functions 
# -------------------------------------------------------------------------

# Check if the RD-simulations are to be run
echo "Should the RD-simulations be run (y/n)"
echo "A cluster is recommended for this"
read run_RD
if [ $run_RD == "y" ]; then
    echo "Will run RD-simulations"
    run_RD="true"
elif [ $run_RD == "n" ]; then
    echo "Will not run RD-simulations"
    run_RD="false"
else
    echo "Wrong input, program will exit"
    exit 1 
fi

# Ask for gmsh
echo "Have you provided the path to the gmsh binary on line 23 (y/n)"
read has_gmsh
if [ $has_gmsh == "n" ]; then
    echo "Add the gmsh binary, program will exit"
    exit 1
fi

# Ask for anaconda
echo "Have you installed the Anaconda environment (y/n)"
read has_environment
if [ $has_environment == "n" ]; then
    echo "Install anaconda, program will exit"
    exit 1
fi

# Ask if the user has R 
echo "Is R installed on the computer with the tidyverse library (y/n)"
read has_R
if [ $has_R == "y" ]; then
    echo "Will run R-part"
    run_R="true"
elif [ $has_R == "n" ]; then
    echo "Will not run R-part"
    echo "Note, a cluster is not required for this part."
    echo "So you can run it on a normal computer"
    run_R="false"
else
    echo "Wrong input, program will exit"
    exit 1 
fi

echo ""

# Check that the script is run from Shell directory
currentDir=${PWD##*/}

if [ ! $currentDir == "Shell" ]; then
    >&2 echo "The script most be run from Code/Shell directory"
    exit 1
fi

# Using the correct anaconda-environment
# Note that Anaconda most be activated to work in a sub-shell 
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate master_thesis

# Do the meshing
cd ../Gmsh/Rectangles/
echo "Creating the different meshes for rectangles"
echo ""
create_mesh Zero_holes.geo
create_mesh Five_holes.geo
create_mesh Seven_holes.geo
create_mesh Twenty_holes.geo

cd ../Circles
echo "Creating the different meshes for circles"
echo ""
create_mesh Zero_holes.geo
create_mesh Five_holes.geo
create_mesh Seven_holes.geo
create_mesh Twenty_holes.geo

# Move to local root
cd ../../..

# Remove files that might be present from earlier run
echo "Removing files from earlier runs"
echo ""
cd Intermediate/
rm -f -r Gierer_files
rm -f -r Schankenberg_files
rm -f -r Circles_mesh
rm -f -r Rectangles_mesh
cd ../Result
rm -f -r Gierer
rm -f -r Schankenberg
cd ..

# Run the PDE:s
cd Code/Python/RD_models
if [ $run_RD == "true" ]; then
   echo "Solving the PDE:s, this will take a while"
   echo "a while > 60 hours"
   echo ""
   ./Run_RD_simulations.py
fi

echo "Solving the illustration case"
./Simulate_one_dim_illustration.py

# Plot the Turing-space
cd ../Parameters
./Gierer_parameters.py
./Schnakenberg_parameters.py

# Move back to local root
cd ../../../

# Process the PDE-result
cd Code/R
if [ $run_R == "true" ]; then
    echo "Hello"
    Rscript ./Explore_data_sets.R 2> /dev/null
fi 
