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
create_mesh Twenty_holes.geo

cd ../Circles
echo "Creating the different meshes for circles"
echo ""
create_mesh Zero_holes.geo
create_mesh Five_holes.geo
create_mesh Twenty_holes.geo

# Move to local root
cd ../../..

# Remove files that might be present from earlier run
echo "Removing files from earlier runs"
echo ""
cd Intermediate/
rm -f -r Rectangles
rm -f -r Circles
cd ../Result
rm -f -r Rectangles_pwd_files
rm -f -r Circles_pwd_files
cd ..

# Run the PDe
echo "Solving the PDE:s, this will take a while"
echo ""
cd Code/Python/Schnakenberg
./Schnakenberg_holes_gmsh_domain.py

# Move back to local root
cd ../../../
