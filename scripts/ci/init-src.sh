#!/bin/bash

# This script creates the src directory if it does not exist.  
# The directory is populated with a fresh clone of a git repository with 
# several remotes from defrost & anatoscope.  

# Exit on error
set -o errexit


### Checks
usage() {
    echo "Usage: init-src.sh <src-dir> <build-commit> <src-commit>"
}

if [[ "$#" = 3 ]]; then
    src_dir="$1"
    build_commit="$2"
    src_commit="$3"
else
    usage; exit 1
fi

echo "Working directory: "$src_dir 

doClone=0
if [ -d $src_dir ]; then
	echo "Check this is a git repository pointing to sofa framework."
	if git status > /dev/null ; then 
		echo "IN"
		doClone=0
	else
		echo "OUT"
		doClone=1
	fi
else
	doClone=1
fi

### Clone the git repository 
if (( doClone == 1 )); then
    echo "Creating a new repository at: "$src_dir
    git clone https://github.com/sofa-framework/sofa.git $src_dir
    
    if [[ "$(uname)" != "Darwin" && "$(uname)" != "Linux" ]]; then
	    echo "Copying dependency pack in the source tree."
	    /bin/cp -rf /c/dependencies/* $branch    
    fi
fi

echo "Fetch and checkout..."
if ! git config remote.anatoscope.url > /dev/null; then
    git remote add anatoscope https://github.com/anatoscope/sofa.git
fi
git fetch --all
git fetch origin +refs/pull/*:refs/remotes/origin/pr/*
git checkout --force $build_commit

