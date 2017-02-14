#!/bin/bash

# This script initialize the src directory.
# If it does not exist it is created and populated with a fresh clone of a git repository with 
# several remotes from defrost & anatoscope.  
#
# It then checkout the src-commit number or sha1. 
# A second build-commit is the number needed by github to reference the build.


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
	cd $src_dir
	echo "Checking that '$src_dir' this is a valid git repository pointing to sofa framework."
	if git status > /dev/null ; then 
		echo "...YES (we can continue)"
		doClone=0
	else
		echo "...NO (we must first do a fresh clone"
		doClone=1
	fi
	cd ..
else
	doClone=1
fi


### Clone the git repository 
if (( doClone == 1 )); then
    if [ -d $src_dir ]; then
        rm -rf $src_dir
    fi

    echo "Creating a new repository at: "$src_dir
    git clone https://github.com/sofa-framework/sofa.git $src_dir
    
    if [[ "$(uname)" != "Darwin" && "$(uname)" != "Linux" ]]; then
	    echo "Copying dependency pack in the source tree."
	    /bin/cp -rf /c/dependencies/* $branch    
    fi
fi

cd "$src_dir"

echo "The clone is ready...fetch and checkout..."
if ! git config remote.anatoscope.url > /dev/null; then
    git remote add anatoscope https://github.com/anatoscope/sofa.git
fi
git fetch --all
git fetch origin +refs/pull/*:refs/remotes/origin/pr/*
git checkout --force $build_commit

cd "../"
