#!/bin/bash

# This script computes the amound of remaining space. If it is below a
# given thresold it will remove the build directories from the oldest to 
# the newest 

# Exit on error
set -o errexit


### Checks

usage() {
    echo "Usage: clean-old-dir.sh <builds-dir> <current-build-dir> <needed size GB> <total size in GB>" 
}

toGB=$((1024*1024*1024))
toGB=$((1024))

if [[ "$#" = 4 ]]; then
    builds_dir="$1"
    current_dir="$2"
    neededspace=$(($toGB*$3))
    allocatedspace=$(($toGB*$4))
else
    usage; exit 1
fi

declare -i spaceused=`du -s $builds_dir | cut -d$'\t' -f1`
declare -i alreadyused=0

if ((spaceused > allocatedspace)); then
	spaceused=$allocatedspace
fi

if [[ -d "$current_dir" ]]; then
   declare -i alreadyused=`du -s $current_dir | cut -d$'\t' -f1`
   neededspace=$(($neededspace-$alreadyused)) 
   if ((neededspace < 0)); then
	neededspace=0
   fi
fi

echo "Building in: '"$builds_dir"'"
echo "- allocated space: "$allocatedspace" Gb"
echo "- currentsize    : "$spaceused" Gb"

availablespace=$(($allocatedspace-$spaceused)) 
if ((availablespace < 0)); then
	availablespace=0
fi
echo "- available space: "$availablespace" Gb"
echo "- needed space   : "$neededspace" Gb"
    		
freeed=0
for file in `ls -cr $builds_dir` ; do
    dir=$builds_dir/$file

    if [[ -d "$dir" ]]; then
    	if [[ "$dir" != "$current_dir" ]]; then
     		realspaceused=`du -s $builds_dir | cut -d$'\t' -f1`
     		thisdirspace=`du -s $dir | cut -d$'\t' -f1`
    		availablespace=$(($allocatedspace-$realspaceused)) 
   		if ((availablespace <= neededspace)); then
  		  	echo "Remove directory: $dir of "$realspaceused" Gb"
    			rm -rf $dir
		fi
	fi
    fi
done

### Real computation of the amount of free space after the cleaning.
realspaceused=`du -s $builds_dir | cut -d$'\t' -f1`
availablespace=$(($allocatedspace-$realspaceused))     		

echo "- free space     : "$availablespace" Gb"
if ((availablespace < neededspace)); then
	echo "Unable to free enough space. Please report the problem to ci@sofa-framework.org"
fi



