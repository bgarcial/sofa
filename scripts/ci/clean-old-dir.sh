#!/bin/bash

# This script computes the amound of remaining space. If it is below a
# given thresold it will remove the build directories from the oldest to 
# the newest 

# Exit on error
set -o errexit


### Checks

usage() {
    echo "Usage: clean-old-dir.sh <builds-dir> <needed size GB> <total size in GB>" 
}

#toGB=$((1024*1024*1024))
toGB=1

if [[ "$#" = 3 ]]; then
    builds_dir="$1"
    neededspace=$(($toGB*$2))
    allocatedspace=$(($toGB*$3))
else
    usage; exit 1
fi

declare -i spaceused=`du -s $builds_dir | cut -d$'\t' -f1`

echo "Building in: '"$builds_dir"'"
echo "- allocated space: "$allocatedspace
echo "- needed space   : "$neededspace
echo "- currentsize    : "$spaceused

freeed=0
for file in `ls -cr $builds_dir` ; do
    dir=$builds_dir/$file

    if [[ -d "$dir" ]]; then
    	spaceused=`du -s $builds_dir | cut -d$'\t' -f1`
    	availablespace=$(($allocatedspace-$spaceused)) 
    	echo "   available space: "$availablespace
    	if ((availablespace < neededspace)); then
  	  	echo "Remove first directory: $dir"
    		rm -rf $dir
    		#declare -i sf=`du -s $dir | cut -d$'\t' -f1`
		#echo "SF: $freeed + [$sf]"
		#freeed=$(($freeed + $sf))	    
	fi
    fi
done
if ((availablespace < neededspace)); then
	echo "Unable to free enough space. Please report the problem to ci@sofa-framework.org"
fi


