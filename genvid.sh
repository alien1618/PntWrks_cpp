#!/bin/bash

echo "Simulation video generation script started..."
i=1;
for var in "$@" 
do
    ffmpeg -i out/"$var"_jpg/"$var"%d.jpg -vcodec mpeg4 out/"$var".avi
    i=$((i + 1));
done
echo "Simulation video generation script complete..."
