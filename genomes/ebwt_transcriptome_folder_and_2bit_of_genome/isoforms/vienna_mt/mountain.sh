#!/bin/bash  
for filename in $1/*.ps; do  
    ./mountain.pl < $filename > "$2/$(basename "$filename" _dp.ps).mt"  
done
