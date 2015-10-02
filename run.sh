#!/bin/bash

bootstraps=1000;
A=2; # Alphabet Size
n=500; # Sequence Length

replacement_percentage=100;
replacement_step=0.01;
replacements=$(echo "$replacement_percentage*$replacement_step*$n" | bc);
replacements=${replacements/\.*}

orig_file_name="origseq_n"$n".dat.gz"

echo "Replacements: " $replacements
echo "Sequence Length: " $n
echo "Bootstraps: " $bootstraps
echo "Alphabet Size: "$A
echo "Sequence File: " $orig_file_name

./main $orig_file_name $bootstraps $n $replacements $A
