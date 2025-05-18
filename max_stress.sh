#!/bin/bash

# find the highest absolute value of the stress on atoms
# you must enetr the path of the OUTCAR file, the nuber of atoms, and the file name of the temp file you want to use

path_to_file='chiralDNTT/OUTCAR'
n_atoms=192
filename=temp2



grep 'total drift' -B $(( $n_atoms+3 )) $path_to_file  | tail -n $(( $n_atoms+2 ))  > temp

awk 'NR < "'"$n_atoms"'" {print $4; $5; $6} ' temp > $filename


max=0

for line in $(cat $filename) ;  do

     temp_value="$line"

     abs_temp_value=${temp_value#-}

     if  (( $( echo "$abs_temp_value > $max" | bc -l  ) ))  ; then
          max=$abs_temp_value
     fi

done

echo "Maximum absolute value: $max"

grep -n "$max" temp
