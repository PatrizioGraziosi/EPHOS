#!/bin/bash

n_q=3
n_modes=300
n_atoms=100
filename='qpoints_amplitude_0.005_800eV.yaml'

for iq in `seq 1 $n_q ` ; do
        for im in `seq 1 $n_modes ` ; do
		nn=$(( 9+3*$im-1+($im-1)*(4*$n_atoms)+($iq-1)*(3+$n_modes*(3+4*$n_atoms)) ))
		sed -n "$nn"p $filename > temp
		awk '{printf("% 18.12f  \n",$2)}' < temp  >> freq.dat
        done
done
