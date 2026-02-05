#!/bin/bash
#SBATCH --time=01:59:00
#SBATCH --job-name=DPsTIPS

n_q=3
n_modes=3
n_atoms=100
total_modes=300
modes_to_skip=0
fileName_bxsf="'TIPS.bxsf'"
material_name="'TIPS'"
fileName="modes_info.m"

printf " n_q = $n_q; \n n_modes = $n_modes; \n n_atoms = $n_atoms; \n total_modes = $total_modes; \n modes_to_skip = $modes_to_skip; \n fileName = $fileName_bxsf; \n material_name = $material_name; \n" > $fileName


module load profile/eng
# module load intel/oneapi-2021--binary
# module load intelmpi/oneapi-2021--binary
module load matlab/r2022b



sed -n "3,5p" ../MPOSCAR-orig > A_matrix.txt

line_POSCAR=$(( $n_atoms+8 ))
sed -n "9,${line_POSCAR}p" ../MPOSCAR-orig > coord_initial.txt


# the matlab script that generates the line for the setting....conf file runs, for each q-point, all the modes
counter=0
for iq in `seq 1 $n_q ` ; do
	for im in `seq 1 $n_modes ` ; do

        	counter=$(( $counter+1 ))
		if (( $counter <= 9 )); then
			mod_dir="mod-00$counter"
		elif (( $counter <= 99 )); then
                        mod_dir="mod-0$counter"
		else
			mod_dir="mod-$counter"
                fi
		cd "$mod_dir"
			cp ../bxsf_to_ELECTRA_shifting_option.m ./
			cp ../shifting_bands.m ./
                	sed -n "9,${line_POSCAR}p" POSCAR > coord.txt
                cd ../
	done
done
cd mod-all
     cp ../bxsf_to_ELECTRA_shifting_option.m ./
     cp ../shifting_bands.m ./
     sed -n "9,${line_POSCAR}p" POSCAR > coord.txt
cd ../


cp mod-orig/*bxsf ./

matlab -nosplash <transfer_integrals_mod_extract_DefPot_script_L_v5.m> DefPot.log
matlab -nosplash <average_DPs.m>DP.out
