#!/bin/bash -l
#SBATCH --time=03:14:59                 # Walltime in hh:mm:ss
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=4             # Number of MPI ranks per node
#SBATCH --cpus-per-task=1               # Number of OpenMP threads for each MPI process/rank
#SBATCH --gpus-per-task=1               # number of gpu per task
#SBATCH --gpus-per-node=4
#SBATCH --mem=470000                    # Per nodes memory request (MB)
#SBATCH --partition=gpu                 # partition
#SBATCH --account=p201002               # project account
#SBATCH --qos=default                   # SLURM qos

#SBATCH --job-name=chiTTF

#SBATCH --gres=gpu:4

#SBATCH --partition=gpu
#SBATCH --output=vasp_gpu_test.%j.out
#SBATCH --error=vasp_gpu_test.%j.err


module load env/release/2024.1
module use /project/home/p201002/install_lxp/modules/all
module load VASP/6.5.1-NVHPC-25.1-CUDA-12.6.0
export LD_LIBRARY_PATH="${EBROOTNVHPC}/Linux_x86_64/25.1/compilers/extras/qd/lib/:${EBROOTNVHPC}/Linux_x86_64/25.1/math_libs/11.8/lib64:${EBROOTNVHPC}/Linux_x86_64/25.1/cuda/11.8/targets/x86_64-linux/lib:${LD_LIBRARY_PATH}"


cp INCAR_mod INCAR

#srun vasp_ncl > scf.out
srun vasp_std > scf.out

cp OUTCAR OUTCAR_scf
cp INCAR_mod_nscf INCAR

#mpirun vasp_ncl > nscf.out
srun vasp_std > nscf.out

/home/users/u103342/c2x_2.42a/c2x --bxsf  EIGENVAL chiral_DNTT_single.bxsf



rm WAVECAR
rm CHG*
