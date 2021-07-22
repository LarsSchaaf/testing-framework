#!/bin/bash
#$ -pe smp 16 # number of cores requested
#$ -l h_rt=72:00:00  # time requested in HH:MM:SS format
#$ -S /bin/bash      # shell to run the job ingit s
#$ -N relax_all            # name of job (will appear in output of qstat)
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q any
#$ -cwd              # execute job in directory from which it was submitted

hostname
echo 'Cutoff MD - Locality Testing' # Title
echo 'quippy'                       # Type eg. gap_fit
echo 'Done on Edvins 110 strucutre' # Descirption
date +"%d/%m/%Y %H:%M:%S"
echo '--------------------------------------------------------' # For visual apeal

export OMP_NUM_THREADS=1
export PYTHONPATH="/opt/womble/QUIP/2021_02_08/linux_x86_64_gfortran/:$PYTHONPATH"  # "/opt/womble/QUIP/latest/linux_x86_64_gfortran/:$PYTHONPATH"
export PYTHONPATH=${PYTHONPATH}:/opt/womble/QUIP/2021_02_08/linux_x86_64_gfortran_openmp/quippy
export PATH="/opt/womble/QUIP/2021_02_08/bin:$PATH"
# Ensures quippy knows it can only use Nslots. Otherwise I reserve 16 and then use 32. Ie from another users reserved spot.

eval "$(/home/lls34/miniconda3/bin/conda shell.bash hook)" # Works
conda activate phd-code

python relax.py $1 $2
