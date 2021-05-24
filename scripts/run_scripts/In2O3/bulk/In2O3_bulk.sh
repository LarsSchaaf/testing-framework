#!/bin/bash
#$ -pe smp 32 # number of cores requested
#$ -l h_rt=72:00:00  # time requested in HH:MM:SS format
#$ -S /bin/bash      # shell to run the job ingit s
#$ -N bulk_In           # name of job (will appear in output of qstat)
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q any
#$ -cwd              # execute job in directory from which it was submitted

hostname
echo 'Bulk In - Testing Framework' # Title
echo 'QE'                       # Type eg. gap_fit
echo 'Running Testing_Framework on womble with md' # Descirption
date +"%d/%m/%Y %H:%M:%S"
echo '--------------------------------------------------------' # For visual apeal

export PYTHONPATH="/opt/womble/QUIP/2021_02_08/linux_x86_64_gfortran/:$PYTHONPATH" # "/opt/womble/QUIP/latest/linux_x86_64_gfortran/:$PYTHONPATH"
export PYTHONPATH=${PYTHONPATH}:/opt/womble/QUIP/2021_02_08/linux_x86_64_gfortran_openmp/quippy
export PATH="/opt/womble/QUIP/2021_02_08/bin:$PATH"

export OMP_NUM_THREADS=1
# Ensures quippy knows it can only use Nslots. Otherwise I reserve 16 and then use 32. Ie from another users reserved spot.

eval "$(/home/lls34/miniconda3/bin/conda shell.bash hook)" # Works
conda activate phd-code

cd /home/lls34/GitHub/01_PhD/PhD_Code/submodules/testing-framework/testing-framework/example_models/In2O3/run_dir

# python ../../../scripts/run-model-test.py  -Nl In2O3 DFT_QE ../../../../tests/In2O3/bulk_In 
# python ../../../scripts/run-model-test.py  -Nl In2O3 DFT_QE ../../../../tests/In2O3/bulk_In2O3_Ia3
# python ../../../scripts/run-model-test.py  -Nl In2O3 DFT_QE ../../../../tests/In2O3/bulk_In2O3_Pbca
# python ../../../scripts/run-model-test.py  -Nl In2O3 DFT_QE ../../../../tests/In2O3/bulk_In2O3_Pbcn
# python ../../../scripts/run-model-test.py  -Nl In2O3 DFT_QE ../../../../tests/In2O3/bulk_In2O3_R3c
# python ../../../scripts/run-model-test.py  -Nl In2O3 DFT_QE ../../../../tests/In2O3/phonons_In2O3_Ia3
python ../../../scripts/run-model-test.py  -Nl In2O3 DFT_QE ../../../../tests/In2O3/phonons_In2O3_Ia3_relaxed
# python ../../../scripts/run-model-test.py  -Nl In2O3 GAP_it1_50s ../../../../tests/In2O3/phonons_In2O3_Ia3
#



