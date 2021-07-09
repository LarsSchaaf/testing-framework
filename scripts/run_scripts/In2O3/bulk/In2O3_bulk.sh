#!/bin/bash
#$ -pe smp 24 # number of cores requested
#$ -l h_rt=99:00:00  # time requested in HH:MM:SS format
#$ -S /bin/bash      # shell to run the job ingit s
#$ -N batch_vac           # name of job (will appear in output of qstat)
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q any
#$ -cwd              # execute job in directory from which it was submitted

hostname
echo 'Bulk In - Phonons Ia3' # Title
echo 'testing-framework'                       # Type eg. gap_fit
echo 'Phonons' # Descirption
date +"%d/%m/%Y %H:%M:%S"
echo '--------------------------------------------------------' # For visual apeal

export PYTHONPATH="/opt/womble/QUIP/2021_02_08/linux_x86_64_gfortran/:$PYTHONPATH" # "/opt/womble/QUIP/latest/linux_x86_64_gfortran/:$PYTHONPATH"
export PYTHONPATH=${PYTHONPATH}:/opt/womble/QUIP/2021_02_08/linux_x86_64_gfortran_openmp/quippy
export PATH="/opt/womble/QUIP/2021_02_08/bin:$PATH"

export OMP_NUM_THREADS=${NSLOTS}
# Ensures quippy knows it can only use Nslots. Otherwise I reserve 16 and then use 32. Ie from another users reserved spot.

eval "$(/home/lls34/miniconda3/bin/conda shell.bash hook)" # Works
conda activate phd-code

cd /home/lls34/GitHub/01_PhD/PhD_Code/submodules/testing-framework/testing-framework/example_models/In2O3/run_dir

# model=DFT_QE
# model=GAP_it1_50s
# model=GAP_itB_1_c
model=GAP_itB_1_a_2
model=GAP_itB_3_b_2 # GAP_itB_2_b_1, GAP_itB_1_c_2 , GAP_itB_1_b_2
model=GAP_itB_4_b_3
model=GAP_inF_1
model=GAP_itB_1_a
model=GAP_inA_50s
model=GAP_itB_3_b_4 
model=GAP_itA_1_a
model=GAP_itB_3_b_5
model=GAP_itB_1_b_4_virials2
model=GAP_inA_50s_virials2
model=GAP_itB_3_b_7_virials2
#model=DFT_QE









args=''

echo model: $model
echo args: $args

echo Analysing Model@ $model

# date +"%T" >> time.txt
# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/bulk_In2O3_Ia3
# date +"%T" >> time.txt
# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/bulk_In2O3_R3c
# date +"%T" >> time.txt
# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/bulk_In2O3_Pbca
# date +"%T" >> time.txt
# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/bulk_In2O3_Pbcn
# date +"%T" >> time.txt
# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/phonons_In2O3_Ia3
# date +"%T" >> time.txt
python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/phonons_In2O3_R3c
# date +"%T" >> time.txt
# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/bulk_In 
# date +"%T" >> time.txt



# python ../../../scripts/run-model-test.py  -Nl In2O3 $model ../../../../tests/In2O3/phonons_In2O3_Ia3_relaxed
# 'GAP_inA_50s_no_virials' , 'GAP_itB_3_b_4_no_virials','GAP_itB_3_b_4_virials2','GAP_itB_3_b_4_virials3'  

# python ../../../scripts/run-model-test.py  -Nl In2O3 GAP_itB_1_b_4_virials2 $args ../../../../tests/In2O3/bulk_In2O3_Ia3
# python ../../../scripts/run-model-test.py  -Nl In2O3 GAP_itB_1_b_4_no_virials $args ../../../../tests/In2O3/bulk_In2O3_Ia3
# python ../../../scripts/run-model-test.py  -Nl In2O3 GAP_itB_3_b_4_virials2 $args ../../../../tests/In2O3/bulk_In2O3_Ia3
# python ../../../scripts/run-model-test.py  -Nl In2O3 GAP_itB_3_b_4_virials3 $args ../../../../tests/In2O3/bulk_In2O3_Ia3


#python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/vacancy_In2O3_singleunit_Ia3

# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/batch_vacancy_In2O3
# python ../../../scripts/run-model-test.py  -Nl In2O3 $model $args ../../../../tests/In2O3/vacancy_In2O3_Ia3_111



#vacancy_In2O3_Ia3_111
#vacancy_In2O3_R3c_222