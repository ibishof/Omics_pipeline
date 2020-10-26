cd /home/wangq13/pipeline_for_meeting/scripts/diagnosis/

# Sometimes this causes problems if it is not done, not sure why...
dos2unix seeds.txt

for j in {0..119}
do

if [ -e diagnosis_ratios_forest_iter${j}_.R ]
then
rm diagnosis_ratios_forest_iter${j}_.R
fi

cp diagnosis_ratios_forest.R diagnosis_ratios_forest_iter${j}_.R

# --PATH SOMETIMES NEEDS TO BE CHANGED
echo y | Rswarm -rfile=diagnosis_ratios_forest_iter${j}_.R --sfile=seeds.txt --path=/home/wangq13/pipeline_for_meeting/scripts/diagnosis --start=0 --reps=10 --sims=1 --ext1=.rds &> /dev/null
rm diagnosis_ratios_forest_iter${j}_.sw

if [ -e swarmR_iter${j}.sw ]
then
rm swarmR_iter${j}.sw 
fi

touch swarmR_iter${j}.sw

for i in {1..10}
do
echo "R --vanilla --args var_path \"./scripts/diagnosis/iter${j}_ratios.txt\" remove_outliers FALSE < /home/wangq13/pipeline_for_meeting/scripts/diagnosis/diagnosis_ratios_forest_iter${j}_${i}.R" >> swarmR_iter${j}.sw
done

if [ -e results_batch_iter${j} ]
then
rm results_batch_iter${j} 
fi

i=$((j+1))

touch results_batch_iter${j}

echo "#!/bin/bash
# This file is runR
# 
#SBATCH -J runR
date
module load R
R --vanilla --args remove_outliers FALSE path1 "./scripts/diagnosis/diagnosis_ratios_forest_iter${j}_" pull_iter "./scripts/diagnosis/iter${j}" push_iter "./scripts/diagnosis/iter${i}" < /home/wangq13/pipeline_for_meeting/scripts/diagnosis/diagnosis_ratios_forest_results.R" >> results_batch_iter${j}

done