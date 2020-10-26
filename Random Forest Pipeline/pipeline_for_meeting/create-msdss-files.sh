cd /data/barbourcr/pipeline_for_meeting/scripts/msdss/

# Sometimes this causes problems if it is not done, not sure why...
dos2unix seeds.txt

for j in {0..119}
do

if [ -e msdss_ratios_forest_iter${j}_.R ]
then
rm msdss_ratios_forest_iter${j}_.R
fi

cp msdss_ratios_forest.R msdss_ratios_forest_iter${j}_.R

# --PATH SOMETIMES NEEDS TO BE CHANGED
echo y | Rswarm -rfile=msdss_ratios_forest_iter${j}_.R --sfile=seeds.txt --path=/data/barbourcr/pipeline_for_meeting/scripts/msdss --start=0 --reps=10 --sims=1 --ext1=.rds &> /dev/null
rm msdss_ratios_forest_iter${j}_.sw

if [ -e swarmR_iter${j}.sw ]
then
rm swarmR_iter${j}.sw 
fi

touch swarmR_iter${j}.sw

for i in {1..10}
do
echo "R --vanilla --args var_path \"./scripts/msdss/iter${j}_ratios.txt\" remove_outliers FALSE < /data/barbourcr/pipeline_for_meeting/scripts/msdss/msdss_ratios_forest_iter${j}_${i}.R" >> swarmR_iter${j}.sw
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
R --vanilla --args remove_outliers FALSE path1 "./scripts/msdss/msdss_ratios_forest_iter${j}_" pull_iter "./scripts/msdss/iter${j}" push_iter "./scripts/msdss/iter${i}" < /data/barbourcr/pipeline_for_meeting/scripts/msdss/msdss_ratios_forest_results.R" >> results_batch_iter${j}

done