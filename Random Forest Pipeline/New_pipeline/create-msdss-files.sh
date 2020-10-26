##############################################################################################################################################################################################################
# This shell script is designed to create all the R scripts, batch scripts, and swarm scripts that will be used to run the entire pipeline. 
# It relies on use of the Rswarm utility that is available on the Biowulf HPC. In all of the following, "username" will need to replaced by the username
# of individual whos Biowulf account is being used
##############################################################################################################################################################################################################

##############################################################################################################################################################################################################
# ACROSS ALL OF THESE COMMANDS WHERE A PATH IS INVOLVED, MAKE SURE TO DOUBLE CHECK THAT THEY ARE CORRECTLY SPECIFIED. 90% OF THE TIME
# THAT THERE WAS AN ERROR FOR ME WAS BECAUSE THESE DID NOT MATCH UP. IF FEELING AMBITIOUS, WOULD BE GOOD TO ADD THE PATH AS A VARIABLE
# AT THE TOP OF THE SCRIPT AND THEN USE THROUGHOUT FOR CONSISTENCY
##############################################################################################################################################################################################################

##############################################################################################################################################################################################################
# Set the working directory to be inside the folder where the template R scripts are and where the results will be saved.
##############################################################################################################################################################################################################
cd /data/(username)/pipeline_for_meeting/scripts/msdss/

##############################################################################################################################################################################################################
# Sometimes problems can happen if the next line is not included... 
##############################################################################################################################################################################################################
dos2unix seeds.txt

##############################################################################################################################################################################################################
# Looping through all of the iterations of the pipeline, creating the proper files for each iteration
##############################################################################################################################################################################################################
for j in {0..119}
do

##############################################################################################################################################################################################################
# The next lines simply copy the R script that will be used to create the forests, adding an "_iter(iteration number)" to the end of the file
# This will then be passed into the Rswarm utility in the next step. Code first checks if the file exists, deleting it if so. 
##############################################################################################################################################################################################################
if [ -e msdss_ratios_forest_iter${j}_.R ]
then
rm msdss_ratios_forest_iter${j}_.R
fi
cp msdss_ratios_forest.R msdss_ratios_forest_iter${j}_.R

##############################################################################################################################################################################################################
# The next lines runs the Rswarm utility to create 10 (reps argument) R scripts with individualized seeds and output file names. The swarm file that this utility creates is removed, as
# we will need to pass in additional command line argument via the batch R library/
##############################################################################################################################################################################################################
echo y | Rswarm -rfile=msdss_ratios_forest_iter${j}_.R --sfile=seeds.txt --path=/data/(username)/pipeline_for_meeting/scripts/msdss --start=0 --reps=10 --sims=1 --ext1=.rds &> /dev/null
rm msdss_ratios_forest_iter${j}_.sw

##############################################################################################################################################################################################################
# Check if file exists, delete if so
##############################################################################################################################################################################################################
if [ -e swarmR_iter${j}.sw ]
then
rm swarmR_iter${j}.sw 
fi

##############################################################################################################################################################################################################
# Create swarm file that allows command line arguments (using the --args statement at the beginning) to be passed into the R jobs. 
##############################################################################################################################################################################################################
touch swarmR_iter${j}.sw

##############################################################################################################################################################################################################
# NOTE: THE END OF THIS NEXT INDEX (10 IN THIS CASE) MUST MATCH THE --reps ARGUMENT IN THE RSWARM STATEMENT ABOVE
##############################################################################################################################################################################################################
for i in {1..10} 
do
echo "R --vanilla --args var_path \"./scripts/msdss/iter${j}_ratios.txt\" remove_outliers FALSE < /data/(username)/pipeline_for_meeting/scripts/msdss/msdss_ratios_forest_iter${j}_${i}.R" >> swarmR_iter${j}.sw
done

##############################################################################################################################################################################################################
# Check if file exists, delete if so
##############################################################################################################################################################################################################
if [ -e results_batch_iter${j} ]
then
rm results_batch_iter${j} 
fi

##############################################################################################################################################################################################################
# Create an index for the next iteration (which the next set of ratios/variables will use)
##############################################################################################################################################################################################################
i=$((j+1))

touch results_batch_iter${j}

##############################################################################################################################################################################################################
# Create a batch file that will be used for calculating model results at each iteration, using command line arguments via the batch R package
##############################################################################################################################################################################################################
echo "#!/bin/bash
# This file is runR
# 
#SBATCH -J runR
date
module load R
R --vanilla --args remove_outliers FALSE path1 "./scripts/msdss/msdss_ratios_forest_iter${j}_" pull_iter "./scripts/msdss/iter${j}" push_iter "./scripts/msdss/iter${i}" < /data/(username)/pipeline_for_meeting/scripts/msdss/msdss_ratios_forest_results.R" >> results_batch_iter${j}

done