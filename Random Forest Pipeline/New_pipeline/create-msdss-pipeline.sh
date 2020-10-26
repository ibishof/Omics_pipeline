##############################################################################################################################################################################################################
# This shell script creates another shell script (titled "msdss.sh" in this example) that, when run, submits all swarm and batch jobs for the pipeline with the correct dependency structure.
# These files are created using the corresponding create-*-files.sh script  For the random forest models submitted via swarm, the cpus and memory can be allocated using the -t and -g arguments. 
# In all of the following, "username" will need to replaced by the username of individual whos Biowulf account is being used
##############################################################################################################################################################################################################

##############################################################################################################################################################################################################
# ACROSS ALL OF THESE COMMANDS WHERE A PATH IS INVOLVED, MAKE SURE TO DOUBLE CHECK THAT THEY ARE CORRECTLY SPECIFIED. 90% OF THE TIME
# THAT THERE WAS AN ERROR FOR ME WAS BECAUSE THESE DID NOT MATCH UP. IF FEELING AMBITIOUS, WOULD BE GOOD TO ADD THE PATH AS A VARIABLE
# AT THE TOP OF THE SCRIPT AND THEN USE THROUGHOUT FOR CONSISTENCY
##############################################################################################################################################################################################################

if [ -e msdss.sh ]
then
rm msdss.sh
fi

touch msdss.sh
echo "jid_01=\$(swarm -f /data/(username)/pipeline_for_meeting/scripts/msdss/swarmR_iter0.sw -t 16 -g 12 --partition=quick,norm --job-name msdss --module R)" >> msdss.sh
echo "jid_02=\$(sbatch --depend=afterok:\$jid_01 --mem=6g --partition=quick,norm --job-name=msdss /data/(username)/pipeline_for_meeting/scripts/msdss/results_batch_iter0)" >> msdss.sh

for j in {1..119}
do

i=$((j-1))

echo "jid_${j}1=\$(swarm -f /data/(username)/pipeline_for_meeting/scripts/msdss/swarmR_iter${j}.sw -t 16 -g 12 --partition=quick,norm --depend=afterok:\$jid_${i}2 --job-name msdss --module R)" >> msdss.sh
echo "jid_${j}2=\$(sbatch --depend=afterok:\$jid_${j}1 --mem=6g --partition=quick,norm --job-name=msdss /data/(username)/pipeline_for_meeting/scripts/msdss/results_batch_iter${j})" >> msdss.sh
done