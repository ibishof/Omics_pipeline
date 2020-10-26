############################################################
# Script to be run after all pipelines have finished, which
# will clear out all of the temporary files

# Inherint to root directory
rm *.e
rm *.o
rm *.out
rm msdss.sh

# Pipeline specific, additional statements like so can be added for each iteration
cd /data/(username)/pipeline_for_meeting/scripts/msdss/

rm *.sw
rm msdss_ratios_forest_iter*
rm results_batch_iter*