if [ -e msdss.sh ]
then
rm msdss.sh
fi

touch msdss.sh
echo "jid_01=\$(swarm -f /data/barbourcr/pipeline_for_meeting/scripts/msdss/swarmR_iter0.sw -t 16 -g 12 --partition=quick,norm --job-name msdss --module R)" >> msdss.sh
echo "jid_02=\$(sbatch --depend=afterok:\$jid_01 --mem=6g --partition=quick,norm --job-name=msdss /data/barbourcr/pipeline_for_meeting/scripts/msdss/results_batch_iter0)" >> msdss.sh

for j in {1..119}
do

i=$((j-1))

echo "jid_${j}1=\$(swarm -f /data/barbourcr/pipeline_for_meeting/scripts/msdss/swarmR_iter${j}.sw -t 16 -g 12 --partition=quick,norm --depend=afterok:\$jid_${i}2 --job-name msdss --module R)" >> msdss.sh
echo "jid_${j}2=\$(sbatch --depend=afterok:\$jid_${j}1 --mem=6g --partition=quick,norm --job-name=msdss /data/barbourcr/pipeline_for_meeting/scripts/msdss/results_batch_iter${j})" >> msdss.sh
done