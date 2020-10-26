if [ -e diagnosis.sh ]
then
rm diagnosis.sh
fi

touch diagnosis.sh
echo "jid_01=\$(swarm -f /home/wangq13/pipeline_for_meeting/scripts/diagnosis/swarmR_iter0.sw -t 16 -g 12 --partition=quick,norm --job-name diagnosis --module R)" >> diagnosis.sh
echo "jid_02=\$(sbatch --depend=afterok:\$jid_01 --mem=6g --partition=quick,norm --job-name=diagnosis /home/wangq13/pipeline_for_meeting/scripts/diagnosis/results_batch_iter0)" >> diagnosis.sh

for j in {1..119}
do

i=$((j-1))

echo "jid_${j}1=\$(swarm -f /home/wangq13/pipeline_for_meeting/scripts/diagnosis/swarmR_iter${j}.sw -t 16 -g 12 --partition=quick,norm --depend=afterok:\$jid_${i}2 --job-name diagnosis --module R)" >> diagnosis.sh
echo "jid_${j}2=\$(sbatch --depend=afterok:\$jid_${j}1 --mem=6g --partition=quick,norm --job-name=diagnosis /home/wangq13/pipeline_for_meeting/scripts/diagnosis/results_batch_iter${j})" >> diagnosis.sh
done