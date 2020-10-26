#!/bin/sh

# Job Name
#$ -N phos_distant
#$ -cwd
#$ -j y
#$ -l h_vmem=8G
#$ -V
#$ -S /bin/sh
#$ -o Result/
#$ -pe threaded 12

# Send mail when the job is submitted, and when the job completes
#$ -m abe

#  Specify an email address to use
#$ -M isaac.bishof@nih.gov

module load mono/5.12.0.233-foss-2016b


mono /hpcdata/bcbb/bishof/UQ_67/MaxQuant/bin/MaxQuantCmd.exe  mqpar.xml
