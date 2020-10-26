# Random Forest Pipeline for NIH HPC

This is a random forest pipeline create for the NIH HPC enviroement. The ratios between all features (proteins) is created. These ratios are then concatednated to the protein abundance table. This list of proteins and ratios in then read into the ranger random forest. The AUC is calculated and the bottom 10% of features are then removed until only one feaure is left. Each round 10 random seeds are used to create 10 models. The aveage AUC is then calculated.
