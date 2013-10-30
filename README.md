MaLTE
=====
Machine Learning of Transcript Expression

News:
0.2-5: EXPERIMENTAL (forked as 'quantreg' branch)
- Incorporated quantregForest for interval estimation of predictions
- Incorporated gene-specific tuning; modified TT.Params to have a slot for the
	tuned OOB Pearson correlation for comparison (without the p-value)
- Fixed a bug that affected OOB filtering
- Delete models immediately after tuning: does this save memory?

0.2-4:
- OOB filter for genes now takes into account significant p-values (<=0.05)

0.2-3:
- Bug: '*NA' in samples.txt fails (FIXED)

0.2-2:
10-10-2013
- Bug: avoid computing correlations when the number of test samples <=2

0.2-1:
10-09-2013
- Added a utility function to estimate within-sample correlations + 
	documentation

0.2:
08-09-2013
- Added TT.Seq.[Gene/Tx] methods: 'cor.P()' and 'cor.S()' to get Pearson and 
	Spearman correlations with RNA-Seq expression
- Added utility functions: 'cors.P()' and 'cors.S()' to collate vectors of 
	correlation
- Added utility funtions to compare MaLTE with microarray summarisation results:
	'compare.correlations()'
- Added necessary documentation

0.1-8:
03-09-2013
- Experimental: modified 'get.predictions()' to work on non-filtered tt.seq 
	lists
- Experimental: added a 'get.trues()' function that can work on non-filtered 
	tt.seq.lists

0.1-7:
22-08-2013
- made 'get.predictions()' much faster using unlist( list, F, F ) and *apply 
	functions instead of a for loop

17-08-2013
- fixed a bug that prevented 'get.predictions()' for transcript isoform 
	expression with more than one isoform filtered (DONE)

0.1-6:
15-08-2013
- refactored 'prepare_txs_data.py'
- fixed a bug in the OOB filter function that prevented filtering of genes with 
	single transcripts that didn't pass OOB
- added a feature to include transcript correlations in the rare event that both 
	RNA-Seq and array data is available for test samples

16-08-2013
- fixed the above bug on filtering transcripts

0.1-5: 
14-08-2013
- refactored 'prepare_data.py'
- corrected stats produced by 'prepare_data.py'

0.1-4:
13-08-2013
- included parallel execution to create resource files used to build training 
	and test data

0.1-3:
05-08-2013
- major fixes to documentation except documenting datasets



Checklist:
- 80 characters
- comments
- replace 'NA' by NULL appropriately
- make sure user cannot enter samples.txt with samples that don't exist
- make sure that the return status of the python script are appropriately set 
	under all conditions of either success or failure
- a data package?
- what if you try to filter and the other is not oob? (DONE)
- consistency in show() methods for all classes (DONE)
- handle warnings to do with correlation computations e.g. 
	run( tt.ready.txs[[10]], tt.params )
- add keywords in the documentation (DONE)
- replace 'T' with 'TRUE' (DONE)

Feature list:
- leave-one-out cross-validation
- direct comparison to median-polish et al. (DONE)
- add ComBat stuff
- make prepare_data parallel (DONE)


- Later: provide a graphical description of which the best-performing probes are
	 (I'm thinking of integrating the gene annotation model information together 
	with a color-key-map that shows the reliablity of a probe in say a particular 
	tissue)
