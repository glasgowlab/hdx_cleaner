Command line inputs:

--t, --table: Path(s) to input csv. Required
--o, --out: path to output csv. Default: 'CLEAN.csv'
--p, --plots: directory for output plots. Default: './csvplots'
--ov, --overlap: flag specifying how to treat overlapping duplicates. Options: 'keep', 'drop', 'select' (or any other). Default: 'select'
--r, --rangeslist: path to output rangeslist.
--RTC: RT cutoff for duplicate peptides. Default: 0.5
--MZC: m/z cutoff for duplicate peptides. Default: 0.1
--scoreC: Score threshold for provisional cut for trendline fit. Default: 150
--ppmC: m/z deviation cutoff (in ppm) for cut about trendline. Default: 7
--maxfev: maxfev for curve_fit for trendline. 

Default inputs:

--t (list of peptide pool csvs)
--b (table for batched data. Each row contains 'name' string for protein, then all pool csvs corresponding to that pool

Default outputs:

--o (batched and cleaned peptide pool csv): CLEAN.csv
--r (rangeslist table for merge in post-PIGEON)
in --plots:
	score histograms:
		/hist-scores-all.png: histogram of score for all pooled data
		/hist-scores-cut.png histogram of score for high scoring data used for curve fit
		/hist-scores-ppm.png histogram of score after cut at ppmC from trendline
		/hist-scores-topscores.png histogram of score after dropping all but best match for each charge state
		/hist-scores-flagged.png histogram of score for non-overlapping duplicates
		/hist-scores-yellow.png histogram of score for overlapping duplicates
		/hist-scores-clean.png histogram of score for final cleaned data
	scatterplots -- ppm m/z deviation vs m/z measured, colorbar Score:	
		/scatter-30ppm.png: all data before cuts
		/scatter-hscoring.png: high scoring data used for curve fit
		/scatter-7ppm.png: all data after cut at ppmC from trendline
		/scatter-topscore.png: data after dropping all but best match for each charge state
		/scatter-flagged.png: data dropped for being duplicates, non-overlapping
		/scatter-yellow.png: data dropped for being duplicates, overlapping
		/scatter-clean.png: final cleaned data
	ppm histogram -- hist of abs(distance from trendline) for all data, showing location of cut
		/hist-ppms.png

Example input: python biotools-cleaner.py --t pool-1.csv pool-2.csv pool-3.csv --p csvplots --ov keep --o keep-pooled.csv --r rangeslist.csv --MZC 0.05 --RTC 0.25 --scoreC 200 --ppmC 5 --maxfev 2000
