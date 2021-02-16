Guide through the analysis performed on the data obtained from our untargeted metabolomics analysis. 

Raw data = bee_August2017_strain_specificity data.xls

#Initial file conversions

Raw intensities data was saved from excel-sheet in tab-delimited format. The data is spread over two sheets, so there are two files:

raw_intensities_sheet1.txt
raw_intensities_sheet2.txt

For futher processing, the backslash should be removed, and the first three
headers are also not needed at this point:

								cut -f4- raw_intensities_sheet1.txt > techrep_sheet1.txt
	sed 's/ \/ /_/g' techrep_sheet1.txt > temp
	mv temp techrep_sheet1.txt

	cut -f4- raw_intensities_sheet2.txt > techrep_sheet2.txt
	sed 's/ \/ /_/g' techrep_sheet2.txt > temp
	mv temp techrep_sheet2.txt


#Calculating mean ion intensity for technical replicates, order the data

This was done using the perl script parse_intensities.pl.

Usage:

	perl parse_intensities.pl raw_intensities_sheet1.txt raw_intensities_sheet2.txt > mean_intensities.txt

From here we started to process the data using the R script = metabolomics_untargeted.R


#Technical reproducibility

The intensity values corresponding to the two injections for each sample are plotted against each other. 

Overall the reproducibility looks good so no sample was excluded from the analysis.  


#Determine which ions are pollen-derived

This is particularly important for the ions that will be considered as substrates (log2FC<=-2 between the 16h time-point and the 0h time-point). 

We included in our untargeted metabolomics analysis 3 replicates (PE_A, PE_B, PE_C) of a dilution series of our pollen extracts (2x dilutions for 10 times). 

We plotted each ion intensity along the dilution series and we extrapolated the R2 of each plot. We calculated the log2FC between the undiluted pollen samples (PE_1X) and the H2O injections. We plotted the R2 vs the log2FC and we set a threshold to define pollen-derived ions = log2FC>=2 and R2>0.75. 

We then started to calculate the log2FC between the two time-points for each strain
#T-tests, fold-change and Volcano plots for individual strains

Generation of files containing log2FC and p-value information:

F5_183_FC.txt
F5_184_FC.txt
F5_185_FC.txt
F5_186_FC.txt

We combined these files into the table all_log2FC.xlsx and used them to produce volcano plots.

We selected the ions that significantly changed in intensity (p<0.01) between the two time-points in at least one strain with log2FC values <=-1 or >=1 and we visualised them using dotplots.

Finally we also used the log2FC values to plot a PCA. 







