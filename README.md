# WW-CoV2-project
Takes excel file output from qPCR and ddPCR and makes neat plots with appropriate labels and calculations using metadata from other google sheets'

## Git organization
1. Each branch of this repo is the code run on a file named after the experiment (ex: WW21) with the desired plotting formats and extra items specific to the run
2. As the methods become standardized, all data can be run just from the master branch without making new branches

## Readme step by step guide
### Current workflow for qPCR data analysis

#### Sample naming in the 'calculations (lab notebook)' google sheet
	1. Set sample names using the automatic labeller to the right side of the plate
		a. Target_Sample category-tube name
		
#### Quantstudio
	1. Open the .eds file in Quantstudio (applied biosystems)
	2. Check amplification curves if everything seems right
		a. Check for systematic amplifications in non technical controls (NTC) - Or stochastic amplifications that seem problematic; by looking at raw data (multicomponent plots)
	3. Check that the threshold is at the desired value for each target (by selecting each target): Setting thresholds to 0.04 right now  
	4. Export excel (.xls) file from Quantstudio with the same name as the qPCR file (maybe exclude the date)

#### Sourcetree (version control)
	1. Navigate and doubleclick on the master branch
	2. Click on BRANCH and make a new branch with the experiment name (shortened)
	3. Now you automatically are on the new branch (in R and your actual folder on windows as well)

<Make sure you mirror the directory structure for the excel files and qPCR analysis folders>

#### Rstudio
	1. Open the Rproject file on Rstudio - this will load from the current directory 
	2. If running a standard curve on the same plate, open the Standard_curve.R file; change the file name and title for the plot and run
	3. Open the inputs_for_analysis script
		a. Change the file name (flnm) to the excel file name
		b. Verify and change the template_volume loaded in the qPCR well
		c. Check the standard curve parameters (std_par) and update it from the equation on the plot produced by the standard_curve.R
		d. Save file
	4. Open analysis.R file
		a. Source the file (click save with the check mark on source on save option
		b. This produces a plot (not saved) and dumps the data with the appropriate sample labels in a google sheet 'qPCR data dump'
	5. Open 'make_html_qPCR.Rmd' file
		a. Change the read_these_sheets to the 1 or more sheets to plot together in the 'qPCR data dump' and change the title_name
		b. Check the HA concentration factor, and spike virus conc variable to match what was used this week
		- Spike virus conc should be manually taken from the qPCR run that measures it - from qPCR data dump file
		c. Change bb_sheets to the appropriate week(s) from the biobot Sample IDs and check that the names of the last 3 columns match that of previous week (with a successful run)
		d. If you want to exclude any particular categories in the plot enter at extra_categories (then you can make separate plots for those categories by changing the switch exclude_sample (control + F for this variable in the Rmd code)
		e. Click Knit
In case you see any errors, follow the line number in the error and do extensive google searches - That's how I learn't it
