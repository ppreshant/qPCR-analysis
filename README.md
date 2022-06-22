# qPCR-analysis
Takes excel file output from qPCR, sample names layout and optional standard curve parameters. Plots Cq data, absolute concentrations and Tms in publication format.

_standard curve_ : Include Stdxx (xx = numbering) in the excel file name to process standard curve from data within the file

*Rationale*: Coding this trivial analysis and plotting is especially relevant because I run tons of qPCR and I want all the plots to have the same formatting for easy comparision of data between runs. I wouldn't want to plot all these on excel!

## How to run
1. Clone the respository using git onto your computer. This brings in the code files and also mimics the folder structure (check if all the necessary folders are present)
2. Open the scripts in R by clicking on the `qPCR analysis.Rproj` file. _This ensures that R's current working directory is at the base of this project. 
	- If you haven't you need to install R, Rstudio, git; and install the packages listed with the library commands in `0-general_functions_main.R` script using `install.packages('..')`_
3. *Set metadata sheets*: Make a googlesheet with sheets named `qPCR plate layouts` and `qPCR Std curves` following the template from [my sheet](https://docs.google.com/spreadsheets/d/1RffyflHCQ_GzlRHbeH3bAkiYo4zNlnFWx4FXo7xkUt8/edit#gid=0)
	- Copy the url of your googlesheet and replace the one in `sheeturls = list(plate_layouts_PK = ..` in the script `0-general_functions_main.R`
4. Open `analysis.R` script. Change the `flnm` and `title_name` variables. and other options
5. You can try running the template file `q33_RNA stability-3_20-6-22` within the `excel files` folder by copying the template data for q33 and Std30 standard data from [my sheet](https://docs.google.com/spreadsheets/d/1RffyflHCQ_GzlRHbeH3bAkiYo4zNlnFWx4FXo7xkUt8/edit#gid=0)

**Sample naming scheme**
Sample names are written in a google sheet (contact me for access to see the template) in this particular format
`targetname-samplecategory_templatename.biologicalreplicatenumber`
- targetname = name of the amplicon, _ex: 16s, U64, spliced_
- samplecategory could stand for test, positive, negative (controls) or NoRT control etc.
- templatename = name of the DNA or RNA sample loaded as the template
- biologicalreplicatenumber = 1,2,3 ..
Examples: old16sU64-Test_295.1,	old16sU64-NoRT_89.1

![image](https://user-images.githubusercontent.com/14856479/113488074-6cae8200-9481-11eb-9d82-e97033b72e2e.png)

The output image looks like this
- primerpairname becomes facets or panels (boxes on the top); 
- template on the x axis (named `assay_variable`) and 
- samplecategory as colour (named as `Sample_name`)
<img src = 'https://user-images.githubusercontent.com/14856479/113488826-1859d100-9486-11eb-8384-1ad17afea737.png' width = "600">



**Git organization**

1. *main* branch holds the latest version of the code for general purpose plotting for any qPCR output file with sample names in the above way
2. For customized plots for individual experiments - there will be an ad-hoc branch with specific code files that modify these base plots. Or for publications there will be a different branch named after the experiment (ex: S024) with the highly customized plotting formats and extra items/tables output for supplmentary data (this will make it easy to provide the code along with the publication)
3. There are different main branches for each of the major tasks the experiments fall under a. Time_series_master (ex: S019) b. master (ex: S024) (for generic qPCR) C. Dose_response (ex: S03) : If branches with these names don't exist yet, you can start with the branches named after the experiments
