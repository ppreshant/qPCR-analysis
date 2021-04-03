# qPCR-analysis
Takes excel file output from qPCR and plots CT curve and Tms in publication format

*Rationale*: Coding this trivial analysis and plotting is especially relevant because I run tons of qPCR and I want all the plots to have the same formatting for easy comparision of data between runs. I wouldn't want to plot all these on excel!

**Sample nameing scheme**
Sample names are written in a google sheet (contact me for access to see the template) in this particular format
`primerpairname-overall name_templatename.biological replicate number`
overall name could stand for test, positive, negative (controls) or NoRT control etc.
Examples: old16sU64-Test_295.1,	old16sU64-NoRT_89.1
![image](https://user-images.githubusercontent.com/14856479/113488074-6cae8200-9481-11eb-9d82-e97033b72e2e.png)

The output image looks like this
- primerpairname becomes facets or panels; 
- template on the x axis and 
- overall name as colour
![image](https://user-images.githubusercontent.com/14856479/113488285-a03ddc00-9482-11eb-810a-4756b09fc10e.png)


**Git organization**

1. *main* branch holds the latest version of the code for general purpose plotting for any qPCR output file with sample names in the above way
2. For customized plots for individual experiments - there will be an ad-hoc branch with specific code files that modify these base plots. Or for publications there might be a different branch named after the experiment (ex: S024) with the highly customized plotting formats and extra items/tables output for supplmentary data (this will make it easy to provide the code along with the publication)
3. There are different main branches for each of the major tasks the experiments fall under a. Time_series_master (ex: S019) b. master (ex: S024) (for generic qPCR) C. Dose_response (ex: S03) : If branches with these names don't exist yet, you can start with the branches named after the experiments
