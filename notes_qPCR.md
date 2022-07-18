Notes on qPCR
Prashant K

## General stuff

(_abandon_) Tried to generalize the 1-loading_files_funs to work with latest Quantstudio v2.6
- Chemistry field is missing -- make it work with the dye column instead
- Sample Setup sheet not exported, could just use Results
- Currently there is not much advantage to using the new software, just use V1.5.2 till these are changed

### Streamlining changes
- (_implement_) Make a conditional to plate layout reader if decimal point should be read vs biological replicates -- or (_alternatively_) Remove the biological replicates section and delimiter `.` from the plate layout. Use automatic inference like in the plate reader? 
	- _Problem with completely removing: the same template is used for more than 1 targets but the order is not consistent.., and we can't trust that the replicates assigned using row_numbers() will be in the correct order; relevant for ratio, normalization or for any matched statistical analysis_
  - _Use case:_ It makes it harder to type in *decimal dilution values* in the `assay_variable` spot. 
	  - Can make it backward compatible - could have a regex with `.[:digit:]$` detection in case replicates are specified :: How will this differentiate with a single decimal value after the point?

 
 Refactoring code to modularly work with both cq or linregpcr data..?
- [x] Need to generalize the `Sample Name` splitting and `Target Name` reassignment for TAQMAN to happen to the plate layout before merging with data - `Results` or `Amplification data` sheets
		- This might be an issue when not using assay mode. Since we are only doing assay mode and see no need for any custom mode, we will bother about bringing the whole sample name stuff into the assay mode `if` loop later

### Standard curve workflow
- [ ] Make a function to find the best fit set of dilutions to calculate std curve parameters - _by iteratively truncating the last 1,2,3 dilutions, with a tolerance of 0.990 R2.
- [x] Record the ID of the standard curve used in the `-processed` dataset for future lookup and easy reference. _output it into the html file for now_
- [x] Add a column for master mix in std curve data output (to be filled manually)
	- [ ] Can grab the master mix type from above the template in the future _(but is complicated when there is a mix of things on the plate)_  
- [x] Save the std curve into the bad std folder if rejected
	- [ ] Could also save the raw data, but add a column to mark std accepted or rejected
- [ ] How to work with multiple targets on plate which require multiple standard curve numbers? _Currently can rename them to the same number by hand and note in alt name their original name_
- [ ] (*feature*) Route rejected std curves into the bad standards folder
- [ ] _(Ignore the tiny value after decimal..)_: Decimal problem : decimals in Std curve quantity clashes with the biological replicate field. 
	- Solve by removing replicates from layout -- chcek if NA in replicates for multiple ones causes problems in the pipeline
	- Whatever comes into the biological replicate field should be added back to the Quantity with a decimal 

### features
- [x] Need to plot pseudolabels next to numbering for horizontal plot, ~~what is the best way to do this without replotting?~~ 
- [ ] Output qPCR std curve parameters at the begining of the html output
	- Here: ![[Pasted image 20220608105107.png]]

### Bugs
- [x] x axis title is being removed in plot_facetted_assay (example: _q32_copies-w line_)
	- [ ] Check if horizontal plots are ok now
- [ ] Check if the 4 ul template volume is being taken into account in std curve processing or regular sample processing. _Since labels say copies per ul template, this needs to be accurate_
- [ ] Make `.xaxis.label.custom` and optional variable in `12-plot_facetted_assay.R`


 ### Method/literature
 
- [ ] What is the problem with re-using a standard curve from a different run? _Considering the inter-run variation paper from Ruijter, they suggest a scaling between runs, which should preserve the efficiency from the standard curve, but the intercept could change?_
- [ ] what is the meaning of the intercept in the STD curve, why is it so variable across targets?

## RDML-linregPCR  

Verdict: **Linregpcr is abandoned for RAM project** due to baseline errors caused by early amplifications in 16s (high abundance) and no plateau for U64 probe assays
- Could try to use for SYBR dye assays though? (19/5/22) 

### Methodology questions

- When is linregpcr workable?
	- Plateau present
	- Not too early amplifications
	- Reasonable (?) baseline to plateau distance - [Tuomi, Jari Michael, et al. "Bias in the Cq value observed with hydrolysis probe based quantitative PCR can be corrected with the estimated PCR efficiency value." _Methods_ 50.4 (2010): 313-322.](https://www-sciencedirect-com.ezproxy.rice.edu/science/article/pii/S1046202310000538?via%3Dihub)
	> An important factor influencing the baseline estimation is the baseline-to-plateau distance (Δ_Rn_) in the raw, i.e. not baseline corrected, data. This distance is determined by intrinsic properties of the fluorescent reporter and on the primer and probe concentrations .. Δ_Rn_ values recorded with hydrolysis probes, however, are reported to be significantly lower due to the use of FRET based chemistries that display inefficient quenching. For these chemistries it is important to choose a reporter/quencher pair that gives the largest baseline-to-plateau difference. This will facilitate the estimation of the correct baseline value and thus reduce the variation in observed PCR efficiency values
	
	> . In optimized DNA dye-based assays, baseline fluorescence is below 1% of the fluorescence at the end of the PCR run; in probe-based assays the baseline may still be as high as 10% of the observed fluorescence. [Web based linregpcr, BMCbioinformatics, 2021](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04306-1)

Standards
- Can standards be combined post processing or does threshold need to be the same?
	- Standards will have their own efficiency (_expected due to DNA context, different contaminant concentrations etc._). So it needs to be analyzed separately anyways/the algorithm takes care if identified as calibrant? _I don't see any mathematical problem with using different threshold -- it is set after wol is identified, in the window of the linear region anyway right?__
- Should the standard calibrant be diluted or not? _dilution introduces errors, not diluting will increase chances of NTC for future runs_
	- Author: Dr. Maurice Van Der Hoff clarified that 1e3 copies of standard at 5 replicates needs to be used (let's use 6 like the paper says)
	- [Paper](https://academic-oup-com.ezproxy.rice.edu/clinchem/article/67/6/829/6247760?login=true)  recommends `undiluted calibrant` - which means use a single dilution in the qPCR..
	
	> PCR efficiency is affected by (_a_) systematic dilution errors in preparing a standard curve, (_b_) random pipetting errors in the standard curve samples, (_c_) the sequence and context of the target, and (_d_) unknown components inherent to the biological sample. Therefore, unbiased efficiency-corrected absolute quantification would benefit from a protocol in which **dilution of the standard is avoided** and the actual PCR efficiencies of the standard and unknown reactions are used in the calculations
	
	- mentions that you need `6` replicates for diluent of 1,000 copies from poisson stats -- that means it is diluted to a 1,000 right?
	- Higher concentrations also reduce subsampling error (poisson).

- [x] Does linregPCR only fit a single run? _Yes, it looks like that is the valid thing to do, there is a paper recommending overlapping standards to merge data from multiple runs_. What if I merge the calibration data from another experiment or run, how to make linregPCR process them together?

_tips for RDML:_ 
1. use `tojson()` to view RDML methods or outputs in plain text.
2. `fl.keys()` to see only the keys of all fields
3. use `fl.__getitem__('fieldname')` to quickly view subelements within the rdml fields -- works even in the subsets like target, dye etc..


### Quantstudio exported RDML format
- [x] Has dyes missing in the RDML file, _fixed in python using `fl.new_dye()` method_
- How to attach target_names in the RDML file: Currently do it manually in quantstudio
- [x] Do we need to set the `target chemistry` (shows up in the results file) to hydrolysis probe, since it is defaulting to `non-saturating DNA binding dye` all the time?
	- Not finding this option in RDML-Ninja
	- Is an attribute of `dye`, with name `dyeChemistry`. Set it to `"hydrolysis probe"` or `non-saturating DNA binding dye` otherwise. 
	- [x] More importanly, does this improve the LinregPCR analysis in any way? _rdml-tools help says it does_
- [ ] (_failed_)Also need to indicate `sample_type` in RDML to be ntc for negative control - this will remove errors about no amplification in those samples. clue: `cl.samples()[x].types()[0]['type']`
- [ ] (_last try_) Remove melt curve data `mdp` for probe data, if you want to open the file as valid rdml in online version. 
	-  Melt curve junk data is not exported if there is no SYBR on the plate, ex `q22 triplex` and `q24` but fails for `q25`  -- (I thought wrong -> )_Save rdml from quantstudio using `split items into individual files` option, this will prevent writing junk values to the melt data_
	- Expt -> run -> react -> data -> mdp. It is hard to locate react data onwards. clue: `_get_all_children(self._node, "react")`

**Problems with rdmlpython** : Is giving baseline error in a few samples, I cannot figure out the reason. 

### amplification analysis, qc
Works for S019_25-11-19 file (SYBR) and q25_S037_RAM repression_14-2-22 (Probe) with recent changes pulled on _18/2/22_
> (old, before pulling) Does not work for q25 (TAQMAN) file with `rawFluor[rowCount, cyc] = float(fluor) ; IndexError: index 72 is out of bounds for axis 0 with size 72`

Probe based assays with late amplification (ex: `U64`) is showing very low efficiency (_due to plateau phase not present_) -- Is there some tweak that can rescue the analysis?
- [ ] check if the data makes sense, low efficiencies are reliable
	- quick amplifications of 16s are not being analyzed due to `baseline error`, _asked on the github page if there is a straightforward fix_ 
	- Considering that we cannot dilute the samples, verdict will be to ABORT LINREGPCR mission 
	
	Improving analysis for flagged samples
	- [ ] How to troubleshoot baseline error due to very early amplifications?
	- [ ] I see error `Cq < 10; N0 unreliable` while analyzing q12 data for 16s. Why is this unreliable? 
	
Understanding baseline error
![[baseline-error_q20b.png]]
Othwerise very similar curves with early amplification have problems with baseline identification if they lack the initial ~ flatness. 
They also inaccurately say `no plateau` since that is linked to finding a baseline?

### python `linregpcr-analyze.py`
- [ ] directory issues : _Sometimes re-running the script fixes it, if not restart spyder_. Error: `FileNotFoundError: [Errno 2] No such file or directory: 'RDML_files/q27_328 330_Std27_9-3-22.rdml'`. There is something wrong with the script's interaction with directory and files. 
	- The script works for one rdml file but does not find the file when run again or with other rdml files
	- Any new files added after running the script are not found. Calling `os.chdir` inside the function seems to freeze it's directory state at when the script was run 
	- [ ] Try running this from the outside directory., after moving the g14 script
	
### melt curve analysis

Melt curve does not return anything for the S019_25-11-19 file both in python and online
- Identified that this could be since the meltingTemperature for each target was not present in the rdml file - `fl.targets()[1].tojson()['meltingTemperature']`. This is not shown in the rdmlninja. Verified with `sample-online.rdml` example file downloaded from online rdml website
- [ ] _Read from an excel sheet of target-meltingtemps and add it in using python with a switch-case kind of statement_
	- Use `fl.targets()[0]__setitem__('meltingTemperature', '70.0')`
- meltingTemperature inside `target` only seems to be available in `rdml 1.3` ☹. Tried `fl.migrate_version_1_2_to_1_3` to fix

Minor things
- [ ] Output processed linregpcr data from `linreg-post-processing.R` script
- [ ] Add N0 = 'Copies per sample, extrapolated' to labelling helper


### R-python integration
- Make the test_linregpcr.py into a function and call it using `reticulate::source_python("...py")` [documentation](https://rstudio.github.io/reticulate/#sourcing-python-scripts)


### /obsolete/: RDML-R-attach names-export
Analysis needs _well names_, and _target names_. Attach the names from the plate layout just before plotting. To get target names for each well for SYBR assays -  _Name manually in Quantstudio_ until more automation is explored
- R's RDML package's `$AxXML` exported `.rdml` file with  cannot be validated by rdmlpython and shows as almost empty file in RDML-ninja
	- Looks like the actual data is being written to a temporary folder : `C:\Users\new\AppData\Local\Temp`
	- View AsXML() method in the RDML class as a learning experience
- _ignore_/ Need to get R's RDML into a table format, so that sample names from the google sheet can be merged by well
  - How do you change the sample names in both the data `($getFdata)` and metadata `($sample)`


## Sliwin using R
_This was abandoned in favour of the rdmlpython workflow since there were lot of problems with fitting when no amplifications happen etc. All this effort was already done by the rdml folks in python_

### pcrfit()
[x] figure out the input format for pcrfit; why is model not working without loading qPCR package..
	- [x] qPCR is blocking `dplyr::select()` : _just override select_
- Negative Rn values cause error (ex: E10 in q22)
	 `qpcr_fit = map(data, ~pcrfit(.x, cyc = 2, fluo = 3)) x 0 (non-NA) cases`
	- How to filter automatically or figure out the cause 
	-[x] _Can have an if_else switch inside the mapping function, with a flag column..?_

### sliwin()	
- Sliwin throws same error and stops at samples that don't amplify. very dumb behavior
	- [ ] identify a variable to flag using fit variablts d (min) and c(max) - if they are too close to each other or some other arbitrary critirea ☹
- Also the Rn values need to be background subtracted before values are found, **does sliwin do this?** 
Abandon and move to rdmlpython until this issue is figured out. There seems to be too many things needed to make sliwin work - since it is a function written by third party folks, they didn't put that much effort.

_26/2/22 : **Update**: re-visiting sliwin to verify if the rdmlpython individual/group efficiencies for U64/flipped are justifiably very low (~ 1.4), since sliwin shows curves and it's easier to get things to plot with R for me_ 
- Figure out how rdmlpython treats things as no amplification. Follow `vecNoAmplification` around line 11030

		# Check to detect the negative slopes and the PCR reactions that have an
		# amplification less than seven the minimum fluorescence
				if slopeAmp[oRow] < 0 or minSlopeAmp[oRow] < (np.log10(7.0) / minFluCountSum[oRow]):
					vecNoAmplification[oRow] = True
		-----------------
		# There must be an increase in fluorescence after the amplification.
		did something where mean of fluors check to increase by 1.2, 
		cycle 8,9 / cycle 0,1 < 1.2
		
- Decided on an arbitrary cutoff of 3 fold of final / initial signal to call amplification (for probe q25 data)
- sliwin works now, but unable to augment fitted data --
- just create a line with slopes taken from the excel file parallel to initial data to compare?

-------------
----- old stuff ----------

## 17/1/22
- Practicing loading rdml files into 
[x] rdmlpython : Validation fails - 

```
Schema validation result:	False	RDML file is not valid.
Schema validation error:	False	
Line 110, Column 0: Element '{http://www.rdml.org}dyeId': 
No match found for key-sequence ['SYBR'] of keyref '{http://www.rdml.org}dyeKeyRef'.  
```

Similar error for FAM and VIC in q20 data

- Unable to view the raw RDML file in notepad/libre calc..
- online tools also don't show the dyeID field
- saving the file from R AsXML cannot be read as RDML by rdmltools

[x] RDML - R : dyes list is empty, validation function does not exist in the package

- Study the schema to understand the error. A snippet from the schema documentation [RDML V1.2](http://rdml.org/files/rdml/RDML_v1_2_REC.xsd)

```
<xs:key name="documentationKey">
<xs:selector xpath="./rdml:documentation"/>
<xs:field xpath="@id"/>
</xs:key>
<xs:keyref name="**dyeKeyRef**" refer="rdml:dyeKey">
<xs:selector xpath="./rdml:target/**rdml:dyeId**"/>
<xs:field xpath="@id"/>
</xs:keyref>
<xs:key name="dyeKey">
<xs:selector xpath="./rdml:dye"/>
<xs:field xpath="@id"/>
</xs:key>
<xs:key name="experimentIdKey">
<xs:selector xpath="./rdml:experiment"/>
<xs:field xpath="@id"/>
```

## 21/1/22

- Figured out that this is the only missing line for each dye, and can be added in RDML-ninja software
`<dye id="SYBR"/>`
[ ] Find out what setting is missing in quantstudio to get this automatically
[ ] see if there is a commandline way to add this line (for all dyes) into the rdml file

Found new versions of quantstudio 1.5.2 and 2.6
- v1.5.2 output RDML has the same validation error as the older 1.5.0
- v2.6 : new validation error: Not writing `description` and `type`  fields under each target. Error : `Element dyeID is not defined in this scope`
    - can't figure out if there is some field that need to be filled in the Quantstudio software that can add these fields in the exported RDML file

```
<target id="flipped">
  <description>None</description>
  <type>toi</type>
  <dyeId id="SYBR"/>
```


- Melt curve analysis on the online portal does not work. Try in the python version. Could be that some parameters need to be changed

## 26/1/22 : exponential fitting done in 2G

2G_stability data : did exponential fitting
- change normalization so that the mean at first data point is 1 for all targets
- Get rate constants from the data and put it on the plot
- Check if the fitting can work with fixed parameters of initial and final asymptotes -- noticed singular gradient error..
