Notes on qPCR
Prashant K

## General stuff

Tried to generalize the 1-loading_files_funs to work with Quantstudio v2.6
- Chemistry field is missing -- make it work with the dye column instead
- Sample Setup sheet not exported, could just use Results
- Currently there is not much advantage to using the new software, just use V1.5.2 till these are changed
 
[x] Need to generalize the `Sample Name` splitting and `Target Name` reassignment for TAQMAN to happen to the plate layout before merging with data - `Results` or `Amplification data` sheets
	- This might be an issue when not using assay mode. Since we are only doing assay mode and see no need for any custom mode, we will bother about bringing the whole sample name stuff into the assay mode `if` loop later

## Sliwin work
### pcrfit()
[x] figure out the input format for pcrfit; why is model not working without loading qPCR package..
	- [x] qPCR is blocking `dplyr::select()` : _just override select_
- Negative Rn values cause error (ex: E10 in q22)
	 `qpcr_fit = map(data, ~pcrfit(.x, cyc = 2, fluo = 3)) x 0 (non-NA) cases`
	- How to filter automatically or figure out the cause 
	- _Can have an if_else switch inside the mapping function, with a flag column..?_
- Sliwin throws same error and stops at samples that don't amplify. very dumb behavior
	- [ ] identify a variable to flag using fit variablts d (min) and c(max) - if they are too close to each other or some other arbitrary critirea ☹
Abandon and move to rdmlpython until this issue is figured out. There seems to be too many things needed to make sliwin work - since it is a function written by third party folks, they didn't put that much effort.

## RDML goals 

[ ] Need to get R's RDML into a table format, so that sample names from the google sheet can be merged by well
  - How do you change the sample names in both the data ($getFdata) and metadata ($sample)
  - Couldn't validate the RDML file saved from R using `$AsXML` as [documented](https://pcruniversum.github.io/RDML/articles/RDML.html) : So this won't open in python. It means that sample names cannot be attached in R before handing it off to python in another RDML file
  	- Could figure out drawing the sample names and attaching it in python
  	- Could try to get rdmlpython to work using a dataframe instead of an RDML file

### RDML-R-attach names-export
- R's RDML package's `$AxXML` exported `.rdml` file with  cannot be validated by rdmlpython and shows as almost empty file in RDML-ninja
	- Looks like the actual data is being written to a temporary folder : `C:\Users\new\AppData\Local\Temp`

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
