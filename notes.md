Notes on qPCR
Prashant K

## Goals 

[ ] Need to get R's RDML into a table format, so that sample names from the google sheet can be merged by well
  - How do you change the sample names in both the data ($getFdata) and metadata ($sample)

## 17/1/22
- Practicing loading rdml files into 
[x] rdmlpython : Validation fails - 

``` Schema validation result:	False	RDML file is not valid.
Schema validation error:	False	
Line 110, Column 0: Element '{http://www.rdml.org}dyeId': 
No match found for key-sequence ['SYBR'] of keyref '{http://www.rdml.org}dyeKeyRef'.``` 

Similar error for FAM and VIC in q20 data

- Unable to view the RDML file's like 110 in notepad/libre calc..
- online tools also don't show the dyeID field
- saving the file from R AsXML cannot be read as RDML by rdmltools

[x] RDML - R : dyes list is empty, validation function does not exist in the package

- Study the schema to understand the error. A snippet from the schema documentation [RDML V1.2](http://rdml.org/files/rdml/RDML_v1_2_REC.xsd)

```<xs:key name="documentationKey">
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
<xs:field xpath="@id"/>```

## 21/1/22

- Figured out that this is the only missing line and can be added in RDML-ninja software
`<dye id="SYBR"/>`
[ ] Find out what setting is missing in quantstudio to get this automatically
[ ] see if there is a commandline way to add this line (for all dyes) into the rdml file

- Found new versions of quantstudio 1.5.2 and 2.6
    - Tested v2.6 : Not writing description and type fields causing validation error
    - can't figure out if there is some field that need to be filled in the software..

```<target id="flipped">
  <description>None</description>
  <type>toi</type>
  <dyeId id="SYBR"/>```


- Melt curve analysis on the online portal does not work. Try in the python version. Could be that some parameters need to be changed
