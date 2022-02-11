# -*- coding: utf-8 -*-
"""
Playing with rdml file input output : rdmlpython package
https://github.com/RDML-consortium/rdmlpython

Created on Mon Jan 17 15:20:18 2022

@author: Prashant Kalvapalle
inspiration # from : http://rdml.org/referenceImplement.html
"""

import rdmlpython as rdml

# User input
filename_rdml = 'S019_25-11-19_v152'
dirname = 'RDML_files/'

flpath = dirname + filename_rdml + '.rdml'

# Read file
fl = rdml.Rdml(flpath)


# Fix the dyes issue for validation purposes
# hint : fl.dyes()[1].tojson()['id']

# get all dues used in the RDML file
dyeIds = set() # initialize empty set. Prevents duplicates so better than list
for targetfield in fl.targets():
    dyeIds.add(targetfield.tojson()['dyeId'])

# add all these dyes into the rdml file
for dye in dyeIds:
    fl.new_dye(dye)

# Validate a file
print(fl.validate())
# doesn't validate here but works in RDML-ninja

# Run LinRegPCR

# unpack the rdml file, 1st experiment, 1st run
qpcr_result = fl.experiments()[0].runs()[0].linRegPCR(updateRDML=True, saveResultsCSV=True, timeRun=True, verbose=True)
    

# save results to RDML
fl.save(dirname + "results/" + filename_rdml + "fixpython.rdml")

# save results as csv file
with open('excel files/' + filename_rdml + ".csv", "w") as cli_f:
    cli_f.write(qpcr_result["resultsCSV"])
        
