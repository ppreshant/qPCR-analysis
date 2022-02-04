# -*- coding: utf-8 -*-
"""
Playing with rdml file input output : rdmlpython package
https://github.com/RDML-consortium/rdmlpython

Created on Mon Jan 17 15:20:18 2022

@author: Prashant Kalvapalle
"""

import rdmlpython as rdml

# User input
filename_rdml = 'S019 25-11-19'
dirname = 'excel files/'

flpath = dirname + filename_rdml + '.rdml'

# Read file
fl = rdml.Rdml(flpath)


# Validate a file
print(fl.validate())
# from : http://rdml.org/referenceImplement.html
    
# Run LinRegPCR

# Need to edit this
cli_linRegPCR = rdml.Rdml("test_2_raw_data.rdml")
cli_expList = cli_linRegPCR.experiments()
cli_exp = cli_expList[0]
cli_runList = cli_exp.runs()
cli_run = cli_runList[0]
cli_result = cli_run.linRegPCR(
    updateRDML=True, saveResultsCSV=True, timeRun=True, verbose=True)
cli_linRegPCR.save("result.rdml")
with open("result.csv", "w") as cli_f:
    cli_f.write(cli_result["resultsCSV"])
