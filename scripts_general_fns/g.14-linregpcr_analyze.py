# -*- coding: utf-8 -*-
"""
Playing with rdml file input output : rdmlpython package
https://github.com/RDML-consortium/rdmlpython

Created on Mon Jan 17 15:20:18 2022

@author: Prashant Kalvapalle
inspiration # from : http://rdml.org/referenceImplement.html
"""

def linregpcr_analyze(filename_rdml = 'S019_25-11-19_v152'):

    import rdmlpython as rdml # load package
    
    # User input
    #filename_rdml = 'S019_25-11-19_v152'
    # directory -- relative might be an issue..
    dirname = 'RDML_files/'
    
    flpath = dirname + filename_rdml + '.rdml'
    
    # Read file
    fl = rdml.Rdml(flpath)
    
    
    # Fix the dyes issue for validation purposes
    
    # get all dyes used in the RDML file
    dyeIds = set() # initialize empty set. Prevents duplicates so better than list
    for targetfield in fl.targets():
        dyeIds.add(targetfield.tojson()['dyeId'])
    
    # add all these dyes into the rdml file
    for dye in dyeIds:
        fl.new_dye(dye)
    
    # Validate a file
    # print(fl.validate())
    # doesn't validate here but works in RDML-ninja
    
    # Run LinRegPCR
    
    # unpack the rdml file, 1st experiment, 1st run
    first_run_data = fl.experiments()[0].runs()[0]
    
    # amplification data - run linregpcr
    qpcr_result = first_run_data.\
        linRegPCR(updateRDML=True,
                  saveResultsCSV=True, 
                  timeRun=True, 
                  verbose=True)
        
    # run melt curve only if dye chemistry (FAM or SUN absent in dyelist)
    if('FAM' or 'SUN' not in dyeIds):      
        melt_curve_result = first_run_data.\
            meltCurveAnalysis(updateRDML=True,
                              saveResultsCSV=True,
                              verbose=True)
     # melt curve not working with default parameters -- same as online   
    
    # save results to RDML
    fl.save(dirname + "results/" + filename_rdml + "fixpython.rdml")
    
    # save amplification analysis results as csv file
    with open('excel files/' + filename_rdml + ".csv", "w") as filehandle:
        filehandle.write(qpcr_result["resultsCSV"])
            
    # save melt curve analysis results as csv file: if a non empty result exists
    if('melt_curve_result' in locals() and len(melt_curve_result)):  
        with open('excel files/' + filename_rdml + "-melt.csv", "w") as filehandle:
            filehandle.write(melt_curve_result["resultsCSV"])

