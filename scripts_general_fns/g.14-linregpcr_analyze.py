# -*- coding: utf-8 -*-
"""
Playing with rdml file input output : rdmlpython package
https://github.com/RDML-consortium/rdmlpython

Created on Mon Jan 17 15:20:18 2022

@author: Prashant Kalvapalle
inspiration # from : http://rdml.org/referenceImplement.html
"""

def linregpcr_analyze(filename_rdml = 'S019_25-11-19_v152'):

    # Module for running as a python script (for debugging)    
    import os
    os.chdir('..') # change to the root directory : '/qPCR'

    # Import rdml
    import rdmlpython as rdml # load package
    
    # User input
    #filename_rdml = 'S019_25-11-19_v152'
    # directory -- relative path might be an issue..
    dirname = 'RDML_files/'
    
    flpath = dirname + filename_rdml + '.rdml'
    
    # Read file
    fl = rdml.Rdml(flpath)
    
    # Change rdml version 1.2 to 1.3 (enables changing dyeChemistry)
    fl.migrate_version_1_2_to_1_3()
    
    
    
    # Fix the dyes issue for validation purposes  
    
    # get all dyes used in the RDML file
    dyeIds = set() # initialize empty set. Prevents duplicates so better than list
    for targetfield in fl.targets():
        dyeIds.add(targetfield.tojson()['dyeId'])
    
    # create a dictionary to check the chemistry type based on the dyes present
    # determines if melt curve was done or not and if plateau exclusion is lenient
    dictChemistry = {
        **dict.fromkeys(['FAM', 'SUN', 'VIC', 'TAMRAC'], 'hydrolysis probe'),
        **dict.fromkeys(['SYBR'], 'non-saturating DNA binding dye')
        }
    
    
    # check the overall chemistry type based on the dyes present
    overall_chemistry_type = 'hydrolysis probe' if('FAM' or 'SUN' or 'VIC' in dyeIds)\
        else 'non-saturating DNA binding dye' 
        
    
    # add all these dyes into the rdml file
    for dye in dyeIds:
        fl.new_dye(dye)
        
    # Update the dye chemistry appropriately for each dye
    for dyeField in fl.dyes():
        dyeID = dyeField.__getitem__('id') # get the name of the dye
        dyeField.__setitem__('dyeChemistry', 
                             dictChemistry.get(dyeID, 'non-saturating DNA binding dye')) # assign the appropriate chemistry
    
    # Change sample_type to ntc when sample name has it 
    # for sample in fl.samples():
    #     sampleid = sample.__getitem__('id')
    #     if ('NTC' or 'ntc') in sampleid:
    #         sample.types()[0].__setitem__('type', 'ntc') # make sample type ntc if named so
    # Does not work :(          
    
    
    # Validate a file
    # print(fl.validate())
    # doesn't validate here but works in RDML-ninja
    # line 0 column 0 element 'dye' not expected
      
    
    # Run LinRegPCR
    
    # unpack the rdml file, 1st experiment, 1st run
    first_run_data = fl.experiments()[0].runs()[0]
    
    
        
    # amplification data - run linregpcr
    
    # default parameters for dye chemistry : If probes are in the run or not
    if(overall_chemistry_type == 'non-saturating DNA binding dye'):
        
        qpcr_result = first_run_data.\
            linRegPCR(updateRDML=True,
                      saveResultsCSV=True, 
                      timeRun=True, 
                      verbose=True)
            
        # run melt curve only if dye chemistry    
        melt_curve_result = first_run_data.\
            meltCurveAnalysis(updateRDML=True,
                              saveResultsCSV=True,
                              verbose=True)
        # melt curve not working with default parameters -- same as online   
    else:
        # Different parameters for probe based data : linregpcr
        # Include plateau phase
        
        qpcr_result = first_run_data.\
            linRegPCR(updateRDML=True,
                      saveResultsCSV=True, 
                      timeRun=True, 
                      verbose=True,
                      excludeNoPlateau = False)
    
        
    # save results to RDML
    fl.save(dirname + "results/" + filename_rdml + "fixpython.rdml")
    
    # save amplification analysis results as csv file
    with open('excel files/' + filename_rdml + "-linreg.csv", "w") as filehandle:
        filehandle.write(qpcr_result["resultsCSV"])
            
    # save melt curve analysis results as csv file: if a non empty result exists
    if('melt_curve_result' in locals() and len(melt_curve_result)):  
        with open('excel files/' + filename_rdml + "-melt.csv", "w") as filehandle:
            filehandle.write(melt_curve_result["resultsCSV"])

# usage :
# linregpcr_analyze('q22_fgfp-triplex-gradient_24-1-22')

