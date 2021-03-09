# qPCR processing: Calculate copy number from Cq and attach sample labels from template table 


process_qpcr <- function(flnm = flnm.here, std_override = NULL, baylor_wells = 'none')
  # enter the file name, standard curve mentioned within filename is used unless override is provided
{ # baylor_wells = # choose : none, '.*' for all, '[A-H]([1-9]$|10)' for columns 1 - 10; '.*(?!3).$' for everything except 3rd column etc.  
  #(will append /baylor to target name; Ad hoc - marking the samples from baylor)
  
  # Data input ----
  
  # Preperation steps
  flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
  
  # dynamic update of standard curve parameters : default parameters in the inputs_for_analysis.R file
  std_to_retrieve <- if_else(is.null(std_override), str_match(flnm, 'Std[:alnum:]*'), std_override)
  
  std_par_update <- read_sheet(sheeturls$data_dump, sheet = 'Standard curves', range = 'A:G', col_types = 'ccnnnnn') %>% 
    filter(str_detect(ID, std_to_retrieve ))
  
  # substitute the new std curve parameters in the old matrix
  std_par %<>% filter(!str_detect(Target,
                                  std_par_update$Target %>% str_c(collapse = "|"))) %>% 
    bind_rows(std_par_update)
  
  # error catching for repeated standard curve names for the same target
  if(std_par_update %>% unique() %>% group_by(Target) %>% summarize(count = n()) %>% pull(count) %>% {. > 1} %>% any())
  {
    stop( str_c('Duplicate entries for the same Target found for standard curve: ', std_to_retrieve, '. \n
  check the Standard curves sheet here : https://docs.google.com/spreadsheets/d/1ouk-kCJHERRhOMNP07lXfiC3aGB4wtWXpnYf5-b2CI4/edit#gid=1980064476'))
  }
  
  # Read in qPCR data and labels from plate template
  fl <- readqpcr(flpath) # read excel file exported by Quantstudio
  plate_template <- get_template_for(flnm, sheeturls$templates)
  
  sample_order = columnwise_index(fl) # this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order
  
  # Load desired qPCR result sheet and columns
  bring_results <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, starts_with('Tm'),`Target Name`) %>% rename(Target = `Target Name`) %>%  .[sample_order,] %>%  
    
    # select only the results used for plotting, calculations etc. and arrange them according to sample order
    select(-`Sample Name`) %>% right_join(plate_template, by = 'Well Position') %>%  # Incorporate samples names from the google sheet by matching well position
    mutate_at('Target', ~str_replace(., 'BSRV', 'BRSV')) %>% 
    filter(!is.na(Target))
  
  # Remove unneccesary data
  rm(fl)  # remove old data for sparsity
  
  # Data polishing ----
  
  
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  
  # isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
  polished_results <- bring_results %>% separate(`Sample_name`,c(NA, 'Sample_name'),'-') %>% separate(`Sample_name`,c('Sample_name','Tube ID'),'_') %>% 
    mutate(`Tube ID` = if_else(`Sample_name` == 'NTC', '0', `Tube ID`)) %>% 
    separate(`Tube ID`, c('assay_variable', 'biological_replicates'), remove = F) %>%  # Separate out biological replicates 
    unite('Tube ID', c(assay_variable, biological_replicates), sep = '.', remove = F, na.rm = T) %>% # remaking Tube ID - removes spaces after 'dot'
    arrange(assay_variable, biological_replicates) %>% mutate_if(is.character,as_factor) # Factorise the sample name and rearrange in column order of appearance on the plate (for plotting)
  
  # select samples to plot (or to exclude write a similar command)
  results_relevant <- polished_results %>% filter(str_detect(`Sample_name`, paste('^', plot_select_facet, sep = ''))) %>%  
    
    # Include only desired facets : str_detect will find for regular expression; ^x => starting with x
    filter(!str_detect(`Sample_name`, plot_exclude_facet)) %>%  # exclude unwanted facets (sample_name) 
    filter(!str_detect(assay_variable, plot_exclude_assay_variable)) %>%  # excluding unwanted x axis variables from assay_variable
    
    # Adding tag to target for baylor smaples
    { if(!str_detect(baylor_wells, 'none|None')) { 
      mutate_at(., 'Target', as.character) %>% 
        mutate_cond(str_detect(`Well Position`, baylor_wells), Target = str_c(Target, '/Baylor'))
    } else .
    }
  
  # Computation ----
  
  
  # Computing copy number from standard curve linear fit information
  results_abs <- results_relevant %>% group_by(Target) %>% do(., absolute_backcalc(., std_par)) %>%  # iteratively calculates copy #'s from standard curve parameters of each Target
    mutate(`Copy #` = `Copy #`/template_volume_qpcr) # Normalizing copy number per micro litre of template in the reaction
  
  # Finding mean and standard deviation within replicates (both technical and biological)
  
  summary_results <- results_abs %>%  group_by(`Sample_name`, Target, assay_variable) %>% summarise_at(vars(`Copy #`), lst(mean, sd), na.rm = T) # find mean and SD of individual copy #s for each replicate
  results_abs$`Copy #` %<>% replace_na(0) # make unamplified values 0 for plotting
  
  plt <- results_abs %>% ggplot(aes(x = `Tube ID`, y = `Copy #`, color = Target)) + ylab('Copies/ul RNA extract') +    # Specify the plotting variables 
    geom_point(size = 2) + facet_grid(~`Sample_name`, scales = 'free_x', space = 'free_x') + # plot points and facetting
    ggtitle(flnm) + xlab(plot_assay_variable)
  plt.formatted <- plt %>% format_classic(.) %>% format_logscale_y() # formatting plot, axes labels, title and logcale plotting
  
  print(plt.formatted)
  
  # Data output ----
  
  # Check for pre-existing file and write. Ask for overwrite permission
  check_ok_and_write(results_abs, sheeturls$data_dump, flnm)
  
  # write_sheet(results_abs, sheeturls$data_dump, sheet = flnm) # save results to a google sheet
  # ggsave('qPCR analysis/', WW1_Baylor-bovine_pilot.png', plot = plt.formatted, width = 8, height = 4)
  
  
  # Saving vaccine data into Vaccines sheet in data dump: For easy book keeping
  vaccine_data <- results_abs %>% filter(str_detect(Sample_name, 'Vaccine')) %>%
    mutate('Prepared on' = '',
           Week = str_extract(flnm, '[:digit:]{3,4}') %>% unlist() %>% str_c(collapse = ', '),
           Vaccine_ID = assay_variable, 
           .before = 1) %>% 
    mutate(Run_ID = str_extract(flnm, 'WW[:digit:]*'))
  
  # Add to existing sheet
  if(vaccine_data %>% plyr::empty() %>% {!.}) sheet_append(sheeturls$data_dump, vaccine_data, 'Vaccines')
  
  # Mean of vaccine data
  vaccine_data.mean <- vaccine_data %>% ungroup() %>% 
    select(1:3, Target, `Copy #`, Run_ID) %>% group_by(across(-`Copy #`)) %>% 
    summarise(across(`Copy #`, list(Mean_qPCR = mean, SD_qPCR = sd), na.rm = T), .groups = 'keep') %>% 
    mutate('[Stock conc.] copies/ul' = `Copy #_Mean_qPCR` * 50/20,
           'Estimated factor' = '',
           Comments = '',
           'Conc normalized to estimated factor' = '') %>% 
    relocate(Run_ID, .after = last_col()) %>% 
    mutate('x' = '', .before = 1)
  
  # Add to existing sheet in Vaccine_summary
  if(vaccine_data.mean %>% plyr::empty() %>% {!.}) sheet_append(sheeturls$data_dump, vaccine_data.mean, 'Vaccine_summary')
}