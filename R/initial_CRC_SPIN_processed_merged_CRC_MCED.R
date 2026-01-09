####################################################
# Function for inital_process_CRC_SPIN function.
####################################################
build_crc_data_revised <- function(crc_spin_scenario, CRC_SPIN_dir) {
  
  if (!crc_spin_scenario %in% c("control", "colonoscopy", "annual_fit",
                                "mced_annual", "mced_biennial", "mced_triennial",
                                "CRC_single_cancer_annual", "CRC_single_cancer_biennial", "CRC_single_cancer_triennial",
                                "MCED_annual_no_CRC", "MCED_biennial_no_CRC", "MCED_triennial_no_CRC")) {
    stop("Invalid crc_spin_scenario")
  }
  
  
  # CONTROL (no screening)
  if (crc_spin_scenario == "control") {
    control_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))
    
    crc_data_revised <- control_data %>%
      rename(start_age = age_at_index,
             other_cause_death_time = age_death_oc,
             clinical_diagnosis_time = age_clindet,
             cancer_death_time_no_screen = age_death_crc_noscreen) %>%
      mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
             clinical_diagnosis_stage = ifelse(stage_clindet %in% c(1, 2), "Early", "Late"),
             cancer_site = "Colorectal",
             # age at other-cause death, binned by 5-year intervals
             age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0,150,5))) %>%
      select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time,
             cancer_death_time_no_screen, clinical_diagnosis_stage, cancer_site, age_OC_death_cat)
    
    
    # =========================
    # Colonoscopy (scenarioID 2)
    # =========================
  } else if (crc_spin_scenario == "colonoscopy") {
    colonoscopy_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_2_person.rds"))
    
    crc_data_revised <- colonoscopy_data %>%
      rename( start_age = age_at_index,
              other_cause_death_time = age_death_oc,
              clinical_diagnosis_time = age_screendet,
              cancer_death_time_no_screen = age_death_crc_screendet ) %>%
      mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
             clinical_diagnosis_stage = ifelse(stage_screendet %in% c(1, 2), "Early", "Late"),
             cancer_site = "Colorectal",
             age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5)) ) %>%
      select(sex, p_id, start_age, other_cause_death_time,
             clinical_diagnosis_time, cancer_death_time_no_screen,
             clinical_diagnosis_stage, cancer_site, age_OC_death_cat)
    
    #==================================================
    # Scenario 3: Annual FIT screening (crc_spin_scenarioID = 3 )
    #   - CRC-SPIN FIT strategy
    #   - Screen diagnosis time = age_screendet
    #   - Cancer death time = age_death_crc_screendet
    #==================================================
  } else if (crc_spin_scenario == "annual_fit") {
    annual_fit_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_3_person.rds"))
    
    crc_data_revised <- annual_fit_data %>%
      rename(start_age = age_at_index,
             other_cause_death_time = age_death_oc,
             clinical_diagnosis_time = age_screendet,
             cancer_death_time_no_screen = age_death_crc_screendet) %>%
      mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
             clinical_diagnosis_stage = ifelse(stage_screendet %in% c(1, 2), "Early", "Late"),
             cancer_site = "Colorectal",
             age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
      select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time,
             cancer_death_time_no_screen, clinical_diagnosis_stage, cancer_site, age_OC_death_cat)
    
    
    #==================================================
    # Scenarios 4–6: MCED including CRC
    #   mced_annual / mced_biennial / mced_triennial
    #   - For each person, we keep:
    #       • control trajectory (no-screen CRC-SPIN)
    #       • MCED screen trajectory (MCED+CRC-SPIN)
    #   - Join on p_id so both paths are available
    #==================================================
    ############## Scenario 4: MCED + CRC, annual  ##############
  } else if (crc_spin_scenario %in% c("mced_annual", "mced_biennial", "mced_triennial")) {
    
    screen_id <- case_when(crc_spin_scenario == "mced_annual" ~ 4L,
                           crc_spin_scenario == "mced_biennial" ~ 5L,
                           TRUE ~ 6L)
    
    control_data <- readRDS( file.path(CRC_SPIN_dir,
                                       "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))
    
    screen_data <- readRDS( file.path(CRC_SPIN_dir,
                                      paste0("scr_parmID_0_blockID_0_scenarioID_",  screen_id, "_person.rds")))
    
    
    # control group (no screening) – used for comparison
    control_data_revised <- control_data %>%
      rename(start_age = age_at_index,
             other_cause_death_time = age_death_oc,
             clinical_diagnosis_time = age_clindet,
             cancer_death_time_no_screen = age_death_crc_noscreen) %>%
      mutate(sex = ifelse(fem==FALSE,"Male","Female"),
             clinical_diagnosis_stage = ifelse(stage_clindet %in% c(1,2),"Early","Late"),
             cancer_site="Colorectal",
             age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0, 150, 5))) %>%
      select(sex,p_id,start_age,other_cause_death_time,clinical_diagnosis_time,
             cancer_death_time_no_screen,clinical_diagnosis_stage,cancer_site,age_OC_death_cat)
    
    # MCED annual screen trajectory, then join in control path by p_id
    crc_data_revised <- screen_data %>%
      rename(start_age=age_at_index,
             other_cause_death_time=age_death_oc,
             screen_diagnosis_time=age_screendet,
             cancer_death_time_screen=age_death_crc_screendet) %>%
      mutate(sex=ifelse(fem==FALSE,"Male","Female"),
             screen_diagnosis_stage=ifelse(stage_screendet %in% c(1,2),"Early","Late"),
             cancer_site="Colorectal",
             age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0, 150, 5))) %>%
      select(sex,p_id,start_age,other_cause_death_time,screen_diagnosis_time,
             cancer_death_time_screen,screen_diagnosis_stage,cancer_site,age_OC_death_cat) %>%
      left_join(select(control_data_revised,p_id,clinical_diagnosis_time,cancer_death_time_no_screen), by="p_id")
    
    
    ## MCED scenarios with CRC only (single-cancer MCED) ########
    #==================================================
    # Scenarios 7–9: CRC_single_cancer (annual/biennial/triennial) 
    #   - CRC-only MCED tests (no multi-cancer detection)
    #   - Here CRC-SPIN scenario 4/5/6 already encode the CRC-only MCED strategy
    #==================================================
  } else if (crc_spin_scenario %in% c("CRC_single_cancer_annual", "CRC_single_cancer_biennial", "CRC_single_cancer_triennial")) {
    
    screen_id <- case_when(crc_spin_scenario == "CRC_single_cancer_annual" ~ 4L,
                           crc_spin_scenario == "CRC_single_cancer_biennial" ~ 5L,
                           TRUE ~ 6L)  
    
    CRC_single_cancer_data <- readRDS( file.path(CRC_SPIN_dir,paste0("scr_parmID_0_blockID_0_scenarioID_",  screen_id, "_person.rds")))
    
    crc_data_revised <- CRC_single_cancer_data %>%
      rename(start_age = age_at_index,
             other_cause_death_time = age_death_oc,
             clinical_diagnosis_time = age_screendet,
             cancer_death_time_no_screen = age_death_crc_screendet) %>%
      mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
             clinical_diagnosis_stage = ifelse(stage_screendet %in% c(1, 2), "Early", "Late"),
             cancer_site = "Colorectal",
             age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
      select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time,
             cancer_death_time_no_screen, clinical_diagnosis_stage, cancer_site, age_OC_death_cat)
    
    
    #==================================================
    # Scenarios 10–12: MCED_no_CRC (annual/biennial/triennial)
    #  Interpretation:
    #   • MCED is used to detect NON-CRC cancers.
    
    # Implementation:
    #   1) control_data_revised
    #         - clinical_diagnosis_time      = age_clindet
    #         - cancer_death_time_no_screen  = age_death_crc_noscreen
    #  2) MCED_no_CRC_data
    #             screen_diagnosis_time     = age_clindet
    #             cancer_death_time_screen  = age_death_crc_noscreen
    # 3) We left_join() the MCED path with the control path by p_id so 
    #          that both sets of CRC variables are available to downstream code.
    #==================================================
  } else if (crc_spin_scenario %in% c("MCED_annual_no_CRC",  "MCED_biennial_no_CRC", "MCED_triennial_no_CRC")) {
    
    screen_id <- case_when( crc_spin_scenario == "MCED_annual_no_CRC" ~ 4L,
                            crc_spin_scenario == "MCED_biennial_no_CRC" ~ 5L,
                            TRUE ~ 6L)  
    
    control_data <- readRDS( file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))
    MCED_annual_no_CRC_data <- readRDS( file.path(CRC_SPIN_dir, 
                                                  paste0("scr_parmID_0_blockID_0_scenarioID_",  screen_id, "_person.rds")))
    
    # control path 
    control_data_revised <- control_data %>%
      rename(start_age = age_at_index,
             other_cause_death_time = age_death_oc,
             clinical_diagnosis_time = age_clindet,
             cancer_death_time_no_screen = age_death_crc_noscreen) %>%
      mutate(sex = ifelse(fem==FALSE,"Male","Female"),
             clinical_diagnosis_stage = ifelse(stage_clindet %in% c(1,2),"Early","Late"),
             cancer_site="Colorectal",
             age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0, 150, 5))) %>%
      select(sex,p_id,start_age,other_cause_death_time,clinical_diagnosis_time,
             cancer_death_time_no_screen,clinical_diagnosis_stage,cancer_site,age_OC_death_cat)
    
    
    # MCED annual for non-CRC cancers
    crc_data_revised <- MCED_annual_no_CRC_data %>%
      rename(start_age=age_at_index,
             other_cause_death_time=age_death_oc,
             screen_diagnosis_time=age_clindet,
             cancer_death_time_screen=age_death_crc_noscreen) %>%
      mutate(sex=ifelse(fem==FALSE,"Male","Female"),
             screen_diagnosis_stage=ifelse(stage_clindet %in% c(1,2),"Early","Late"),
             cancer_site="Colorectal",
             age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0, 150, 5))) %>%
      select(sex,p_id,start_age,other_cause_death_time,screen_diagnosis_time,
             cancer_death_time_screen,screen_diagnosis_stage,cancer_site,age_OC_death_cat) %>%
      left_join(select(control_data_revised,p_id,clinical_diagnosis_time,cancer_death_time_no_screen), by="p_id")
    
  }
  
  return(crc_data_revised)
}


# ==================================
# Function to merged CRC-SPIN and MCED data
# ==================================
process_mced_crc_data <- function(scenario_name, mced_data_dir, crc_data_path) {
  
  # Pick which MCED results file to load based on scenario name
  if (scenario_name %in% c("mced_biennial", "MCED_biennial_no_CRC","CRC_single_cancer_biennial")) {
    mced_id <- 2L
  } else if (scenario_name %in% c("mced_triennial", "MCED_triennial_no_CRC","CRC_single_cancer_triennial")) {
    mced_id <- 3L
  } else {
    mced_id <- 1L
  }
  
  # Load MCED annual/biennial/triennial results
  load(file.path(mced_data_dir, paste0("mced_scenario_", mced_id, "_combined_results.Rdata")))
  if (!exists("combined_first_results")) stop("combined_first_results not found in MCED file")
  
  # Load CRC list (all 12 scenarios)
  load(crc_data_path)
  
  # Define which scenarios use clinical vs screen variables
  clinical_scenarios <- c("control", "colonoscopy", "annual_fit", "CRC_single_cancer_annual", "CRC_single_cancer_biennial","CRC_single_cancer_triennial")
  screen_scenarios <- c("mced_annual", "mced_biennial", "mced_triennial","MCED_annual_no_CRC", "MCED_biennial_no_CRC", "MCED_triennial_no_CRC")
  
  # Determine scenario type
  if (scenario_name %in% clinical_scenarios) {
    scenario_type <- "clinical"
  } else if (scenario_name %in% screen_scenarios) {
    scenario_type <- "screen"
  } else {
    stop("Invalid scenario_name: ", scenario_name, ". Valid options are: ", 
         paste(c(clinical_scenarios, screen_scenarios), collapse = ", "))
  }
  
  #============================================
  # Merge processed CRC-SPIN data with MCED data for scenarios without MCED screening for non-CRC cancers
  # Scenarios include control, colonoscopy, annual_fit, CRC_single_cancer
  # The processed CRC-SPIN specifies CRC cancer diagnoses with "clinical" labeled variables, because these are used
  # when converting to MCEDsim variables in absence of MCED screening for non-CRC cancers.  
  # Natural history assumption: fast-fast-optimistic
  #At the end of the script, we have two files: excess cancers that includes non-CRC cancers (storing time of clinical diagnosis)
  #                                             joined_data that consists of the merged CRC-SPIN and MCEDsim files
  #============================================  
  if (scenario_type == "clinical") {
    
    # Get CRC data from crc_list
    CRC_data <- crc_list[[scenario_name]]
    
    #Load processed CRC data
    CRC_data <- CRC_data %>% mutate(fem = ifelse(sex == "Female", TRUE, FALSE)) %>%
      select(c("fem","sex","p_id","start_age", "other_cause_death_time", "clinical_diagnosis_time",
               "cancer_death_time_no_screen", "clinical_diagnosis_stage", "cancer_site"))
    
    
    # Sort by sex and set new ids in both datasets
    CRC_data <- CRC_data %>% arrange(sex, .by_group = FALSE) %>% mutate(new_ID = seq(1, dim(CRC_data)[1]))
    combined_first_results <- combined_first_results %>% arrange(sex, .by_group = FALSE) %>% mutate(new_ID = seq(1, dim(combined_first_results)[1]))
    
    # Join datasets by new id
    # Note: variable names ending in .x are from MCEDsim .y are from CRC-spin
    joined_data <- combined_first_results %>% 
      left_join(select(CRC_data, c("fem","sex","p_id","new_ID","other_cause_death_time",
                                   "clinical_diagnosis_time", "cancer_death_time_no_screen",
                                   "clinical_diagnosis_stage","cancer_site")),  by = c("new_ID" = "new_ID"))
    
    
    #If sex allocations are not aligned, remove individuals with mismatched sex in joined data--may also have to remove these people from combined additional cancers
    #Also set other-cause death time to be CRC-SPIN OC death time
    joined_data <- joined_data %>%  filter(sex.x == sex.y) %>% mutate(other_cause_death_time = other_cause_death_time.y, sex = sex.x) 
    
    
    #If individual has CRC and no primary cancer before OC death, assign CRC to be the cancer site of primary cancer
    joined_data <- joined_data %>% mutate(cancer_site = ifelse(!is.na(clinical_diagnosis_time.y) & clinical_diagnosis_time.y < other_cause_death_time,
                                                               cancer_site.y, as.character(cancer_site.x)))
    
    #If individual has both CRC and one or more MCED cancers before OC death: move the MCED cancer to an additional cancer file
    joined_data <- joined_data %>% mutate(MCED_and_CRC = ifelse(cancer_site == "Colorectal" & onset_time < other_cause_death_time, TRUE, FALSE))
    
    
    # Extract individuals with both MCED and CRC cancers to a separate dataset
    # These MCED cancers will be treated as "additional cancers" rather than primary
    excess_cancers <- subset(joined_data, MCED_and_CRC == TRUE) %>% 
      mutate(clinical_diagnosis_time = clinical_diagnosis_time.x,
             cancer_death_time_no_screen = cancer_death_time_no_screen.x,
             clinical_diagnosis_stage = clinical_diagnosis_stage.x) %>%
      select(c(names(combined_first_results)))
    
    
    # Set variables in joined_data depending on if they have CRC or MCED as their primary cancer
    # If CRC is the primary cancer, use CRC values (.y), otherwise use MCED values (.x)
    joined_data <- joined_data %>% mutate(clinical_diagnosis_time = ifelse(cancer_site == "Colorectal", clinical_diagnosis_time.y, clinical_diagnosis_time.x),
                                          cancer_death_time_no_screen = ifelse(cancer_site == "Colorectal", cancer_death_time_no_screen.y, cancer_death_time_no_screen.x),
                                          clinical_diagnosis_stage = ifelse(cancer_site == "Colorectal", clinical_diagnosis_stage.y, clinical_diagnosis_stage.x)) %>% 
      select(-ends_with(".x"), -ends_with(".y"))
    
    #============================================    
    # Merge processed CRC-SPIN data with MCED data for scenarios with MCED screening for non-CRC cancers
    # Scenarios include mced, MCED_no_CRC
    # The processed CRC-SPIN specifies CRC cancer diagnoses with "screen" labeled variables, because these are used
    # when converting to MCEDsim variables in presence of MCED screening for non-CRC cancers.  
    # Natural history assumption: fast-fast-optimistic
    #At the end of the script, we have two files: excess cancers that includes non-CRC cancers (storing time of screened diagnosis)
    #                                             joined_data that consists of the merged CRC-SPIN and MCEDsim files
    #============================================    
  } else if (scenario_type == "screen") {
    
    # Get CRC data from crc_list  
    CRC_data <- crc_list[[scenario_name]]  
    
    # Load processed CRC data    
    CRC_data <- CRC_data %>% mutate(fem = ifelse(sex == "Female", TRUE, FALSE)) %>%
      select(c("fem","sex","p_id","start_age", "other_cause_death_time", "screen_diagnosis_time",
               "cancer_death_time_screen", "screen_diagnosis_stage", "cancer_site"))
    
    
    # sort by sex and set new ids in both datasets
    CRC_data <- CRC_data %>% arrange(sex, .by_group = FALSE) %>% mutate(new_ID = seq(1, dim(CRC_data)[1])) 
    combined_first_results <- combined_first_results %>% arrange(sex, .by_group = FALSE) %>% mutate(new_ID = seq(1, dim(combined_first_results)[1]))
    
    # Join datasets by new id
    # Note: variable names ending in .x are from MCEDsim .y are from CRC-spin
    joined_data <- combined_first_results %>% 
      left_join(select(CRC_data, c("fem","sex","p_id","new_ID","other_cause_death_time", 
                                   "screen_diagnosis_time", "cancer_death_time_screen",
                                   "screen_diagnosis_stage","cancer_site")), 
                by = c("new_ID" = "new_ID"))  
    
    #If sex allocations are not aligned, remove individuals with mismatched sex in joined data--may also have to remove these people from combined additional cancers
    #Also set other-cause death time to be CRC-SPIN OC death time
    joined_data <- joined_data %>% filter(sex.x == sex.y) %>% mutate(other_cause_death_time = other_cause_death_time.y, sex = sex.x)
    
    #If individual has CRC and no primary cancer before OC death, assign CRC to be the cancer site of primary cancer.
    joined_data <- joined_data %>% mutate(cancer_site = ifelse(!is.na(screen_diagnosis_time.y) & screen_diagnosis_time.y < other_cause_death_time,
                                                               cancer_site.y, as.character(cancer_site.x)))
    
    #If individual has both CRC and one or more MCED cancers before OC death: move the MCED cancer to an additional cancer file
    joined_data <- joined_data %>% mutate(MCED_and_CRC = ifelse(cancer_site == "Colorectal" & onset_time < other_cause_death_time, 
                                                                TRUE, FALSE)) 
    
    # Extract individuals with both MCED and CRC cancers to a separate dataset
    # These MCED cancers will be treated as "additional cancers" rather than primary
    excess_cancers <- subset(joined_data, MCED_and_CRC == TRUE) %>% mutate(screen_diagnosis_time = screen_diagnosis_time.x,
                                                                           cancer_death_time_screen = cancer_death_time_screen.x,
                                                                           screen_diagnosis_stage = screen_diagnosis_stage.x) %>%
      select(c(names(combined_first_results)))
    
    # Set variables in joined_data depending on if they have CRC or MCED as their primary cancer
    # If CRC is the primary cancer, use CRC values (.y), otherwise use MCED values (.x)
    joined_data <- joined_data %>% mutate(screen_diagnosis_time = ifelse(cancer_site == "Colorectal", screen_diagnosis_time.y, screen_diagnosis_time.x),
                                          cancer_death_time_screen = ifelse(cancer_site == "Colorectal",cancer_death_time_screen.y,cancer_death_time_screen.x),
                                          screen_diagnosis_stage = ifelse(cancer_site == "Colorectal", screen_diagnosis_stage.y, screen_diagnosis_stage.x)) %>%
      select(-ends_with(".x"), -ends_with(".y"))
  }
  
  return(list(joined_data = joined_data,
              excess_cancers = excess_cancers,
              scenario = scenario_name)
  )
  
}  






