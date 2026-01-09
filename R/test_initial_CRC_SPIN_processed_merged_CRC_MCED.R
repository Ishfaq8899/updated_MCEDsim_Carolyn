library(readr)
library(dplyr)

# ============================================================
# Test inital processed CRC SPIN data
# Build all 12 scenarios + save each one as its own .RData file
# ============================================================
data_dir <- "/Users/ahmadish/Desktop/OHSU_Projects/Ishfaq_folder/MCED_Colorectal_Project/MCEDsim_Carolyn/data"
CRC_SPIN_dir <- file.path(data_dir, "CRC_SPIN_data/CRCSPIN_Sample_Results_New_Survival")

# ---- Build scenarios ----
crc_control <- build_crc_data_revised("control", CRC_SPIN_dir)
crc_colon <- build_crc_data_revised("colonoscopy", CRC_SPIN_dir)
crc_fit <- build_crc_data_revised("annual_fit", CRC_SPIN_dir)
crc_mced_annual <- build_crc_data_revised("mced_annual", CRC_SPIN_dir)
crc_mced_biennial <- build_crc_data_revised("mced_biennial", CRC_SPIN_dir)
crc_mced_triennial <- build_crc_data_revised("mced_triennial", CRC_SPIN_dir)
crc_single_cancer_annual <- build_crc_data_revised("CRC_single_cancer_annual", CRC_SPIN_dir)
crc_single_cancer_biennial <- build_crc_data_revised("CRC_single_cancer_biennial", CRC_SPIN_dir)
crc_single_cancer_triennial <- build_crc_data_revised("CRC_single_cancer_triennial", CRC_SPIN_dir)
mced_annual_no_crc <- build_crc_data_revised("MCED_annual_no_CRC", CRC_SPIN_dir)
mced_biennial_no_crc <- build_crc_data_revised("MCED_biennial_no_CRC", CRC_SPIN_dir)
mced_triennial_no_crc <- build_crc_data_revised("MCED_triennial_no_CRC", CRC_SPIN_dir)


# Build the 12 scenarios
crc_list <- list(control = crc_control,
                 colonoscopy = crc_colon,
                 annual_fit = crc_fit,
                 mced_annual = crc_mced_annual,
                 mced_biennial = crc_mced_biennial,
                 mced_triennial = crc_mced_triennial,
                 CRC_single_cancer_annual = crc_single_cancer_annual,
                 CRC_single_cancer_biennial = crc_single_cancer_biennial,
                 CRC_single_cancer_triennial = crc_single_cancer_triennial,
                 MCED_annual_no_CRC = mced_annual_no_crc,
                 MCED_biennial_no_CRC = mced_biennial_no_crc,
                 MCED_triennial_no_CRC = mced_triennial_no_crc)


# Output folder
out_dir <- file.path(data_dir, "initial_process_CRC_SPIN_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save(crc_list, file = file.path(out_dir, "crc_spin_initial_process_data.RData"))


########################################
# Test mereged CRC-SPIN MCED data
########################################

all_scenarios <- c("control", "colonoscopy", "annual_fit",
                   "CRC_single_cancer_annual", "CRC_single_cancer_biennial", "CRC_single_cancer_triennial",
                   "mced_annual", "mced_biennial", "mced_triennial",
                   "MCED_annual_no_CRC", "MCED_biennial_no_CRC", "MCED_triennial_no_CRC")

all_results <- list()

mced_data_dir <- "/Users/ahmadish/Desktop/OHSU_Projects/Ishfaq_folder/MCED_Colorectal_Project/MCEDsim_Carolyn/data/MCED_data"
crc_data_path <- "/Users/ahmadish/Desktop/OHSU_Projects/Ishfaq_folder/MCED_Colorectal_Project/MCEDsim_Carolyn/data/initial_process_CRC_SPIN_data/crc_spin_initial_process_data.RData"


for (scenario in all_scenarios){
  all_results[[scenario]] <- process_mced_crc_data(scenario_name = scenario,
                                                   mced_data_dir = mced_data_dir,
                                                   crc_data_path = crc_data_path
  )
}

output_dir <- "/Users/ahmadish/Desktop/OHSU_Projects/Ishfaq_folder/MCED_Colorectal_Project/MCEDsim_Carolyn/data"

for (scenario in all_scenarios) {
  joined_data <- all_results[[scenario]]$joined_data
  excess_cancers <- all_results[[scenario]]$excess_cancers

  save(joined_data, excess_cancers,
       file = file.path(output_dir, paste0(scenario, "_processed_data.RData")))
}









