if (!exists("iph")) source("00_load_data.R")
source("./00_Data_Processing/remove_low_prevalence.R")
source("./01_Regression/lasso.R")

################################################################################
## Settings
################################################################################

# The !exists check makes it easier to source this script from another location
# control the settings

# Create an rds file with the vector of variables which meet the prevalence cutoff?
if (!exists("save_min_prev_vars_flag")) save_min_prev_vars_flag <- TRUE

# Create rds files for the model formulas?
if (!exists("save_formula_flag")) save_formula_flag <- TRUE

# Directory to save the model formulas to
if (!exists("formula_save_dir")) formula_save_dir <- "./Work_Products/Model_Formulas/"

# Create a csv with with a row for every column in the data set and columns 
# indicating whether it was considered in the model, whether it has low prevalence, 
# and whether it is included in the model 
if (!exists("save_full_variable_summary_flag")) save_full_variable_summary_flag <- FALSE

################################################################################
## Variables for consideration
################################################################################

demo_vars <-  c("age",
                "sexBinary",
                "dateInterval",
                paste0("race___", 0:4),
                "ethnicityRecoded")

state_vars <- c("stateCatMixWOther")

country_vars <- c("countryGrouped", 
                  "CountryIncomeGrp")

country_indicators <- iph %>% select(starts_with("countryInd")) %>% colnames


patient_char_vars <- c("diagnosisRecoded",
                       "pgaUnknownRecoded",
                       paste0("comorbid___", 1:11),
                       "heightCombinedCm",
                       "weightBMIAbv60Rm",
                       "bmiAbv60Rm",
                       "ComorbidCount")

med_subcat_vars <- c( "asa___1",
                      "asa___2",
                      "Bud_flg", 
                       "Other_CS_flg",
                       "TNF_flg",
                       "INT_flg",
                       "IL_flg",
                       "ImmodMet",
                       "Immod6mpOrAza")

med_dose_vars <- c("preddose_daily",
                   "azadose",
                   "mcpdose")

med_dose_weight_vars <- c("preddose_dailyByWt",
                          "mcpdoseByWt", 
                          "azadoseByWt")

med_cat_vars <-  c(paste0("med_cat___", 1:7))

# Comorbidity interactions were manually cut off September 23 2020
# At least 30 with both comorbidities

interaction_terms_wo_country <- c("age:sexBinary",
                       "age:ComorbidCount",
                       "age:bmiAbv60Rm",
                       "comorbid___11:ComorbidCount",
                       "comorbid___1:comorbid___2",
                       "comorbid___1:comorbid___6",
                       "comorbid___2:comorbid___6",
                       "comorbid___6:comorbid___9",
                       paste0("med_cat___1:med_cat___", 2:4),
                       paste0("med_cat___2:med_cat___", 3:4),
                       "med_cat___3:med_cat___4",
                       "Other_CS_flg:age",
                       "Other_CS_flg:sexBinary",
                       "Other_CS_flg:bmiAbv60Rm",
                       "I(age^2)",
                       "I(ComorbidCount^2)",
                       "I(bmiAbv60Rm^2)",
                       "I(preddose_daily^2)",
                       "I(azadose^2)",
                       "I(mcpdose^2)")

interaction_terms_country <- c("dateInterval:countryGrouped")

interaction_terms_country_ind <- paste0("dateInterval:", country_indicators)

interaction_terms <- c(interaction_terms_wo_country, interaction_terms_country)

candidate_vars <- c(demo_vars, state_vars, country_vars,
                    patient_char_vars, med_subcat_vars, 
                    med_dose_vars,
                    med_dose_weight_vars, med_cat_vars)


################################################################################
## Determine what's included
################################################################################
count_cutoff <- 25
candidate_vars_with_min_prev <- DetermineVariablesForConsideration(ibd.dat = iph, 
                                                                   candidate.vars = candidate_vars, 
                                                                   count.cutoff = count_cutoff)

if (save_min_prev_vars_flag == TRUE) saveRDS(candidate_vars_with_min_prev, 
                                             paste0(formula_save_dir, "vars_w_min_prev.rds"))

low_prevalence_vars <- candidate_vars[!candidate_vars %in% candidate_vars_with_min_prev]

################################################################################
## Create formulas and save them
################################################################################

hosp_formula <- ConstructLassoFormula("HospOrDeath", 
                                      covariate.names = candidate_vars_with_min_prev,
                                      interaction.terms = interaction_terms)

icu_formula <- ConstructLassoFormula("IcuVentDeath", 
                                     covariate.names = candidate_vars_with_min_prev,
                                     interaction.terms = interaction_terms)

death_formula <- ConstructLassoFormula("death", 
                                       covariate.names = candidate_vars_with_min_prev,
                                       interaction.terms = interaction_terms)

ord_formula <- ConstructLassoFormula("HospIcuDthOrd", 
                                       covariate.names = candidate_vars_with_min_prev,
                                       interaction.terms = interaction_terms)

if (save_formula_flag == TRUE) {
  saveRDS(hosp_formula, paste0(formula_save_dir, "hosp_formula.rds"))
  saveRDS(icu_formula,  paste0(formula_save_dir, "icu_formula.rds"))
  saveRDS(death_formula, paste0(formula_save_dir, "death_formula.rds"))
  saveRDS(ord_formula, paste0(formula_save_dir, "ord_formula.rds"))
}

################################################################################
## Create variable inclusion summary
################################################################################

iph_raw <- ReadIPHData(path = paste0(data_path, analytic_filename), 
                       col_spec_df = colspec_data,
                       clean = FALSE)

variable_summary <- tibble(Variable = colnames(iph)) %>% 
  mutate(Considered = factor(Variable %in% candidate_vars,
                             levels = c(TRUE, FALSE),
                             labels = c("Yes", "No")),
         Source = if_else(Variable %in% colnames(iph_raw), "SECURE-IBD", "BIOS"),
         LowPrevalence = factor(Variable %in% low_prevalence_vars,
                                levels = c(TRUE, FALSE),
                                labels = c("Yes", "No")),
         IncludedInModel = factor(Variable %in% candidate_vars_with_min_prev,
                                  levels = c(TRUE, FALSE),
                                  labels = c("Yes", "No")))

if (save_full_variable_summary_flag == TRUE) {
  write_csv(variable_summary, 
            path = paste0(formula_save_dir, "variable_inclusion_summary.csv"))
}