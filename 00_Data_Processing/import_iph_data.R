################################################################################
## Functions for reading in the registry data and applying data cleaning steps
##
################################################################################

#' A function for reading in the analytic data set
#' @param path 
#' @param col_spec An optional specification of the column types
#' @param clean A logical variable specifying if the cleaning steps should be applied. Defaults to TRUE
#' @return A tibble containing the cleaned study data
#' 
ReadIPHData <- function(path, col_spec_df = NULL, clean = TRUE, 
                        state.map.path = NULL,
                        hardcode.dir = NULL,
                        country.inc.dir = NULL) {
  if (is_logical(clean) == FALSE) stop("clean must be a logical variable")
  if (clean == TRUE & is.null(state.map.path)) stop("Must specify the path to the state region mapping")
  if (clean == TRUE & is.null(hardcode.dir)) stop("Must specify the path to the hard code directory")
  if (clean == TRUE & is.null(country.inc.dir)) stop("Must specify the path to the directory containing country income information")
  
  if (is.null(col_spec_df) == FALSE) {
    # Noticed that the order of some variables has changed from export to export
    unspec_read <- suppressWarnings(read_csv(path)) 
    cur_var_data <- tibble(Colname = names(unspec_read))
    merged_col_spec <- left_join(cur_var_data, col_spec_df, by = "Colname")
    
    if (all(merged_col_spec$Colname == names(unspec_read)) == FALSE) {
      stop("Column specifications not in the same order as the column names for this data")
    }
    
    if (any(is.na(merged_col_spec$ChrSpec))) {
      missing_col_spec_vars <- merged_col_spec %>% filter(is.na(ChrSpec)) %>% pull(Colname)
      warn_msg_preamble <- "Column specification does not match columns for this data - missing column specs for: \n"
      warn_msg_vars <- paste(missing_col_spec_vars, collapse = "\n") 
      unspec_warn_msg <- paste(warn_msg_preamble, warn_msg_vars, sep = "")
      warning(unspec_warn_msg)
      merged_col_spec$ChrSpec <- replace_na(merged_col_spec$ChrSpec, "?")
    } 
    
    col_spec <- merged_col_spec %>% pull(ChrSpec) %>% paste(., collapse = "")
  }
  
  if (is.null(col_spec_df) == TRUE) col_spec <- NULL

    iph <- read_csv(path, col_types = col_spec)
  
  if(clean == TRUE){
    iph <- iph %>% 
    .ApplyHardcodeChangesCountryVars(., hardcode.dir = hardcode.dir) %>% 
    .CreateAdditionalContinentVars(.) %>%
    .RecodeDemoVars(.) %>%
    .CheckAges(.) %>% 
    .RecodeClinicalVars(.) %>%
    .RecodeOutcomeVars(.)%>% 
    .RecodeCovidVars(.) %>% 
    .DeriveClinicalVars(.) %>% 
    .RecodeDrugVars(.) %>% 
    .DeriveDrugVars(.) %>% 
    .DeriveComorbidVars(.) %>% 
    .DeriveOutcomeVars(.) %>% 
    .CreateIncomeClassVars(., country.inc.dir = country.inc.dir) %>%
    .CreateAdditionalStateVars(., state.map.path = state.map.path) %>% 
    .DeriveCountryBinIndicators(.)
    
  }
  
  return(iph)
}


#' Apply hard code changes to missing country and state variables
#' Changes were identified by Googling the reporting physician and practice
.ApplyHardcodeChangesCountryVars <- function(in.dat, hardcode.dir){
  country_path = paste0(hardcode.dir, "missing_country_hard_codes.csv")
  state_path = paste0(hardcode.dir, "missing_state_hard_codes.csv")
  
  missing_country_hard_codes <- read_csv(country_path)
  missing_state_hard_codes <- read_csv(state_path)
  
  updated.dat <- rows_update(x = in.dat, 
                             y = missing_country_hard_codes, by = "numeric_studyid") %>%
    rows_update(., y = missing_state_hard_codes, by = "numeric_studyid")
  
  return(updated.dat)
}


#' A function for recoding the demographic variables. Examples include replacing "Missing" in country with NA,
#' creating a factor for ethnicity
#' @param path 
#' @param uknown.as.missing 
#' @return A tibble where the demographic variables have been recoded

.RecodeDemoVars <- function(iph_dat) {
  recoded_dat <- iph_dat %>%  
    mutate_at(., .vars = "country", ~na_if(., "Missing")) %>%
    mutate_at(., .vars = "DATE", ~lubridate::ymd(.)) %>%
    mutate(dateInterval = lubridate::interval(min(DATE), DATE)/lubridate::days(1)) %>% 
    # some country have very low frequency and if we seperate into test/training set, 
    # the may not show up in training set at all
    mutate(countryGrouped = as.factor(ifelse(country %in% 
                                               as.vector((as.data.frame(table(iph_dat$country)) %>% filter(Freq > 30))$Var1), 
                                            as.vector(country), 
                                            "others")
                                      )
           ) %>%  #
    mutate(countryImportant = as.factor(ifelse(is.na(country_cat) == TRUE, "others",as.vector(country)))) %>% #
    mutate_at(., .vars = "ethnicity", ~factor(.,levels = 1:3, labels = c("Hispanic", "Not Hispanic", "Unknown"))) %>% 
    mutate_at(., .vars = "sex", ~na_if(., "")) %>% 
    mutate(BMIGeq30 = ifelse(bmi_final >= 30, TRUE, FALSE),
           ethnicityRecoded = ethnicity) %>%
    mutate_at(., .vars = "ethnicityRecoded", ~droplevels(na_if(., "Unknown"))) %>% 
    group_by(numeric_studyid) %>% 
    # If the height is between .9 and 2.25 assume that they entered as meters rather than centimeters
    # If the height is between 36 to 90 assume that they entered as inches rather than centimeters
    # If the height is between 90 and 230 assume it was entered correctly
    # Otherwise assume it was incorrectly entered in a non-obvious way and mark as missing
    # For weight, assume <27 is an error unless they are 12 or younger
    # assume >= 250 and height units inches means they entered weight in pounds instead of kilograms
    mutate(heightCmCleaned = case_when(
      heightcm >= 90 & heightcm <= 230 ~ heightcm,
      heightcm >= 0.9 & heightcm <= 2.25 ~ 100*heightcm,
      heightcm >= 36 & heightcm <= 90 ~ 2.54*heightcm),
      weightCleaned = case_when(
        weight >= 27 & weight <= 250 ~ weight,
        weight <= 27 & age <= 12 ~ weight,
        weight >= 250 & heightunits == 2 & sex == "Male" ~ weight/2.205,
        weight >= 200 & heightunits == 2 & sex == "Female" ~ weight/2.205,
        weight >= 250 & heightunits == 1 ~ weight
      )) %>% 
    mutate(
      heightCombinedCm = case_when(
      is.na(heightCmCleaned) == FALSE ~ heightCmCleaned,
      is.na(heightintocm) == FALSE ~ heightintocm),
      bmiDerived = weightCleaned/((heightCombinedCm/100)^2),
      sexBinary = droplevels(as.factor(na_if(sex, "Other")))) %>% 
    ungroup %>% 
    mutate(bmiAbv60Rm = ifelse(bmiDerived > 60, NA, bmiDerived),
           weightBMIAbv60Rm = ifelse(bmiDerived > 60, NA, weightCleaned)) %>% 
    group_by(countryGrouped) %>% 
    mutate(dateCountryInterval = lubridate::interval(min(DATE), DATE)/lubridate::days(1)) %>% 
    ungroup
    
  return(recoded_dat)
}

#' A function for creating additional state variables
.CreateAdditionalStateVars <- function(iph_dat, state.map.path) {
  recoded_dat <- iph_dat %>%  
    mutate(stateWithNonUSLevel = factor(case_when(
      as.character(countryGrouped) != "United States" ~ "Other country",
      !is.na(state) ~ as.character(state))))
  
  state_mapping <- read_csv(state.map.path, col_types = "ffff")
  
  recoded_dat <- left_join(recoded_dat, state_mapping, by = "state")
  
  recoded_dat$stateCatMixWOther <- if_else(recoded_dat$country != "United States", 
                                     "Other country", 
                                    as.character(recoded_dat$stateCatMix)) %>% 
    factor
  
  recoded_dat$stateCensusDivisionWOther <- if_else(recoded_dat$country != "United States", 
                                           "Other country", 
                                           as.character(recoded_dat$stateCensusDivision)) %>% 
    factor
  
  recoded_dat$stateRegionWOther <- if_else(recoded_dat$country != "United States", 
                                                   "Other country", 
                                                   as.character(recoded_dat$stateRegion)) %>% 
    factor
  
  return(recoded_dat)  
}


#' A function for creating additional continent variables
.CreateAdditionalContinentVars <- function(iph_dat) {
  recoded_dat <- iph_dat
  recoded_dat$country = as.factor(iconv(recoded_dat$country, to='ASCII//TRANSLIT'))
  recoded_dat$country = plyr::revalue(recoded_dat$country, 
                                c("Bolivia, Plurinational State of"="Bolivia"))
  
  recoded_dat$continent = countrycode::countrycode(sourcevar = recoded_dat$country,
                                      origin = "country.name",
                                      destination = "continent")
  
  recoded_dat2 = iph_dat
  recoded_dat2$continent = recoded_dat$continent
  
  return(recoded_dat2)
}

#' A function for creating additional country classify by income variables
.CreateIncomeClassVars <- function(iph_dat, country.inc.dir) {
  incomeClassDat <- read.csv(paste0(country.inc.dir, "incomeclass.csv"),skip=1, header = FALSE)
  colnames(incomeClassDat) <- c("country","income_class")
  
  recoded_dat <- iph_dat
  # change special characters to normal
  recoded_dat$country = as.factor(iconv(recoded_dat$country, to='ASCII//TRANSLIT'))
  # change "XXX,XXX" to "XXX"  trim the name after a comma
  recoded_dat$country = as.factor(gsub("(.*),.*", "\\1", recoded_dat$country))
  # join two datasets
  recoded_dat2 = left_join(recoded_dat,incomeClassDat, by = "country") %>% 
    mutate(income_class = factor(income_class),
           CountryIncomeGrp = fct_collapse(income_class, "Low or lower-middle" = c("low", "lower_middle")))
  return(recoded_dat2)
}

#' A function for recoding the clinical variables. Examples include converting treatment frequency variables from 
#' factors, which represent weeks, to continuous, which represent days
#' @param path
#' @return A tibble where the treatment frequency variables have been recoded
.RecodeClinicalVars <- function(iph_dat) {
  
  # Frequency variables in the dataset
  freq_vars <- c("inflixfreq", "adadosefreq", "vedofreq", "ukfreq", "cpfreq", 
                 "gofreq", "natfreq", "otherfreq", "sulfafreq", "mesafreq", 
                 "mtxfreq", "csfreq", "tacfreq", "azafreq", "mcpfreq", 
                 "budfreq", "predfreq", "ivsfreq", "tofafreq", "othermedfreq")
  
  # Converts treatment frquency factors to number of days between treatment
  recode_freq <- function(x) {
    x_new <- case_when(
      x == 1 ~ 1, x == 2 ~ 3.5, x == 3 ~ 7, x == 4 ~ 14,
      x == 5 ~ 21, x == 6 ~ 28, x == 7 ~ 35, x == 8 ~ 42,
      x == 9 ~ 49, x == 10 ~ 56, x == 11 ~ 63, x == 12 ~ 70
    )
    return(x_new)
  }
  
  recoded_dat <- iph_dat %>%
    mutate_at(.vars = freq_vars,
              .funs = recode_freq) %>% 
    mutate(pgaUnknownRecoded =  pga,
           diagnosisRecoded = diagnosis) %>%
    mutate_at(., .vars = "pgaUnknownRecoded", ~droplevels(na_if(., "Unknown"))) %>%
    mutate_at(., .vars = "pgaUnknownRecoded", 
              ~ordered(., 
                       levels = c("Remission", "Mild", "Moderate", "Severe"))) %>% 
    mutate_at(., .vars = "diagnosisRecoded", ~droplevels(na_if(., "Inflammatory bowel disease unspecified")))
  
    return(recoded_dat)
}


#' A function for recoding the outcome variables. Examples include replacing "Unknown" for death with NA
#' @param path
#' @return A tibble where the outcome variables have been recoded

.RecodeOutcomeVars <- function(iph_dat) {

  int_coded_outcomes_with_unknown <- c("death", "gisx", "clot", "emerg")
  chr_coded_outcomes <- c("icu", "hosp", "vent")
  # For modeling, unknown outcomes should be considered missing (NA)
  # For matching pilot exploratory tables style, unknown should be kept as a category
  recoded_dat <- iph_dat %>%  
    mutate(deathCat = factor(death, levels = 1:3, labels = c("Died", "Lived", "Unknown")),
           emergCat = factor(emerg, levels = 1:3, labels = c("Evaluated in ER", "Not evaluated in ER", "Unknown")),
           clotCat = factor(clot, levels = 1:3, labels = c("Developed any thrombotic complications", 
                                                           "Did not develop any thrombotic complications", "Unknown")),
           gisxCat = factor(gisx, levels = 1:3, labels = c("Developed new gastrointestinal symptoms at infection", 
                                                           "Did not develop new gastrointestinal symptoms at infection", 
                                                           "Unknown")),
           ventCat = factor(vent),
           hospCat = factor(hosp),
           icuCat = factor(icu)) %>% 
    mutate_at(., .vars = int_coded_outcomes_with_unknown, ~na_if(., 3)) %>% 
    mutate_at(., .vars = int_coded_outcomes_with_unknown, ~ifelse(. == 1, TRUE, FALSE)) %>% 
    mutate_at(., .vars = chr_coded_outcomes, ~na_if(., "Unknown")) %>% 
    mutate_at(., .vars = chr_coded_outcomes, ~ifelse(. == "Yes", TRUE, FALSE)) %>% 
    mutate(losAmongHospitalized = los,
           los = ifelse(hosp == FALSE, 0, los),
           vent = ifelse(hosp == FALSE, 0, vent),
           icu = ifelse(hosp == FALSE, 0, icu),
           ventCat = ifelse(hosp == FALSE, "No", ventCat),
           icuCat = ifelse(hosp == FALSE, "No", icuCat)) %>% 
    mutate(HospOrDeath = hosp | death,
           IcuVentDeath = icu | vent | death,
           IcuOrVent = icu | vent)
    
  
  return(recoded_dat)
}

#' Apply data cleaning steps to the drug dosage variables in the model. 
#' Checks that the frequency matches an actual dosage frequency that would be prescribed,
#' and that dosages are in a prescribing range. Also recodes dosages for those not 
#' taking a medication to zero instead of missing because the dosage is known and zero,
#' not missing. 
.RecodeDrugVars <- function(iph_dat){
  ## preddose_daily
  # This is already a cleaned variable from our collaborators
  # Those who don't take steroids should have a dosage of zero, not missing
  # Those who report taking corticosteroids but don't report a dosage are true missings
  
  ## mcpdose
  # Frequency should be daily only
  
  ## azadose
  # Dosage should be between 25 and 250
  # Frequency should be daily only
  recoded_dat <- iph_dat %>% 
    group_by(numeric_studyid) %>% 
    mutate(
      preddose_daily = case_when(
        !is.na(preddose_daily) ~ preddose_daily,
        is.na(preddose_daily) & `med_cat___4` == FALSE ~ 0),
      mcpdose = case_when(
        mcpfreq != 1 ~ NA_real_,
        !is.na(mcpdose) ~ mcpdose,
        is.na(mcpdose) & `immod___5` == FALSE ~ 0),
      azadose = case_when(
        azadose > 250 ~ NA_real_,
        azadose < 25 ~ NA_real_,
        azafreq != 1 ~ NA_real_,
        !is.na(azadose) ~ azadose,
        is.na(azadose) & `immod___4` == FALSE ~ 0)) %>% 
    ungroup
  
  
  return(recoded_dat)
}

#' A function for recoding the COVID-19 variables. Examples include replacing "Unknown" for death with NA
#' @param path
#' @return A tibble where the outcome variables have been recoded
.RecodeCovidVars <- function(iph_dat){
  recoded_dat <- iph_dat %>%  
    mutate_at(., .vars = "covidpcr", 
    ~factor(., levels = 1:3, labels = c("Nasopharyngeal PCR", "Saliva PCR", 
    "Suspected symptom-based diagnosis with confirmatory antibody serology"))) %>% 
    mutate_at(., "sxres", ~factor(., 1:4, c("Yes", "No", "Unknown", "Asymptomatic")))
    
  return(recoded_dat)
}

#' A function for creating derived clinical variables. Includes:
#' a binary remission or not indicator for disease activity
.DeriveClinicalVars <- function(iph_dat){
  to_return <- iph_dat %>% 
    mutate(Remission = ifelse(pga == "Remission", TRUE, FALSE),
           Crohns = ifelse(diagnosis == "Crohn's disease", TRUE, FALSE))
    
  return(to_return)
}

#' Does not derive an indicator for the United States
.DeriveCountryBinIndicators <- function(iph_dat){
  # The prevalence cut-off should cut off those countries with fewer than 30 people in the September data.
  country_count_cutoff <- 30
  
  most_reported_countries <- iph_dat %>% 
    filter(DATE <= "2020-09-15") %>% 
    dplyr::count(country, sort = TRUE) %>% 
    filter(n >= country_count_cutoff) %>% 
    filter(country != "United States") %>% 
    pull(country) %>% as.character
  
  # Create state indicator variables for those with a prevalence below the percentage threshold 
  countryIndHelper <- function(in.data, in.state) {
    state.sym <- sym(paste0("countryInd", gsub(pattern = " ", replacement = "", x = in.state)))
    to.ret <- in.data %>% 
      mutate(!!state.sym := (country == in.state))
    return(to.ret)
  }
  
  for(country in most_reported_countries){
    iph_dat <- countryIndHelper(iph_dat, country)
  }
  
  return(iph_dat)
}

.DeriveComorbidVars <- function(iph_dat){
  # These are the only comorbidities and risk factors 
  comorbid_vars_for_tot <- paste0("comorbid___", c(1:10, 13))
  
  iph_dat$ComorbidCount <- iph_dat %>% 
    select(all_of(comorbid_vars_for_tot)) %>% 
    rowSums(., na.rm = TRUE)
  
  iph_dat %<>% 
    mutate(ComorbidTotalGrp = cut(ComorbidCount, breaks = c(0, 1, 2, Inf), include.lowest = TRUE,
                                  right = FALSE, labels = c("0", "1", "2+")))

  return(iph_dat)
}

#' A function for creating derived drug variables. Includes:
#' Two indicators of immunomodulator subgroups
.DeriveDrugVars <- function(iph_dat){
  to_return <- iph_dat %>% 
    mutate(ImmodMet = immod___1,
           Immod6mpOrAza = immod___4 | immod___5,
           preddose_dailyByWt = preddose_daily/weightBMIAbv60Rm,
           mcpdoseByWt = mcpdose/weightBMIAbv60Rm,
           azadoseByWt = azadose/weightBMIAbv60Rm)
  
  return(to_return)
}

#' A function for creating derived outcome variables. Includes:
#' an ordinal outcome for Hospitalized/ICU/Ventilated/Died
.DeriveOutcomeVars <- function(iph_dat){
  to_return <- iph_dat %>% 
    group_by(numeric_studyid) %>% 
    mutate(HospIcuVentDthOrd = case_when(
      death == TRUE ~ 4,
      vent == TRUE ~ 3,
      icu == TRUE ~ 2,
      hosp == TRUE ~ 1,
      all(c(death, vent, icu, hosp) == FALSE) ~ 0)) %>% 
    mutate_at(., .vars = "HospIcuVentDthOrd",
              ~factor(., levels = 0:4, labels = 
                        c("None", "Hospitalized", "ICU", "Ventilated", "Died"),
                      ordered = TRUE)) %>% 
    ungroup %>% 
    mutate(HospIcuDthOrd = HospIcuVentDthOrd)
  
  levels(to_return$HospIcuDthOrd) <- c("None", "Hospitalized", "ICU", "ICU", "Died")
  
  return(to_return)
}

#' Check that ages are less than or equal to 91 (data entry cap), and apply a hardcode if needed
.CheckAges <- function(iph_dat){
  to_return <- iph_dat %>% 
    mutate(age = if_else(age > 91, NA_integer_, age))
  
  to_return[to_return$numeric_studyid == 2779, "age"] <- 28
  
  return(to_return)
}