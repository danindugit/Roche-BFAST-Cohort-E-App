# turn treatment variables into factors
ADSL$ACTARM <- as.factor(ADSL$ACTARM)
ADSL$ACTARMCD <- as.factor(ADSL$ACTARMCD)
ADSL$ARM <- as.factor(ADSL$ARM)
ADSL$ARMCD <- as.factor(ADSL$ARMCD)
ADSL$TRT01A <- as.factor(ADSL$TRT01A)
ADSL$TRT01P <- as.factor(ADSL$TRT01P)
ADEX$TRT01A <- as.factor(ADEX$TRT01A)
ADEX$TRT01P <- as.factor(ADEX$TRT01P)

# function to merge ADSUB parameters to ADSL
mergeSubParam <- function(paramName, avalName){
  # filter by sub param
  ADSUB_PARAM <- subset(ADSUB, PARAMCD==paramName)
  # remove all columns except USUBJID and avalName
  keeps <- c("USUBJID", avalName)
  ADSUB_PARAM <- ADSUB_PARAM[keeps]
  # rename avalName column to paramName
  names(ADSUB_PARAM)[names(ADSUB_PARAM)==avalName] <- paramName
  # inner merge with ADSL
  inner_join(ADSL, ADSUB_PARAM, by="USUBJID")
}

# Function to merge ADEX parameters to ADEX
mergeExParam <- function(paramName, avalName){
  # filter by sub param
  ADEX_PARAM <- subset(ADEX, PARAMCD==paramName)
  # remove all columns except keys and avalName and aval
  keeps <- c("STUDYID", "USUBJID", "PARCAT1", "PARCAT2", "PARAMCD", "AVISITN", "ASTDTM", "EXSEQ", avalName)
  ADEX_PARAM <- ADEX_PARAM[keeps]
  # rename avalName column to paramName
  names(ADEX_PARAM)[names(ADEX_PARAM)==avalName] <- paramName
  # full merge with ADSL
  full_join(ADEX, ADEX_PARAM, by=c("STUDYID", "USUBJID", "PARCAT1", "PARCAT2", "PARAMCD", "AVISITN", "ASTDTM", "EXSEQ"))
}

mergeExDrugParam <- function(drugName, paramCode, ADEX_NEW){
  drugPrefix = ""
  if(drugName == "ATEZOLIZUMAB"){
    drugPrefix = "ATEZO"
  } else if(drugName == "COBIMETINIB"){
    drugPrefix = "COBI"
  } else{
    drugPrefix = "VEM"
  }

  # filter by drug name
  ADEX_DRUG <- subset(ADEX_NEW, PARCAT2==drugName)

  # remove all columns except keys and paramCode
  keeps <- c("STUDYID", "USUBJID", "PARCAT1", "PARCAT2", "PARAMCD", "AVISITN", "ASTDTM", "EXSEQ", paramCode)
  ADEX_DRUG <- ADEX_DRUG[keeps]

  # rename paramCode to drug prefix + paramCode
  names(ADEX_DRUG)[names(ADEX_DRUG)==paramCode] <- (paste(drugPrefix, paramCode, sep = ""))

  full_join(ADEX, ADEX_DRUG, by=c("STUDYID", "USUBJID", "PARCAT1", "PARCAT2", "PARAMCD", "AVISITN", "ASTDTM", "EXSEQ"))
}

# function to return a dataset that can be used to get AVAL analysis stats
getADEX_AVAL <- function(drugName, paramCode, paramLabel){
  drugPrefix = ""
  if(drugName == "ATEZOLIZUMAB"){
    drugPrefix = "ATEZO"
  } else if(drugName == "COBIMETINIB"){
    drugPrefix = "COBI"
  } else{
    drugPrefix = "VEM"
  }

  varName <- paste(drugPrefix, paramCode, sep = "")

  # merge the parameter to ADEX
  ADEX <- mergeExParam(paramCode, "AVAL")
  ADEX <- mergeExDrugParam(drugName, paramCode, ADEX)

  # keep only the rows where this parameter has data and the imc flag is set
  subset(ADEX, !is.na(ADEX[grepl(varName, colnames(ADEX))]) & ADSL$IMCFL == 'Y')
}

# merge for Baseline height
ADSL <- mergeSubParam("BHGHTCM", "AVAL")
ADSL <- var_relabel(ADSL, BHGHTCM = "Baseline Height (cm)")
# merge for baseline weight
ADSL <- mergeSubParam("BWGHTSI", "AVAL")
ADSL <- var_relabel(ADSL, BWGHTSI = "Baseline Weight (kg)")
#merge for BMI
ADSL <- mergeSubParam("BBMISI", "AVAL")
ADSL <- var_relabel(ADSL, BBMISI = "Body Mass Index (kg/m2)")
#merge for baseline ECOG
ADSL <- mergeSubParam("BECOG", "AVALC")
ADSL <- var_relabel(ADSL, BECOG = "Baseline ECOG")
#merge for tobacco use history
ADSL <- mergeSubParam("BTOBAC", "AVALC")
ADSL <- var_relabel(ADSL, BTOBAC = "Tobacco Use History")
#merge for measurable disease at BL
ADSL <- mergeSubParam("MDBASEFL", "AVALC")
ADSL <- var_relabel(ADSL, MDBASEFL = "Measurable Disease at BL (Y/N)")
#merge for tissue availability at BL
ADSL <- mergeSubParam("TISSUEFL", "AVALC")
ADSL <- var_relabel(ADSL, TISSUEFL = "Tissue Availability at Baseline")

#change description for AGEGR2
ADSL <- var_relabel(ADSL, AGEGR2 = "Age Group")

# re-order to make Y come before N in the demo table
# measurable disease at BL
ADSL$MDBASEFL <- as.factor(ADSL$MDBASEFL)
ADSL$MDBASEFL <- factor(ADSL$MDBASEFL, levels = c("Y", "N"))
# tissue availability at baseline
ADSL$TISSUEFL <- as.factor(ADSL$TISSUEFL)
ADSL$TISSUEFL <- factor(ADSL$TISSUEFL, levels = c("Y", "N"))

# IMC futility analysis flag pre-processing
ADSL$IMCFL <- c("N", "Y", "Y", "Y", "N", "Y", "Y", "Y", "Y", "N", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "N", "Y", "N", "Y", "Y", "Y", "Y")
ADSL <- var_relabel(ADSL, IMCFL = "Futility Analysis Flag")

# AE grade pre-processing
varlabel_ADAE <- var_labels(ADAE, fill = TRUE)

date_vars_ADAE <- names(ADAE)[vapply(ADAE, function(x) inherits(x, c("Date", "POSIXct", "POSIXlt")), logical(1))]

# convert variables to factor; exclude those ending in DTC
ADAE   <- df_explicit_na(data = ADAE,
                         omit_columns = c("USUBJID", "SUBJID", names(ADAE)[grep("*DTC$", names(ADAE))], date_vars_ADAE),
                         char_as_factor = TRUE
) %>% droplevels()

var_labels(ADAE) <- varlabel_ADAE

gr_grp <- list(
  "- Any Grade -" = c("1", "2", "3", "4", "5"),
  "Grade 1-2" = c("1", "2"),
  "Grade 3-4" = c("3", "4"),
  "Grade 5" = "5"
)

raw_table <- function(adae, adsl, gr_grp) {
  split_fun <- trim_levels_in_group
  
  lyt <- basic_table() %>%
    split_cols_by("ACTARM") %>%
    add_colcounts() %>%
    summarize_occurrences_by_grade(
      var = "AETOXGR",
      grade_groups = gr_grp
    ) %>%
    split_rows_by(
      "AEBODSYS",
      child_labels = "visible",
      nested = TRUE,
      indent_mod = -1L,
      split_fun = split_fun("AETOXGR"),
      label_pos = "topleft",
      split_label = obj_label(adae$AEBODSYS)
    ) %>%
    summarize_occurrences_by_grade(
      var = "AETOXGR",
      grade_groups = gr_grp
    ) %>%
    split_rows_by(
      "AEDECOD",
      child_labels = "visible",
      nested = TRUE,
      indent_mod = -1L,
      split_fun = split_fun("AETOXGR"),
      label_pos = "topleft",
      split_label = obj_label(adae$AEDECOD)
    ) %>%
    summarize_num_patients(
      var = "USUBJID",
      .stats = "unique",
      .labels = "- Any Grade -"
    ) %>%
    count_occurrences_by_grade(
      var = "AETOXGR",
      grade_groups = gr_grp[-1],
      .indent_mods = -1L
    ) %>%
    append_varlabels(adae, "AETOXGR", indent = 2L)
  
  result <- lyt %>%
    build_table(adae, alt_counts_df = adsl) %>%
    prune_table() %>%
    sort_at_path(
      path = "AEBODSYS",
      scorefun = cont_n_allcols,
      decreasing = TRUE
    ) %>%
    sort_at_path(
      path = c("AEBODSYS", "*", "AEDECOD"),
      scorefun = cont_n_allcols,
      decreasing = TRUE
    )
  
  result
}

# Simple wrapper to return subset ADAE to a threshold of xx%.
get_adae_trimmed <- function(adsl, adae, cutoff_rate) {
  n_per_arm <- adsl %>%
    dplyr::count(ACTARM)
  
  anl_terms <- adae %>%
    dplyr::group_by(ACTARM, AEBODSYS, AEDECOD) %>%
    dplyr::count(
      unique_terms = n_distinct(USUBJID)
    ) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup()
  
  anl_terms <- dplyr::left_join(
    anl_terms,
    n_per_arm,
    by = "ACTARM"
  ) %>%
    dplyr::mutate(
      ae_rate = unique_terms / n
    ) %>%
    dplyr::filter(ae_rate >= cutoff_rate) %>%
    dplyr::select(AEDECOD) %>%
    unique()
  
  anl <- dplyr::left_join(
    anl_terms,
    adae,
    by = "AEDECOD"
  )
  anl
}


# events summary pre-processing
add_event_flags <- function(dat) {
  dat %>%
    dplyr::mutate(
      TMPFL_REL = AEREL == "Y",
      TMPFL_GR34 = AETOXGR == "3" | AETOXGR == "4",
      TMPFL_GR34_REL = (AETOXGR == "3" | AETOXGR == "4") & AEREL == "Y",
      TMPFL_GR5 = AETOXGR == "5",
      TMPFL_GR5_REL = AETOXGR == "5" & AEREL == "Y",
      TMPFL_SER = AESER == "Y",
      TMPFL_SER_REL = AESER == "Y" & AEREL == "Y",
      TMPFL_WD = AEACN1 == "DRUG WITHDRAWN" | AEACN2 == "DRUG WITHDRAWN" | AEACN3 == "DRUG WITHDRAWN",
      TMPFL_WD1 = AEACN3 == "DRUG WITHDRAWN",
      TMPFL_WD2 = AEACN1 == "DRUG WITHDRAWN",
      TMPFL_WD3 = AEACN2 == "DRUG WITHDRAWN",
      TMPFL_RD = AEACN1 == "DOSE REDUCED" | AEACN2 == "DOSE REDUCED" | AEACN3 == "DOSE REDUCED",
      TMPFL_RD1 = AEACN3 == "DOSE REDUCED",
      TMPFL_RD2 = AEACN1 == "DOSE REDUCED",
      TMPFL_RD3 = AEACN2 == "DOSE REDUCED",
      TMPFL_INT = AEACN1 == "DRUG INTERRUPTED" | AEACN2 == "DRUG INTERRUPTED" | AEACN3 == "DRUG INTERRUPTED",
      TMPFL_INT1 = AEACN3 == "DRUG INTERRUPTED",
      TMPFL_INT2 = AEACN1 == "DRUG INTERRUPTED",
      TMPFL_INT3 = AEACN2 == "DRUG INTERRUPTED",
      TMP_SMQ02 = !is.na(SMQ02NAM),
      TMP_CQ01 = !is.na(CQ01NAM)
    ) %>%
    var_relabel(
      TMPFL_REL = "Treatment-related AE",
      TMPFL_GR34 = "Grade 3-4 AE",
      TMPFL_GR34_REL = "Treatment-related Grade 3-4 AE",
      TMPFL_GR5 = "Grade 5 AE",
      TMPFL_GR5_REL = "Treatment-related Grade 5 AE",
      TMPFL_SER = "Serious AE",
      TMPFL_SER_REL = "Treatment-related Serious AE",
      TMPFL_WD = "AE leading to withdrawal",
      TMPFL_WD1 = "- Atezolizumab",
      TMPFL_WD2 = "- Cobimetinib",
      TMPFL_WD3 = "- Vemurafenib",
      TMPFL_RD = "AE leading to dose reduction",
      TMPFL_RD1 = "- Atezolizumab",
      TMPFL_RD2 = "- Cobimetinib",
      TMPFL_RD3 = "- Vemurafenib",
      TMPFL_INT = "AE leading to dose interruption",
      TMPFL_INT1 = "- Atezolizumab",
      TMPFL_INT2 = "- Cobimetinib",
      TMPFL_INT3 = "- Vemurafenib"
    )
}

# Generating user-defined event flags.
ADAE <- ADAE %>% add_event_flags()

# ADEX pre-processing

ADEX <- ADEX %>%
  mutate(
    AVALCAT1 = case_when(
      PARAMCD == "TNCYC" & AVAL <= 3 ~ "0 to <= 3",
      PARAMCD == "TNCYC" & AVAL > 3 & AVAL <= 6 ~ "> 3 to <= 6",
      PARAMCD == "TNCYC" & AVAL > 6 & AVAL <= 12 ~ "> 6 - <= 12",
      PARAMCD == "TNCYC" & AVAL > 12 ~ "> 12",
      PARAMCD == "TNDOSE" & AVAL <= 3 ~ "0 to <= 3",
      PARAMCD == "TNDOSE" & AVAL > 3 & AVAL <= 6 ~ "> 3 to <= 6",
      PARAMCD == "TNDOSE" & AVAL > 6 & AVAL <= 12 ~ "> 6 - <= 12",
      PARAMCD == "TNDOSE" & AVAL > 12 ~ "> 12",
      PARAMCD == "TDURD" & AVAL <= 3 ~ "0 to <= 3",
      PARAMCD == "TDURD" & AVAL > 3 & AVAL <= 6 ~ "> 3 to <= 6",
      PARAMCD == "TDURD" & AVAL > 6 & AVAL <= 12 ~ "> 6 - <= 12",
      PARAMCD == "TDURD" & AVAL > 12 ~ "> 12",
      PARAMCD == "TIADJ" & AVALC == "Y" ~ "Total number of patients with at least one infusion modification",
      PARAMCD == "TIADJ" & AVALC == "N" ~ "Total number of patients with no infusion modifications"
      # AVALC != "Y" ~ AVALC
    ),
    # AVALCAT1 = na_if(AVALCAT1, "N"),
    AVAL = case_when(
      PARAMCD == "TDURD" ~ AVAL / 30.4375,
      PARAMCD != "TDURD" ~ AVAL
    ),
    PARAM = case_when(
      PARAMCD == "TDURD" ~ "Treatment duration (months)",
      PARAMCD != "TDURD" ~ PARAM
    )
  )

# ordering the yesses before nos
ADEX$AVALC <- as.factor(ADEX$AVALC)
ADEX$AVALC <- factor(ADEX$AVALC, levels = c("Y", "N"))

# ordering all AVALCAT rows
ADEX$AVALCAT1 <- as.factor(ADEX$AVALCAT1)
ADEX$AVALCAT1 <- factor(ADEX$AVALCAT1, levels = c("Total number of patients with at least one infusion modification", "Total number of patients with no infusion modifications", "0 to <= 3", "> 3 to <= 6", "> 6 - <= 12", "> 12"))

ADEX <- var_relabel(ADEX, AVALCAT1 = "Categorical Analysis")
ADEX <- var_relabel(ADEX, AVAL = "Numerical Analysis")
ADEX <- var_relabel(ADEX, PARAM = "Parameter")

# Disposition pre-processing
ADSL <- ADSL %>%
  mutate(
    EOSSTT = ifelse(EOSSTT == "Discontinued", EOSSTT, "Ongoing"),
    EOT1STT = ifelse(EOT1STT == "Discontinued", EOT1STT, "Ongoing"),
    EOT2STT = ifelse(EOT2STT == "Discontinued", EOT2STT, "Ongoing"),
    EOT3STT = ifelse(EOT3STT == "Discontinued", EOT3STT, "Ongoing")
  )

# ordering the OTHER last
ADSL$DCT1RS <- as.factor(ADSL$DCT1RS)
ADSL$DCT2RS <- as.factor(ADSL$DCT2RS)
ADSL$DCT3RS <- as.factor(ADSL$DCT3RS)
ADSL$DCT1RS <- factor(ADSL$DCT1RS, levels = c("ADVERSE EVENT", "DEATH", "PROGRESSIVE DISEASE", "OTHER"))
ADSL$DCT2RS <- factor(ADSL$DCT2RS, levels = c("ADVERSE EVENT", "DEATH", "PROGRESSIVE DISEASE", "WITHDRAWAL BY SUBJECT", "OTHER"))
ADSL$DCT3RS <- factor(ADSL$DCT3RS, levels = c("ADVERSE EVENT", "DEATH", "PROGRESSIVE DISEASE", "WITHDRAWAL BY SUBJECT", "OTHER"))

# relabeling variables
ADSL <- var_relabel(ADSL, EOT1STT = "Atezo End of Treatment Status")
ADSL <- var_relabel(ADSL, EOT2STT = "Cobi End of Treatment Status")
ADSL <- var_relabel(ADSL, EOT3STT = "Vemura End of Treatment Status")
ADSL <- var_relabel(ADSL, DCT1RS = "Atezo Discontinuation Reason")
ADSL <- var_relabel(ADSL, DCT2RS = "Cobi Discontinuation Reason")
ADSL <- var_relabel(ADSL, DCT3RS = "Vemura Discontinuation Reason")

############ Pre-processing ADAE Barplot ############
# Extracting Relevant Variables from ADAE
vars_adae_bp <- c("USUBJID", "ASTDT", "TRTSDT", "AESTDTC", "AESTDY", "AETERM",
                  "AEDECOD", "AEPTCD", "AETOXGR", "AEITOXGR", "ATOXGR")

ADAE_BP <- ADAE[vars_adae_bp]

# Make ASTDT a Date Variable
ADAE_BP$ASTDT <- as.Date(ADAE_BP$ASTDT)

# Fill NA for AESTDY as 0 (Pre-Treatment)
ADAE_BP$AESTDY[is.na(ADAE_BP$AESTDY)] <- 0

# Create Cycle Variable (X-AXIS)
ADAE_BP$CYCLE <- ceiling(ADAE_BP$AESTDY / 21) # 21 comes from the Cycle Length from Protocol

# Proportion of AE per Cycle with Max Toxicity Value per Unique Subject ID
ADAE_BP <- ADAE_BP %>%
  arrange(desc(AETOXGR)) %>%
  distinct(USUBJID, CYCLE, .keep_all = TRUE) %>%
  count(CYCLE, AETOXGR) %>%
  group_by(CYCLE) %>%
  mutate(PROP = sum(n)/26) %>% # 26 comes from 26 Unique USUBJIDs
  mutate(AEPROP = n/sum(n) * PROP)

# Butterfly plot pre-processing
ADAE <- ADAE %>%
  mutate(
    AEREL = ifelse(AEREL == "Y", "Y", "N")
  )
