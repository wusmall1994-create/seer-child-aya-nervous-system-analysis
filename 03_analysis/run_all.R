#!/usr/bin/env Rscript

# Reproducible analysis pipeline for childhood/AYA cancer survivors and
# subsequent malignant brain and other nervous system tumors in SEER MP-SIR data.
# Raw files are read-only. Patient ID is used only in memory and is never exported.

options(stringsAsFactors = FALSE, scipen = 999)

required_packages <- c("ggplot2", "patchwork", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
project_dir <- Sys.getenv("SEER_PROJECT_DIR", unset = dirname(dirname(script_path)))
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "02_raw_exports")
analysis_dir <- file.path(project_dir, "03_analysis")
results_dir <- file.path(project_dir, "04_results")
tables_dir <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
source_dir <- file.path(results_dir, "source_data")
logs_dir <- file.path(project_dir, "05_logs")
for (d in c(analysis_dir, tables_dir, figures_dir, source_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  main = "01_main_first_cancer_latency_mpsir.txt",
  age_main = "02_age_at_diagnosis_latency_mpsir.txt",
  sex = "03_sex_latency_mpsir.txt",
  attained_age = "04_attained_age_latency_mpsir.txt",
  site_age = "05_first_cancer_diag_age_latency_mpsir.txt",
  child_iccc = "06_child_iccc_latency_mpsir.txt",
  aya = "07_aya_recode_latency_mpsir.txt",
  age_12m = "08_age_latency_12m_sensitivity_mpsir.txt",
  age_60m = "09_age_latency_60m_fiveyear_survivor_mpsir.txt",
  year_dx = "10_year_dx_latency_mpsir.txt",
  race = "11_race_latency_mpsir.txt",
  radiation = "12_radiation_latency_mpsir.txt",
  cases = "13_brain_cns_event_case_listing.txt"
)

missing_inputs <- input_files[!file.exists(file.path(raw_dir, input_files))]
if (length(missing_inputs)) stop("Missing input files: ", paste(missing_inputs, collapse = ", "))

read_seer <- function(filename, patient_id_character = FALSE) {
  path <- file.path(raw_dir, filename)
  if (patient_id_character) {
    read.delim(path, check.names = FALSE, na.strings = c("NA"), quote = "\"",
               colClasses = c("Patient ID" = "character"))
  } else {
    read.delim(path, check.names = FALSE, na.strings = c("NA"), quote = "\"")
  }
}

dat <- lapply(names(input_files), function(nm) read_seer(input_files[[nm]], nm == "cases"))
names(dat) <- names(input_files)

write_csv <- function(x, filename) {
  if ("Patient ID" %in% names(x)) stop("Privacy stop: attempted export of Patient ID.")
  write.csv(x, file.path(tables_dir, filename), row.names = FALSE, na = "")
}

write_source <- function(x, filename) {
  if ("Patient ID" %in% names(x)) stop("Privacy stop: attempted source-data export of Patient ID.")
  write.csv(x, file.path(source_dir, filename), row.names = FALSE, na = "")
}

exact_ci <- function(observed, expected, level = 0.95) {
  alpha <- 1 - level
  lower <- rep(NA_real_, length(observed))
  upper <- rep(NA_real_, length(observed))
  ok <- expected > 0
  lower[ok & observed == 0] <- 0
  lower[ok & observed > 0] <- 0.5 * qchisq(alpha / 2, 2 * observed[ok & observed > 0]) / expected[ok & observed > 0]
  upper[ok] <- 0.5 * qchisq(1 - alpha / 2, 2 * (observed[ok] + 1)) / expected[ok]
  data.frame(Lower = lower, Upper = upper)
}

add_sir <- function(x) {
  x$SIR <- ifelse(x$Expected > 0, x$Observed / x$Expected, NA_real_)
  ci <- exact_ci(x$Observed, x$Expected)
  x$CI_Lower <- ci$Lower
  x$CI_Upper <- ci$Upper
  if ("Person Years at Risk" %in% names(x)) {
    x$EAR_per_10000_PY <- ifelse(
      x$`Person Years at Risk` > 0,
      (x$Observed - x$Expected) / x$`Person Years at Risk` * 10000,
      NA_real_
    )
  }
  x
}

aggregate_sir <- function(x, groups) {
  if (!nrow(x)) return(data.frame())
  measure_vars <- c("Observed", "Expected")
  for (v in c("Persons", "Person Years at Risk")) {
    if (v %in% names(x)) measure_vars <- c(measure_vars, v)
  }
  z <- aggregate(x[, measure_vars, drop = FALSE], x[, groups, drop = FALSE], sum, na.rm = TRUE)
  add_sir(z)
}

poisson_heterogeneity <- function(x, group_var) {
  x <- x[is.finite(x$Expected) & x$Expected > 0, , drop = FALSE]
  if (nrow(x) < 2L || length(unique(x[[group_var]])) < 2L) return(NA_real_)
  m0 <- glm(Observed ~ 1, poisson, offset = log(Expected), data = x)
  m1 <- glm(Observed ~ factor(x[[group_var]]), poisson, offset = log(Expected), data = x)
  anova(m0, m1, test = "Chisq")[2, "Pr(>Chi)"]
}

poisson_trend <- function(x, ordinal) {
  x$ordinal_internal <- ordinal
  x <- x[is.finite(x$Expected) & x$Expected > 0, , drop = FALSE]
  if (nrow(x) < 2L) return(NA_real_)
  m <- glm(Observed ~ ordinal_internal, poisson, offset = log(Expected), data = x)
  coef(summary(m))["ordinal_internal", "Pr(>|z|)"]
}

parse_age <- function(x) suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", x)))

diagnosis_age_group <- function(x) {
  a <- parse_age(x)
  out <- ifelse(a <= 14, "0-14", ifelse(a <= 19, "15-19", ifelse(a <= 29, "20-29", ifelse(a <= 39, "30-39", NA))))
  factor(out, levels = c("0-14", "15-19", "20-29", "30-39"))
}

outcome_target <- "Brain and Other Nervous System"
outcome_brain <- "Brain"
latency_order <- c("2-11 months", "12-59 months", "60-119 months", "120+ months", "Total")

# ---- Main outcome and latency table ----
main <- dat$main
main_row_var <- names(main)[3]
main_summary <- main[main[[main_row_var]] == "All Sites", ]
main_summary$Latency <- factor(main_summary$Latency, levels = latency_order)
main_summary <- main_summary[order(main_summary$`Selected Events`, main_summary$Latency), ]
main_summary <- add_sir(main_summary[, c("Selected Events", "Latency", "Observed", "Expected", "Persons", "Person Years at Risk")])
names(main_summary)[1] <- "Outcome"
write_csv(main_summary, "table_01_main_outcome_latency.csv")

# ---- Diagnosis-age analyses ----
make_age_table <- function(d, analysis_name, outcome = outcome_target) {
  x <- d[d$`Selected Events` == outcome, ]
  x$Age_Group <- diagnosis_age_group(x[[3]])
  x <- x[!is.na(x$Age_Group), ]
  total <- aggregate_sir(x[x$Latency == "Total", ], "Age_Group")
  total$Analysis <- analysis_name
  by_latency <- aggregate_sir(x[x$Latency != "Total", ], c("Age_Group", "Latency"))
  by_latency$Analysis <- analysis_name
  list(total = total, by_latency = by_latency)
}

age_main <- make_age_table(dat$age_main, "2-month exclusion")
age_12m <- make_age_table(dat$age_12m, "12-month exclusion")
age_60m <- make_age_table(dat$age_60m, "60-month exclusion")
age_brain <- make_age_table(dat$age_main, "Brain-only (2-month exclusion)", outcome_brain)
age_sensitivity <- rbind(age_main$total, age_12m$total, age_60m$total)
age_sensitivity$Age_Group <- factor(age_sensitivity$Age_Group, levels = c("0-14", "15-19", "20-29", "30-39"))
age_sensitivity <- age_sensitivity[order(age_sensitivity$Analysis, age_sensitivity$Age_Group), ]
write_csv(age_main$total, "table_02_age_main.csv")
write_csv(age_sensitivity, "table_03_age_sensitivity.csv")
write_csv(age_main$by_latency, "table_04_age_by_latency.csv")
write_csv(age_brain$total, "table_25_age_brain_only.csv")
write_csv(age_brain$by_latency, "table_26_age_brain_only_by_latency.csv")

# ---- Sex ----
sex <- dat$sex
sex <- sex[sex$`Selected Events` == outcome_target & sex$Latency == "Total" & sex$Sex %in% c("Male", "Female"), ]
sex_table <- add_sir(sex[, c("Sex", "Observed", "Expected", "Persons", "Person Years at Risk")])
write_csv(sex_table, "table_05_sex.csv")

# ---- Attained age ----
att <- dat$attained_age
att <- att[att$`Selected Events` == outcome_target & att$Latency == "Total", ]
att <- att[att$Expected > 0 | att$Observed > 0, ]
att_table <- add_sir(att[, c("Attained Age", "Observed", "Expected", "Persons", "Person Years at Risk")])
write_csv(att_table, "table_06_attained_age.csv")

# ---- Diagnosis era ----
yd <- dat$year_dx
year_numeric <- suppressWarnings(as.integer(yd$`Year of diagnosis`))
yd <- yd[!is.na(year_numeric) & yd$`Selected Events` == outcome_target, ]
year_numeric <- as.integer(yd$`Year of diagnosis`)
yd$Era <- cut(year_numeric, breaks = c(1974, 1989, 1999, 2009, 2023),
              labels = c("1975-1989", "1990-1999", "2000-2009", "2010-2023"))
era_total <- aggregate_sir(yd[yd$Latency == "Total", ], "Era")
era_latency <- aggregate_sir(yd[yd$Latency != "Total", ], c("Era", "Latency"))
write_csv(era_total, "table_07_diagnosis_era.csv")
write_csv(era_latency, "table_08_diagnosis_era_by_latency.csv")

# ---- Race ----
race_var <- "Race recode (White, Black, Other)"
race <- dat$race
race <- race[race$`Selected Events` == outcome_target & race[[race_var]] != "All races", ]
race_total <- add_sir(race[race$Latency == "Total", c(race_var, "Observed", "Expected", "Persons", "Person Years at Risk")])
race_latency <- add_sir(race[race$Latency != "Total", c(race_var, "Latency", "Observed", "Expected", "Persons", "Person Years at Risk")])
write_csv(race_total, "table_09_race.csv")
write_csv(race_latency, "table_10_race_by_latency.csv")

# ---- Radiation ----
rt <- dat$radiation
rt_value <- rt$`Radiation recode`
rt$RT_Group <- ifelse(rt_value == "None/Unknown", "None/Unknown",
  ifelse(rt_value == "Beam radiation", "Beam radiation",
    ifelse(rt_value %in% c("Refused (1988+)", "Recommended, unknown if administered"),
           "Uncertain", "Other documented RT")))
rt <- rt[rt$`Selected Events` == outcome_target, ]
rt_total <- aggregate_sir(rt[rt$Latency == "Total", ], "RT_Group")
rt$RT_Binary <- ifelse(rt$RT_Group %in% c("Beam radiation", "Other documented RT"),
                       "Any documented RT", rt$RT_Group)
rt_binary <- aggregate_sir(rt[rt$Latency == "Total" & rt$RT_Binary != "Uncertain", ], "RT_Binary")
rt_latency <- aggregate_sir(rt[rt$Latency != "Total" & rt$RT_Binary != "Uncertain", ], c("RT_Binary", "Latency"))
write_csv(rt_total, "table_11_radiation_categories.csv")
write_csv(rt_binary, "table_12_radiation_binary.csv")
write_csv(rt_latency, "table_13_radiation_by_latency.csv")

# ---- Child ICCC top-level categories ----
iccc <- dat$child_iccc
iccc_var <- names(iccc)[3]
iccc_top_regex <- "^(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII) "
iccc$is_top <- grepl(iccc_top_regex, iccc[[iccc_var]])
iccc$is_cns <- grepl("^III ", iccc[[iccc_var]])
iccc <- iccc[iccc$`Selected Events` == outcome_target & iccc$is_top & !iccc$is_cns, ]
iccc_total <- add_sir(iccc[iccc$Latency == "Total", c(iccc_var, "Observed", "Expected", "Persons", "Person Years at Risk")])
names(iccc_total)[1] <- "Primary_Cancer_ICCC"
iccc_latency <- add_sir(iccc[iccc$Latency != "Total", c(iccc_var, "Latency", "Observed", "Expected", "Persons", "Person Years at Risk")])
names(iccc_latency)[1] <- "Primary_Cancer_ICCC"
write_csv(iccc_total, "table_14_child_iccc.csv")
write_csv(iccc_latency, "table_15_child_iccc_by_latency.csv")

# ---- AYA top-level categories ----
aya <- dat$aya
aya_var <- names(aya)[3]
aya$is_top <- grepl("^[0-9]+\\. [A-Z]", aya[[aya_var]])
aya$is_cns <- grepl("^3\\. ", aya[[aya_var]])
aya <- aya[aya$`Selected Events` == outcome_target & aya$is_top & !aya$is_cns, ]
aya_total <- add_sir(aya[aya$Latency == "Total", c(aya_var, "Observed", "Expected", "Persons", "Person Years at Risk")])
names(aya_total)[1] <- "Primary_Cancer_AYA"
aya_latency <- add_sir(aya[aya$Latency != "Total", c(aya_var, "Latency", "Observed", "Expected", "Persons", "Person Years at Risk")])
names(aya_latency)[1] <- "Primary_Cancer_AYA"
write_csv(aya_total, "table_16_aya_primary_cancer.csv")
write_csv(aya_latency, "table_17_aya_primary_cancer_by_latency.csv")

# ---- Cohort composition for supplementary reporting ----
make_composition <- function(x, characteristic, category_var) {
  data.frame(
    Characteristic = characteristic,
    Category = as.character(x[[category_var]]),
    Persons = x$Persons,
    Percent = 100 * x$Persons / sum(x$Persons),
    Person_Years = x$`Person Years at Risk`,
    Mean_Followup_Years = x$`Person Years at Risk` / x$Persons,
    stringsAsFactors = FALSE
  )
}

age_comp <- make_composition(age_main$total, "Age at first cancer", "Age_Group")
sex_comp <- make_composition(sex_table, "Sex", "Sex")
era_comp <- make_composition(era_total, "Diagnosis era", "Era")

child_comp <- data.frame(
  Characteristic = "First-cancer group: age 0-14 years",
  Category = as.character(iccc_total$Primary_Cancer_ICCC),
  Persons = iccc_total$Persons,
  Person_Years = iccc_total$`Person Years at Risk`,
  stringsAsFactors = FALSE
)
child_comp <- rbind(child_comp, data.frame(
  Characteristic = "First-cancer group: age 0-14 years", Category = "Other/unclassified",
  Persons = age_main$total$Persons[age_main$total$Age_Group == "0-14"] - sum(child_comp$Persons),
  Person_Years = age_main$total$`Person Years at Risk`[age_main$total$Age_Group == "0-14"] - sum(child_comp$Person_Years),
  stringsAsFactors = FALSE
))
child_comp$Percent <- 100 * child_comp$Persons / sum(child_comp$Persons)
child_comp$Mean_Followup_Years <- child_comp$Person_Years / child_comp$Persons

aya_comp <- data.frame(
  Characteristic = "First-cancer group: age 15-39 years",
  Category = as.character(aya_total$Primary_Cancer_AYA),
  Persons = aya_total$Persons,
  Person_Years = aya_total$`Person Years at Risk`,
  stringsAsFactors = FALSE
)
aya_target_persons <- sum(age_main$total$Persons[age_main$total$Age_Group %in% c("15-19", "20-29", "30-39")])
aya_target_py <- sum(age_main$total$`Person Years at Risk`[age_main$total$Age_Group %in% c("15-19", "20-29", "30-39")])
aya_comp <- rbind(aya_comp, data.frame(
  Characteristic = "First-cancer group: age 15-39 years", Category = "Other/unclassified",
  Persons = aya_target_persons - sum(aya_comp$Persons),
  Person_Years = aya_target_py - sum(aya_comp$Person_Years),
  stringsAsFactors = FALSE
))
aya_comp$Percent <- 100 * aya_comp$Persons / sum(aya_comp$Persons)
aya_comp$Mean_Followup_Years <- aya_comp$Person_Years / aya_comp$Persons

cohort_composition <- rbind(age_comp, sex_comp, era_comp,
                            child_comp[, names(age_comp)], aya_comp[, names(age_comp)])
write_csv(cohort_composition, "table_27_cohort_composition.csv")

# ---- Age heterogeneity adjusted for mutually exclusive first-cancer site groups ----
site_age <- dat$site_age
site_age_var <- names(site_age)[4]
top_site_groups <- c(
  "Oral Cavity and Pharynx", "Digestive System", "Respiratory System",
  "Bones and Joints", "Soft Tissue including Heart",
  "Skin excluding Basal and Squamous", "Breast", "Female Genital System",
  "Male Genital System", "Urinary System", "Eye and Orbit", "Endocrine System",
  "Lymphoma", "Myeloma", "Leukemia", "Mesothelioma", "Kaposi Sarcoma", "Miscellaneous"
)
site_age_adjust <- site_age[
  site_age$`Selected Events` == outcome_target & site_age$Latency == "Total" &
    site_age[[site_age_var]] %in% top_site_groups,
]
site_age_adjust$Age_Group <- diagnosis_age_group(site_age_adjust[[3]])
site_age_adjust <- site_age_adjust[!is.na(site_age_adjust$Age_Group), ]
site_age_adjust$First_Cancer_Site_Group <- factor(site_age_adjust[[site_age_var]], levels = top_site_groups)
site_age_adjust <- aggregate_sir(site_age_adjust, c("Age_Group", "First_Cancer_Site_Group"))
site_age_adjust <- site_age_adjust[site_age_adjust$Persons > 0, ]
site_age_adjust$Age_Group <- relevel(factor(site_age_adjust$Age_Group), ref = "30-39")
stopifnot(
  sum(site_age_adjust$Observed) == sum(age_main$total$Observed),
  sum(site_age_adjust$Persons) == sum(age_main$total$Persons)
)
write_csv(site_age_adjust, "table_28_age_by_first_cancer_site.csv")

site_age_model_data <- site_age_adjust[site_age_adjust$Expected > 0, ]
site_only_model <- glm(Observed ~ First_Cancer_Site_Group, poisson,
                       offset = log(Expected), data = site_age_model_data)
site_age_model <- glm(Observed ~ First_Cancer_Site_Group + Age_Group, poisson,
                      offset = log(Expected), data = site_age_model_data)
age_adjusted_p <- anova(site_only_model, site_age_model, test = "Chisq")[2, "Pr(>Chi)"]
age_adjusted_coef <- coef(summary(site_age_model))
age_adjusted_rr_ci <- function(term) {
  beta <- age_adjusted_coef[term, "Estimate"]
  se <- age_adjusted_coef[term, "Std. Error"]
  c(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
}
rr_0_14 <- age_adjusted_rr_ci("Age_Group0-14")
rr_15_19 <- age_adjusted_rr_ci("Age_Group15-19")
rr_20_29 <- age_adjusted_rr_ci("Age_Group20-29")
age_adjusted_effects <- data.frame(
  Age_Group = c("0-14", "15-19", "20-29", "30-39"),
  Rate_Ratio_vs_30_39 = c(rr_0_14[1], rr_15_19[1], rr_20_29[1], 1),
  CI_Lower = c(rr_0_14[2], rr_15_19[2], rr_20_29[2], NA_real_),
  CI_Upper = c(rr_0_14[3], rr_15_19[3], rr_20_29[3], NA_real_),
  Observed_Events = as.integer(tapply(site_age_model_data$Observed,
                                      site_age_model_data$Age_Group, sum)[c("0-14", "15-19", "20-29", "30-39")]),
  Model_Persons = as.integer(tapply(site_age_model_data$Persons,
                                    site_age_model_data$Age_Group, sum)[c("0-14", "15-19", "20-29", "30-39")]),
  stringsAsFactors = FALSE
)
age_adjusted_effects$Overall_Age_Heterogeneity_P <- c(age_adjusted_p, rep(NA_real_, 3))
write_csv(age_adjusted_effects, "table_29_age_adjusted_effects.csv")

# ---- Event case listing: aggregate only, never export Patient ID ----
cases <- dat$cases
events <- cases[cases$`Event Number` == "1", ]
indices <- cases[cases$`Event Number` == "Index Record", ]

index_age <- parse_age(indices$`Age recode with single ages and 90+`)
index_age_map <- setNames(as.character(diagnosis_age_group(indices$`Age recode with single ages and 90+`)),
                          as.character(indices$`Patient ID`))
events$Index_Age_Group <- unname(index_age_map[as.character(events$`Patient ID`)])

histology_group <- function(code) {
  ifelse(code %in% c(9440, 9442, 9445), "Glioblastoma/gliosarcoma",
  ifelse(code %in% c(9400, 9401, 9411, 9420, 9421, 9424), "Other astrocytic tumors",
  ifelse(code %in% c(9450, 9451), "Oligodendroglial tumors",
  ifelse(code %in% c(9380, 9381, 9382, 9385), "Mixed/unspecified diffuse glioma",
  ifelse(code %in% c(9390, 9391, 9392), "Ependymal/choroid plexus tumors",
  ifelse(code %in% c(9470, 9473, 9474, 9500, 9508), "Embryonal/neuroblastic tumors",
  ifelse(code %in% c(9530, 9538), "Malignant meningioma",
  ifelse(code %in% c(8680, 9540, 9560, 9561), "Peripheral nerve sheath/paraganglial tumors",
  ifelse(code %in% c(8800, 8801, 8806, 8963, 9120), "Sarcoma/vascular/rhabdoid tumors",
  ifelse(code == 8000, "Malignant neoplasm NOS", "Other"))))))))))
}

events$Histology_Group <- histology_group(as.integer(events$`Histologic Type ICD-O-3`))
hist_detail <- as.data.frame(sort(table(events$`ICD-O-3 Hist/behav`), decreasing = TRUE), stringsAsFactors = FALSE)
names(hist_detail) <- c("Histology", "N")
hist_detail$Percent <- 100 * hist_detail$N / sum(hist_detail$N)
hist_group <- as.data.frame(sort(table(events$Histology_Group), decreasing = TRUE), stringsAsFactors = FALSE)
names(hist_group) <- c("Histology_Group", "N")
hist_group$Percent <- 100 * hist_group$N / sum(hist_group$N)
hist_group_latency <- aggregate(as.numeric(events$`Months Since Index (Calculated)`),
                                list(Histology_Group = events$Histology_Group),
                                function(x) c(N = length(x), Median = median(x), Q1 = unname(quantile(x, .25)), Q3 = unname(quantile(x, .75))))
hist_group_latency <- data.frame(Histology_Group = hist_group_latency$Histology_Group,
                                 N = hist_group_latency$x[, "N"],
                                 Median_Latency_Months = hist_group_latency$x[, "Median"],
                                 Q1_Latency_Months = hist_group_latency$x[, "Q1"],
                                 Q3_Latency_Months = hist_group_latency$x[, "Q3"])
hist_by_index_age <- as.data.frame.matrix(table(events$Histology_Group,
                                                 factor(events$Index_Age_Group, levels = c("0-14", "15-19", "20-29", "30-39"))))
hist_by_index_age$Histology_Group <- rownames(hist_by_index_age)
rownames(hist_by_index_age) <- NULL
hist_by_index_age <- hist_by_index_age[, c("Histology_Group", "0-14", "15-19", "20-29", "30-39")]

site_table <- as.data.frame(sort(table(events$`Primary Site - labeled`), decreasing = TRUE), stringsAsFactors = FALSE)
names(site_table) <- c("Primary_Site", "N")
site_table$Percent <- 100 * site_table$N / sum(site_table$N)

event_age <- parse_age(events$`Age recode with single ages and 90+`)
event_year <- as.numeric(events$`Year of diagnosis`)
event_latency <- as.numeric(events$`Months Since Index (Calculated)`)
event_characteristics <- data.frame(
  Characteristic = c("Unique patients", "Event records", "Brain events", "Cranial nerves/other nervous system events",
                     "Event age, median", "Event age, Q1", "Event age, Q3",
                     "Latency months, median", "Latency months, Q1", "Latency months, Q3",
                     "Event diagnosis year, median"),
  Value = c(length(unique(cases$`Patient ID`)), nrow(events),
            sum(events$`Site recode ICD-O-3/WHO 2008 (for SIRs) (Event Variable)` == "Brain"),
            sum(events$`Site recode ICD-O-3/WHO 2008 (for SIRs) (Event Variable)` == "Cranial Nerves Other Nervous System"),
            median(event_age), unname(quantile(event_age, .25)), unname(quantile(event_age, .75)),
            median(event_latency), unname(quantile(event_latency, .25)), unname(quantile(event_latency, .75)),
            median(event_year))
)

write_csv(hist_detail, "table_18_event_histology_detail.csv")
write_csv(hist_group, "table_19_event_histology_groups.csv")
write_csv(hist_group_latency, "table_20_event_histology_latency.csv")
write_csv(hist_by_index_age, "table_21_event_histology_by_index_age.csv")
write_csv(site_table, "table_22_event_primary_site.csv")
write_csv(event_characteristics, "table_23_event_characteristics.csv")

# ---- Statistical tests ----
tests <- list()
add_test <- function(name, p, method, note = "") {
  tests[[length(tests) + 1L]] <<- data.frame(Test = name, P_value = p, Method = method, Note = note)
}
add_test("Diagnosis age heterogeneity, 2-month exclusion", poisson_heterogeneity(age_main$total, "Age_Group"), "Poisson likelihood-ratio test")
add_test("Diagnosis age heterogeneity, 12-month exclusion", poisson_heterogeneity(age_12m$total, "Age_Group"), "Poisson likelihood-ratio test")
add_test("Diagnosis age heterogeneity, 60-month exclusion", poisson_heterogeneity(age_60m$total, "Age_Group"), "Poisson likelihood-ratio test")
add_test("Diagnosis age heterogeneity, brain-only 2-month exclusion", poisson_heterogeneity(age_brain$total, "Age_Group"), "Poisson likelihood-ratio test", "Sensitivity analysis restricted to Brain sites")
add_test("Diagnosis age heterogeneity adjusted for first-cancer site group", age_adjusted_p, "Poisson likelihood-ratio test", "Adjusted for 18 mutually exclusive SEER first-cancer site groups")
for (lat in c("2-11 months", "12-59 months", "60-119 months", "120+ months")) {
  q <- age_brain$by_latency[age_brain$by_latency$Latency == lat, ]
  add_test(paste0("Diagnosis age heterogeneity, brain-only, ", lat), poisson_heterogeneity(q, "Age_Group"), "Poisson likelihood-ratio test", "Sensitivity analysis restricted to Brain sites")
}
add_test("Sex heterogeneity", poisson_heterogeneity(sex_table, "Sex"), "Poisson likelihood-ratio test")
add_test("Diagnosis-era heterogeneity", poisson_heterogeneity(era_total, "Era"), "Poisson likelihood-ratio test", "Crude across latency distributions")
add_test("Diagnosis-era ordinal trend", poisson_trend(era_total, seq_len(nrow(era_total)) - 1L), "Poisson trend test", "Crude across latency distributions")
for (lat in c("2-11 months", "12-59 months", "60-119 months", "120+ months")) {
  q <- era_latency[era_latency$Latency == lat, ]
  add_test(paste0("Diagnosis-era heterogeneity, ", lat), poisson_heterogeneity(q, "Era"), "Poisson likelihood-ratio test", "Within latency stratum")
  add_test(paste0("Diagnosis-era trend, ", lat), poisson_trend(q, seq_len(nrow(q)) - 1L), "Poisson trend test", "Within latency stratum")
}
race_for_test <- race_total[!race_total[[race_var]] %in% "Unknown", ]
add_test("Race heterogeneity (excluding unknown)", poisson_heterogeneity(race_for_test, race_var), "Poisson likelihood-ratio test", "Descriptive; Other pools AI/AN and API")
add_test("Recorded RT vs None/Unknown", poisson_heterogeneity(rt_binary, "RT_Binary"), "Poisson likelihood-ratio test", "None and unknown cannot be separated")
tests <- do.call(rbind, tests)
write_csv(tests, "table_24_statistical_tests.csv")

# ---- Quality control ----
qc_files <- data.frame(
  Dataset = names(dat),
  File = unname(input_files),
  Rows = vapply(dat, nrow, integer(1)),
  Columns = vapply(dat, ncol, integer(1)),
  stringsAsFactors = FALSE
)

main_overall <- main_summary[main_summary$Outcome == outcome_target & main_summary$Latency == "Total", ]
qc_checks <- data.frame(
  Check = c(
    "Main age strata observed sum matches overall",
    "Main age strata expected sum matches overall",
    "Brain-only age strata observed sum matches Brain-site total",
    "Brain-only age expected sum matches main Brain total",
    "12-month sensitivity events do not exceed primary events",
    "60-month sensitivity events do not exceed 12-month events",
    "Case listing has one unique patient per event row",
    "Case listing event rows match overall observed events",
    "Case listing event-site counts match event rows",
    "All event behaviors are malignant",
    "First-event histology matches event-row histology",
    "Index-age groups in case listing cover all event rows"
  ),
  Observed_Value = c(
    sum(age_main$total$Observed),
    sum(age_main$total$Expected),
    sum(age_brain$total$Observed),
    sum(age_brain$total$Expected),
    sum(age_12m$total$Observed),
    sum(age_60m$total$Observed),
    length(unique(cases$`Patient ID`)),
    nrow(events),
    paste0(sum(events$`Site recode ICD-O-3/WHO 2008 (for SIRs) (Event Variable)` == "Brain"), "+",
           sum(events$`Site recode ICD-O-3/WHO 2008 (for SIRs) (Event Variable)` == "Cranial Nerves Other Nervous System")),
    sum(events$`Behavior code ICD-O-3` == "Malignant"),
    sum(events$`Histologic Type ICD-O-3` == events$`First Event.Histologic Type ICD-O-3`, na.rm = TRUE),
    sum(table(events$Index_Age_Group))
  ),
  Expected_Value = c(
    as.character(main_overall$Observed),
    as.character(main_overall$Expected),
    as.character(sum(main_summary$Observed[main_summary$Outcome == outcome_brain & main_summary$Latency == "Total"])),
    as.character(sum(main_summary$Expected[main_summary$Outcome == outcome_brain & main_summary$Latency == "Total"])),
    paste0("<=", sum(age_main$total$Observed)),
    paste0("<=", sum(age_12m$total$Observed)),
    as.character(nrow(events)),
    as.character(main_overall$Observed),
    as.character(nrow(events)),
    as.character(nrow(events)),
    as.character(nrow(events)),
    as.character(nrow(events))
  ),
  Status = c(
    ifelse(sum(age_main$total$Observed) == main_overall$Observed, "PASS", "FAIL"),
    ifelse(abs(sum(age_main$total$Expected) - main_overall$Expected) < 0.02, "PASS", "FAIL"),
    ifelse(sum(age_brain$total$Observed) == sum(main_summary$Observed[main_summary$Outcome == outcome_brain & main_summary$Latency == "Total"]), "PASS", "FAIL"),
    ifelse(abs(sum(age_brain$total$Expected) - sum(main_summary$Expected[main_summary$Outcome == outcome_brain & main_summary$Latency == "Total"])) < 0.02, "PASS", "FAIL"),
    ifelse(sum(age_12m$total$Observed) <= sum(age_main$total$Observed), "PASS", "FAIL"),
    ifelse(sum(age_60m$total$Observed) <= sum(age_12m$total$Observed), "PASS", "FAIL"),
    ifelse(length(unique(cases$`Patient ID`)) == nrow(events), "PASS", "FAIL"),
    ifelse(nrow(events) == main_overall$Observed, "PASS", "FAIL"),
    ifelse(sum(events$`Site recode ICD-O-3/WHO 2008 (for SIRs) (Event Variable)` %in%
                 c("Brain", "Cranial Nerves Other Nervous System")) == nrow(events), "PASS", "FAIL"),
    ifelse(all(events$`Behavior code ICD-O-3` == "Malignant"), "PASS", "FAIL"),
    ifelse(all(events$`Histologic Type ICD-O-3` == events$`First Event.Histologic Type ICD-O-3`), "PASS", "FAIL"),
    ifelse(sum(table(events$Index_Age_Group)) == nrow(events), "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)
write_csv(qc_files, "qc_01_input_files.csv")
write_csv(qc_checks, "qc_02_validation_checks.csv")
if (any(qc_checks$Status != "PASS")) stop("One or more critical QC checks failed. See qc_02_validation_checks.csv")

# Privacy scan after all tables are written.
derived_csv <- list.files(c(tables_dir, source_dir), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
privacy_scan <- data.frame(File = basename(derived_csv), Contains_Patient_ID = FALSE)
for (k in seq_along(derived_csv)) {
  header <- readLines(derived_csv[k], n = 1L, warn = FALSE)
  privacy_scan$Contains_Patient_ID[k] <- grepl("Patient ID", header, fixed = TRUE)
}
write.csv(privacy_scan, file.path(logs_dir, "privacy_scan.csv"), row.names = FALSE)
if (any(privacy_scan$Contains_Patient_ID)) stop("Privacy QC failed: Patient ID found in derived output.")

# ---- Figure contract and plotting ----
palette <- c(
  navy = "#244A6B", blue = "#4C78A8", teal = "#59A89C",
  orange = "#E39C37", red = "#C44E52", grey = "#7A7A7A", light = "#D9E2E8"
)

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 0.5, hjust = 0),
      strip.text = element_text(face = "bold"),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "#E6E6E6"),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(5, 6, 5, 5)
    )
}

save_pub <- function(plot, basename_out, width_mm = 183, height_mm = 110, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svg_path <- file.path(figures_dir, paste0(basename_out, ".svg"))
  pdf_path <- file.path(figures_dir, paste0(basename_out, ".pdf"))
  tiff_path <- file.path(figures_dir, paste0(basename_out, ".tiff"))
  png_path <- file.path(figures_dir, paste0(basename_out, "_preview.png"))
  svglite::svglite(svg_path, width = w, height = h, bg = "white")
  print(plot); dev.off()
  grDevices::cairo_pdf(pdf_path, width = w, height = h, family = "sans", bg = "white")
  print(plot); dev.off()
  grDevices::tiff(tiff_path, width = w, height = h, units = "in", res = dpi,
                  compression = "lzw", type = "cairo", bg = "white")
  print(plot); dev.off()
  grDevices::png(png_path, width = w, height = h, units = "in", res = 200,
                 type = "cairo", bg = "white")
  print(plot); dev.off()
}

# Figure 1: diagnosis-age forest with sensitivity analyses.
fig1 <- age_sensitivity
fig1$Analysis <- factor(fig1$Analysis,
                        levels = c("2-month exclusion", "12-month exclusion", "60-month exclusion"),
                        labels = c("Primary (2 mo)", "12-month risk set", "60-month risk set"))
fig1$Age_Group <- factor(fig1$Age_Group, levels = rev(c("0-14", "15-19", "20-29", "30-39")))
fig1$Label <- sprintf("%d / %.2f", fig1$Observed, fig1$Expected)
pd <- position_dodge(width = 0.55)
p1 <- ggplot(fig1, aes(x = SIR, y = Age_Group, colour = Analysis, shape = Analysis)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4, colour = palette["grey"]) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), orientation = "y", width = 0,
                linewidth = 0.55, position = pd) +
  geom_point(size = 2.0, stroke = 0.35, position = pd) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10)) +
  coord_cartesian(xlim = c(0.5, 12)) +
  scale_colour_manual(values = unname(c(palette["navy"], palette["teal"], palette["orange"]))) +
  labs(x = "Standardized incidence ratio (log scale)", y = "Age at first cancer (years)",
       title = "Excess brain and nervous system tumor risk is concentrated in younger survivors") +
  theme_pub() + theme(legend.position = "bottom")
save_pub(p1, "Figure_1_age_forest", 183, 95)
write_source(fig1, "Figure_1_source_data.csv")

# Figure 2: latency pattern by diagnosis age.
fig2 <- age_main$by_latency
fig2$Age_Group <- factor(fig2$Age_Group, levels = c("0-14", "15-19", "20-29", "30-39"))
fig2$Latency <- factor(fig2$Latency, levels = c("2-11 months", "12-59 months", "60-119 months", "120+ months"),
                       labels = c("2-11 mo", "1-4 y", "5-9 y", "10+ y"))
fig2$SIR_plot <- ifelse(fig2$Observed == 0, NA_real_, fig2$SIR)
fig2$Lower_plot <- ifelse(fig2$Observed == 0, NA_real_, fig2$CI_Lower)
fig2$Upper_plot <- ifelse(fig2$Observed == 0, NA_real_, fig2$CI_Upper)
fig2_nonzero <- fig2[fig2$Observed > 0, ]
pd2 <- position_dodge(width = 0.42)
p2 <- ggplot(fig2_nonzero, aes(x = Latency, y = SIR_plot, colour = Age_Group)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4, colour = palette["grey"]) +
  geom_errorbar(aes(ymin = Lower_plot, ymax = Upper_plot), width = 0.08, linewidth = 0.45,
                position = pd2) +
  geom_point(size = 1.8, position = pd2) +
  geom_text(aes(label = paste0("n=", Observed)), vjust = -0.85, size = 2.5,
            position = pd2, show.legend = FALSE) +
  geom_text(data = fig2[fig2$Observed == 0, ], aes(y = 0.42, label = "0 events"),
            inherit.aes = TRUE, size = 2.1, show.legend = FALSE) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
  coord_cartesian(ylim = c(0.35, 25)) +
  scale_colour_manual(values = unname(c(palette["red"], palette["orange"], palette["teal"], palette["navy"]))) +
  labs(x = "Latency since first cancer", y = "Standardized incidence ratio (log scale)",
       title = "Relative incidence by age at first cancer and latency") +
  theme_pub(8) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
save_pub(p2, "Figure_2_age_latency", 183, 105)
write_source(fig2, "Figure_2_source_data.csv")

# Figure 3: child and AYA first-cancer categories.
child_plot <- iccc_total[is.finite(iccc_total$SIR) & iccc_total$Expected > 0 & iccc_total$Observed > 0, ]
child_plot$Label <- sub("^(I|II|IV|V|VI|VII|VIII|IX|X|XI|XII) ", "", child_plot$Primary_Cancer_ICCC)
child_plot$Label <- tools::toTitleCase(tolower(child_plot$Label))
child_plot <- child_plot[order(child_plot$SIR), ]
child_plot$Label <- factor(child_plot$Label, levels = child_plot$Label)
p3a <- ggplot(child_plot, aes(SIR, Label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = palette["grey"]) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), orientation = "y", width = 0,
                colour = palette["navy"], linewidth = 0.45) +
  geom_point(aes(size = Observed), colour = palette["navy"], alpha = 0.9) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50)) +
  coord_cartesian(xlim = c(0.4, 50)) +
  scale_size_continuous(range = c(1.2, 3.2)) +
  labs(x = "SIR (log scale)", y = NULL, title = "Childhood (0-14 y)") +
  theme_pub(7.2) + theme(legend.position = "none")

aya_plot <- aya_total[is.finite(aya_total$SIR) & aya_total$Expected > 0 & aya_total$Observed > 0, ]
aya_plot$Label <- sub("^[0-9]+\\. ", "", aya_plot$Primary_Cancer_AYA)
aya_plot <- aya_plot[order(aya_plot$SIR), ]
aya_plot$Label <- factor(aya_plot$Label, levels = aya_plot$Label)
p3b <- ggplot(aya_plot, aes(SIR, Label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = palette["grey"]) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), orientation = "y", width = 0,
                colour = palette["teal"], linewidth = 0.45) +
  geom_point(aes(size = Observed), colour = palette["teal"], alpha = 0.9) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50)) +
  coord_cartesian(xlim = c(0.4, 50)) +
  scale_size_continuous(range = c(1.2, 3.2)) +
  labs(x = "SIR (log scale)", y = NULL, title = "AYA (15-39 y)") +
  theme_pub(7.2) + theme(legend.position = "none")

p3 <- p3a | p3b
p3 <- p3 + plot_annotation(title = "Risk varies by type of first cancer", tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8),
        plot.title = element_text(size = 8, face = "bold"))
save_pub(p3, "Figure_3_primary_cancer_forest", 183, 175)
write_source(child_plot, "Figure_3a_child_source_data.csv")
write_source(aya_plot, "Figure_3b_aya_source_data.csv")

# Figure 4: histology composition of subsequent tumors.
fig4 <- hist_group
fig4 <- fig4[order(fig4$N), ]
fig4$Histology_Group <- factor(fig4$Histology_Group, levels = fig4$Histology_Group)
p4 <- ggplot(fig4, aes(N, Histology_Group)) +
  geom_col(width = 0.68, fill = palette["navy"]) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", N, Percent)), hjust = -0.08, size = 2.3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Number of subsequent malignant brain and nervous system tumors", y = NULL,
       title = "Gliomas dominate the histologic spectrum") +
  theme_pub() + theme(panel.grid.major.x = element_blank())
save_pub(p4, "Figure_4_event_histology", 150, 105)
write_source(fig4, "Figure_4_source_data.csv")

# Figure/table notes and reproducibility log.
figure_contract <- c(
  "Core conclusion: Excess subsequent malignant brain and other nervous system tumor risk is concentrated in survivors diagnosed with cancer at younger ages and persists long term after childhood cancer.",
  "Figure archetype: quantitative grid with a diagnosis-age forest as the hero figure.",
  "Backend: R only (ggplot2, patchwork, svglite, Cairo devices).",
  "Figure 1: primary age-gradient evidence with 12-month and 5-year-survivor sensitivity analyses.",
  "Figure 2: latency-specific validation of persistent childhood risk.",
  "Figure 3: first-cancer category heterogeneity in child and AYA cohorts.",
  "Figure 4: descriptive histology composition of observed subsequent tumors.",
  "Statistics: exact Poisson 95% confidence intervals; Poisson likelihood-ratio heterogeneity tests.",
  "Review risks: sparse subgroup counts; noncausal race and radiation comparisons; unequal follow-up by diagnosis era; historical ICD-O classification changes.",
  "Privacy: all source-data files are aggregate and contain no Patient ID."
)
writeLines(figure_contract, file.path(results_dir, "figure_contract.txt"))

qa_notes <- c(
  paste("Pipeline completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R version:", R.version.string),
  paste("Input files:", length(input_files)),
  paste("Critical QC checks passed:", sum(qc_checks$Status == "PASS"), "of", nrow(qc_checks)),
  paste("Derived CSV files scanned for Patient ID:", nrow(privacy_scan)),
  "All plotting and graphics exports were generated in R.",
  "SVG/PDF retain editable text; TIFF exports use 600 dpi LZW compression.",
  "No individual-level case listing was copied to the results directory."
)
writeLines(qa_notes, file.path(logs_dir, "pipeline_qa_notes.txt"))
capture.output(sessionInfo(), file = file.path(logs_dir, "R_session_info.txt"))

cat("Pipeline completed successfully.\n")
cat("Tables:", length(list.files(tables_dir, pattern = "\\.csv$")), "\n")
cat("Figures:", length(list.files(figures_dir)), "files\n")
cat("Source data:", length(list.files(source_dir, pattern = "\\.csv$")), "\n")
