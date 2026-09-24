options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
root <- normalizePath(Sys.getenv("SEER_PROJECT_DIR", unset = dirname(dirname(script_path))), winslash = "/", mustWork = TRUE)

raw_file <- file.path(root, "02_raw_exports", "15_full_cohort_index_competing_risk.txt")
boundary_file <- file.path(root, "02_raw_exports", "16_full_cohort_all_malignancies_0m.txt")
out_dir <- Sys.getenv("SEER_CIF_OUTPUT_DIR", unset = file.path(root, "04_results", "competing_risk_20260924"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- fread(raw_file, sep = "\t", na.strings = "NA", quote = "\"")
stopifnot(uniqueN(d$`Patient ID`) == 276844L)

idx <- d[`Event Number` == "Index Record", .(
  patient_id = `Patient ID`,
  age_raw = `Age recode with <1 year olds and 90+`,
  index_survival_months_raw = `Survival months`
)]
stopifnot(nrow(idx) == 276844L, uniqueN(idx$patient_id) == 276844L)

ex <- d[trimws(`Reason for Exit`) != "", .(
  patient_id = `Patient ID`,
  reason = trimws(`Reason for Exit`),
  event_site = `Site recode ICD-O-3/WHO 2008 (for SIRs) (Event Variable)`,
  time_years = as.numeric(`Person Time Years (Calculated)`)
)]
eligible_target <- d[`Event Number` == "1"]
stopifnot(nrow(eligible_target) == 345L,
          all(as.numeric(eligible_target$`Months Since Index (Calculated)`) >= 2),
          sum(as.numeric(eligible_target$`Months Since Index (Calculated)`) == 2) == 1L)

# When available, a separate otherwise-matched 0-month case listing verifies
# the 2-month eligibility boundary independently of the primary export.
boundary_qa <- NULL
if (file.exists(boundary_file)) {
  b <- fread(boundary_file, sep = "\t", na.strings = "NA", quote = "\"")
  b_target <- b[`Event Number` == "1" &
                  `Site recode ICD-O-3/WHO 2008 (for SIRs) (Event Variable)` %chin%
                    c("Brain", "Cranial Nerves Other Nervous System")]
  b_target[, latency_months := as.numeric(`Months Since Index (Calculated)`)]
  stopifnot(nrow(b_target) == 366L,
            sum(b_target$latency_months < 2) == 21L,
            sum(b_target$latency_months == 2) == 1L,
            sum(b_target$latency_months > 2) == 344L,
            setequal(eligible_target$`Patient ID`,
                     b_target[latency_months >= 2, `Patient ID`]))
  boundary_qa <- list(
    zero_month_people = uniqueN(b$`Patient ID`),
    target_before_landmark = sum(b_target$latency_months < 2),
    target_at_landmark = sum(b_target$latency_months == 2),
    target_after_landmark = sum(b_target$latency_months > 2)
  )
}
stopifnot(nrow(ex) == 276844L, uniqueN(ex$patient_id) == 276844L, !anyNA(ex$time_years))

cohort <- merge(idx, ex, by = "patient_id", all = TRUE)
cohort[, index_survival_months := suppressWarnings(as.numeric(index_survival_months_raw))]
cohort[, age_num := fifelse(age_raw == "00 years", 0,
                     fifelse(grepl("^01-04", age_raw), 1,
                     fifelse(grepl("^05-09", age_raw), 5,
                     fifelse(grepl("^10-14", age_raw), 10,
                     fifelse(grepl("^15-19", age_raw), 15,
                     fifelse(grepl("^20-24", age_raw), 20,
                     fifelse(grepl("^25-29", age_raw), 25,
                     fifelse(grepl("^30-34", age_raw), 30,
                     fifelse(grepl("^35-39", age_raw), 35, NA_real_)))))))))]
stopifnot(!anyNA(cohort$age_num))
cohort[, age_group := factor(fcase(
  age_num <= 14, "0-14 years",
  age_num <= 19, "15-19 years",
  age_num <= 29, "20-29 years",
  default = "30-39 years"
), levels = c("0-14 years", "15-19 years", "20-29 years", "30-39 years"))]

target_sites <- c("Brain", "Cranial Nerves Other Nervous System")
cohort[, status_chr := fcase(
  reason == "Any Event in Rate File" & event_site %chin% target_sites, "target_cns",
  reason == "Any Event in Rate File", "other_malignancy",
  reason == "Death", "death",
  default = "censor"
)]

# The MP-SIR clock starts after the 2-month exclusion. Subjects known to have
# died or been lost before that landmark are not members of the landmark risk set.
cohort[, pre_entry := time_years == 0 & reason %chin% c("Death", "Lost to Follow-up") &
                         !is.na(index_survival_months) & index_survival_months < 2]
analysis <- cohort[pre_entry == FALSE]

# SEER monthly dates can place eligible events or exits exactly at the landmark.
# A tiny positive time retains these observations in software requiring time > 0.
analysis[, time_analysis := fifelse(time_years == 0, 1e-8, time_years)]
analysis[, status := factor(status_chr,
  levels = c("censor", "target_cns", "other_malignancy", "death"))]

qa <- list(
  raw_rows = nrow(d), raw_people = uniqueN(d$`Patient ID`),
  exit_rows = nrow(ex), excluded_before_landmark = sum(cohort$pre_entry),
  analysis_people = nrow(analysis),
  status_counts = table(analysis$status_chr),
  age_counts = table(analysis$age_group),
  target_counts = table(analysis$age_group[analysis$status_chr == "target_cns"]),
  zero_time_retained = sum(analysis$time_years == 0),
  boundary_validation = boundary_qa
)
capture.output(str(qa), file = file.path(out_dir, "cif_qa.txt"))

fit <- survfit(Surv(time_analysis, status) ~ age_group, data = analysis,
               conf.type = "log-log")

extract_curve <- function(fit) {
  target_col <- match("target_cns", colnames(fit$pstate))
  stopifnot(!is.na(target_col))
  lens <- as.integer(fit$strata)
  groups <- sub("age_group=", "", names(fit$strata), fixed = TRUE)
  starts <- cumsum(c(1L, head(lens, -1L)))
  rbindlist(lapply(seq_along(lens), function(i) {
    z <- starts[i]:(starts[i] + lens[i] - 1L)
    data.table(
      age_group = groups[i], time_years = fit$time[z],
      estimate = fit$pstate[z, target_col],
      se = fit$std.err[z, target_col]
    )
  }))
}
curve <- extract_curve(fit)
curve[, `:=`(lower = NA_real_, upper = NA_real_)]
curve[estimate == 0, `:=`(lower = 0, upper = 0)]
curve[estimate > 0 & estimate < 1, c("lower", "upper") := {
  g <- log(-log(1 - estimate))
  seg <- se / ((1 - estimate) * (-log(1 - estimate)))
  list(1 - exp(-exp(g - qnorm(.975) * seg)),
       1 - exp(-exp(g + qnorm(.975) * seg)))
}]

point_estimate <- function(g, tt) {
  z <- curve[age_group == g & time_years <= tt]
  if (!nrow(z)) return(data.table(estimate = 0, se = 0, lower = 0, upper = 0))
  z[.N, .(estimate, se, lower, upper)]
}
milestones <- rbindlist(lapply(levels(analysis$age_group), function(g)
  rbindlist(lapply(c(5, 10, 20), function(tt)
    cbind(data.table(age_group = g, years = tt), point_estimate(g, tt))))))

# Independent Aalen-Johansen recursion validates the survival package estimates.
manual_aj <- function(x, tt) {
  times <- sort(unique(x$time_analysis[x$time_analysis <= tt]))
  s <- 1; f <- 0
  for (u in times) {
    n <- sum(x$time_analysis >= u)
    dc <- sum(x$time_analysis == u & x$status_chr == "target_cns")
    da <- sum(x$time_analysis == u & x$status_chr != "censor")
    f <- f + s * dc / n
    s <- s * (1 - da / n)
  }
  f
}
manual <- rbindlist(lapply(levels(analysis$age_group), function(g)
  rbindlist(lapply(c(5,10,20), function(tt)
    data.table(age_group=g, years=tt,
               manual_estimate=manual_aj(analysis[age_group==g],tt))))))
milestones <- merge(milestones, manual, by=c("age_group","years"))
stopifnot(max(abs(milestones$estimate-milestones$manual_estimate)) < 1e-12)

risk_counts <- rbindlist(lapply(levels(analysis$age_group), function(g) {
  x <- analysis[age_group == g]
  data.table(age_group = g, years = c(0, 5, 10, 20),
             n_at_risk = sapply(c(0, 5, 10, 20), function(tt) sum(x$time_analysis >= tt)))
}))

# Patient identifiers are used in memory only; no patient-level file is written.
fwrite(curve, file.path(out_dir, "cif_curve_source_data.csv"))
fwrite(milestones, file.path(out_dir, "cif_milestones.csv"))
fwrite(risk_counts, file.path(out_dir, "cif_numbers_at_risk.csv"))

cols <- c("0-14 years" = "#C44E52", "15-19 years" = "#DD8452",
          "20-29 years" = "#4C72B0", "30-39 years" = "#55A868")
plot_curve <- curve[time_years <= 30]
p <- ggplot(plot_curve, aes(time_years, estimate * 100, colour = age_group)) +
  geom_step(linewidth = 0.8) +
  scale_colour_manual(values = cols, name = "Age at first cancer") +
  scale_x_continuous("Years after the 2-month landmark", limits = c(0, 30),
                     breaks = seq(0, 30, 5), expand = expansion(mult = c(0, .01))) +
  scale_y_continuous("Cumulative incidence (%)", limits = c(0, NA),
                     expand = expansion(mult = c(0, .06))) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .35),
        axis.ticks = element_line(linewidth = .35),
        legend.position = c(.19, .78), legend.background = element_blank(),
        legend.key.width = unit(9, "mm"),
        plot.margin = margin(6, 8, 5, 6))

ggsave(file.path(out_dir, "Supplementary_Fig_S2_CIF.svg"), p,
       width = 89, height = 72, units = "mm", device = svglite::svglite)
ggsave(file.path(out_dir, "Supplementary_Fig_S2_CIF.pdf"), p,
       width = 89, height = 72, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "Supplementary_Fig_S2_CIF.tiff"), p,
       width = 89, height = 72, units = "mm", dpi = 600,
       compression = "lzw", device = "tiff")
ggsave(file.path(out_dir, "Supplementary_Fig_S2_CIF_preview.png"), p,
       width = 89, height = 72, units = "mm", dpi = 200)

cat("Analysis complete\n")
print(qa)
print(milestones)
