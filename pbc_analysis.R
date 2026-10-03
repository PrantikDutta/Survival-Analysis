# ============================================================
# SURVIVAL ANALYSIS - PBC DATASET
# ONLY survival PACKAGE REQUIRED
# ============================================================

# Install ONLY if needed
# install.packages("survival")

library(survival)

# Load PBC dataset
data("pbc", package = "survival")

# Check data
head(pbc)
dim(pbc)
summary(pbc)


# ------------------------------------------------------------
# 1. SELECT VARIABLES
# ------------------------------------------------------------

data <- pbc[, c(
  "time",
  "status",
  "trt",
  "age",
  "sex",
  "bili",
  "albumin",
  "edema",
  "stage"
)]


# ------------------------------------------------------------
# 2. CREATE EVENT VARIABLE
# ------------------------------------------------------------
# Original:
# 0 = censored
# 1 = transplant
# 2 = death
#
# Here:
# 0 = censored
# 1 = death

data$status_death <- ifelse(data$status == 2, 1, 0)


# ------------------------------------------------------------
# 3. REMOVE MISSING VALUES
# ------------------------------------------------------------

data <- na.omit(data)

dim(data)

table(data$status_death)

prop.table(table(data$status_death)) * 100


# ------------------------------------------------------------
# 4. CONVERT CATEGORICAL VARIABLES
# ------------------------------------------------------------

data$trt <- factor(
  data$trt,
  levels = c(1, 2),
  labels = c("Treatment 1", "Treatment 2")
)

data$sex <- factor(data$sex)

data$stage <- factor(data$stage)


# ============================================================
# 5. BASIC EXPLORATORY ANALYSIS
# ============================================================

summary(data)
par(mfrow=c(1,3))
hist(
  data$age,
  main = "Distribution of Age",
  xlab = "Age",
  breaks = 20
)

hist(
  data$bili,
  main = "Distribution of Bilirubin",
  xlab = "Bilirubin",
  breaks = 20
)

hist(
  data$albumin,
  main = "Distribution of Albumin",
  xlab = "Albumin",
  breaks = 20
)


# ============================================================
# 6. SURVIVAL OBJECT
# ============================================================

surv_object <- Surv(
  data$time,
  data$status_death
)

head(surv_object)


# ============================================================
# 7. OVERALL KAPLAN-MEIER ANALYSIS
# ============================================================

km_fit <- survfit(
  surv_object ~ 1,
  data = data
)

summary(km_fit)

# Median survival
km_fit$table

# Kaplan-Meier plot
plot(
  km_fit,
  main = "Kaplan-Meier Survival Curve",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  conf.int = TRUE,
  mark.time = TRUE
)

grid()


# ============================================================
# 8. SURVIVAL BY TREATMENT
# ============================================================

km_trt <- survfit(
  Surv(time, status_death) ~ trt,
  data = data
)

summary(km_trt)

plot(
  km_trt,
  col = 1:2,
  lwd = 2,
  main = "Kaplan-Meier Survival by Treatment",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  conf.int = TRUE
)

legend(
  "topright",
  legend = levels(data$trt),
  col = 1:2,
  lwd = 2
)

grid()


# ------------------------------------------------------------
# 9. LOG-RANK TEST: TREATMENT
# ------------------------------------------------------------

logrank_trt <- survdiff(
  Surv(time, status_death) ~ trt,
  data = data
)

print(logrank_trt)

p_trt <- 1 - pchisq(
  logrank_trt$chisq,
  df = length(logrank_trt$n) - 1
)

cat("Treatment log-rank p-value =", p_trt, "\n")


# ============================================================
# 10. SURVIVAL BY SEX
# ============================================================

km_sex <- survfit(
  Surv(time, status_death) ~ sex,
  data = data
)

plot(
  km_sex,
  col = 1:length(levels(data$sex)),
  lwd = 2,
  main = "Kaplan-Meier Survival by Sex",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  conf.int = TRUE
)

legend(
  "topright",
  legend = levels(data$sex),
  col = 1:length(levels(data$sex)),
  lwd = 2
)

grid()


# ------------------------------------------------------------
# 11. LOG-RANK TEST: SEX
# ------------------------------------------------------------

logrank_sex <- survdiff(
  Surv(time, status_death) ~ sex,
  data = data
)

print(logrank_sex)

p_sex <- 1 - pchisq(
  logrank_sex$chisq,
  df = length(logrank_sex$n) - 1
)

cat("Sex log-rank p-value =", p_sex, "\n")


# ============================================================
# 12. SURVIVAL BY DISEASE STAGE
# ============================================================

km_stage <- survfit(
  Surv(time, status_death) ~ stage,
  data = data
)

plot(
  km_stage,
  col = 1:length(levels(data$stage)),
  lwd = 2,
  main = "Kaplan-Meier Survival by Disease Stage",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  conf.int = FALSE
)

legend(
  "topright",
  legend = levels(data$stage),
  col = 1:length(levels(data$stage)),
  lwd = 2
)

grid()


# ------------------------------------------------------------
# 13. LOG-RANK TEST: STAGE
# ------------------------------------------------------------

logrank_stage <- survdiff(
  Surv(time, status_death) ~ stage,
  data = data
)

print(logrank_stage)

p_stage <- 1 - pchisq(
  logrank_stage$chisq,
  df = length(logrank_stage$n) - 1
)

cat("Stage log-rank p-value =", p_stage, "\n")


# ============================================================
# 14. COX PROPORTIONAL HAZARDS MODEL
# ============================================================

cox_model <- coxph(
  Surv(time, status_death) ~
    age + sex + trt + bili + albumin + edema + stage,
  data = data
)

summary(cox_model)


# ============================================================
# 15. HAZARD RATIOS
# ============================================================

HR <- exp(coef(cox_model))

CI <- exp(confint(cox_model))

p_values <- summary(cox_model)$coefficients[, 5]

cox_results <- data.frame(
  Variable = names(HR),
  Hazard_Ratio = HR,
  Lower_95_CI = CI[, 1],
  Upper_95_CI = CI[, 2],
  P_Value = p_values
)

print(cox_results)


# ============================================================
# 16. COX MODEL RESULTS
# ============================================================

print(cox_results)

# Plot estimated hazard ratios with 95% confidence intervals

HR <- cox_results$Hazard_Ratio
lower <- cox_results$Lower_95_CI
upper <- cox_results$Upper_95_CI

plot(
  HR,
  seq_along(HR),
  xlim = range(c(lower, upper)),
  yaxt = "n",
  xlab = "Hazard Ratio",
  ylab = "",
  main = "Cox Proportional Hazards Model",
  pch = 19
)

segments(
  lower,
  seq_along(HR),
  upper,
  seq_along(HR)
)

abline(
  v = 1,
  lty = 2
)

axis(
  2,
  at = seq_along(HR),
  labels = cox_results$Variable,
  las = 1
)

# ============================================================
# 17. PROPORTIONAL HAZARDS ASSUMPTION
# ============================================================

ph_test <- cox.zph(cox_model)

par(mfrow=c(2,4))

plot(
  ph_test,
  main="Schoenfeld Residual Test"
)

par(mfrow=c(1,1))

# ============================================================
# 18. FINAL SUMMARY
# ============================================================

cat("\n====================================\n")
cat("FINAL SURVIVAL ANALYSIS SUMMARY\n")
cat("====================================\n")

cat("Number of patients:", nrow(data), "\n")

cat(
  "Number of deaths:",
  sum(data$status_death == 1),
  "\n"
)

cat(
  "Number censored:",
  sum(data$status_death == 0),
  "\n"
)

cat(
  "Treatment log-rank p-value:",
  p_trt,
  "\n"
)

cat(
  "Sex log-rank p-value:",
  p_sex,
  "\n"
)

cat(
  "Stage log-rank p-value:",
  p_stage,
  "\n"
)

cat("\nCox Regression Results:\n")

print(cox_results)

cat("\nProportional Hazards Test:\n")

print(ph_test)
