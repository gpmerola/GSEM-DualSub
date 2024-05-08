library(corrplot)
library(ggplot2)
library(devtools)
library(shru)
library(forcats)
library(tidyr)
library(dplyr)
library(Matrix)

setwd(file.path("~", "proj", "NewGwasSub"))

dir_path <- readLines("dir.txt", warn = FALSE)
if (length(dir_path) > 0 && dir.exists(dir_path[1])) {
  setwd(dir_path[1])
} else {
  stop("Directory path in dir.txt is invalid or does not exist.")
}

add_significance <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}


difftest.matrix <- function(
    mValues1,
    mStandard_errors1,
    mValues2,
    mStandard_errors2,
    df1 = NULL,
    df2 = NULL,
    mCorrelationEstimate.values = 0.95,
    mCorrelationEstimate.standard_errors = 0.95,
    effectiveNumberOfTests=NULL,
    mValueCovariances=NULL,
    symmetric = T,
    eigenSumLimit = 0.995
){
  
  mValues1 <- symMatrixConform(mValues1)
  mStandard_errors1 <- symMatrixConform(mStandard_errors1)
  
  mValues2 <- symMatrixConform(mValues2)
  mStandard_errors2 <- symMatrixConform(mStandard_errors2)
  
  if(is.null(effectiveNumberOfTests) & !is.null(mValueCovariances)){
    effectiveNumberOfTests <- getEffectiveNumberOfTests(covarianceMatrix = mValueCovariances, symmetric = symmetric, eigenSumLimit = eigenSumLimit)
  }
  
  mTest.values <- abs(mValues1 - mValues2)
  sde <- sqrt((mStandard_errors1^2) + (mStandard_errors2^2) - 2 * mCorrelationEstimate.values * mStandard_errors1 * mStandard_errors2)
  pTest.values <- 2 * pnorm(q = mTest.values, sd = sde, lower.tail = F)
  pTest.values.adj <- matrix(data = p.adjust(pTest.values, method = "fdr", n = effectiveNumberOfTests), nrow = nrow(pTest.values), ncol = ncol(pTest.values))
  rownames(pTest.values.adj) <- rownames(pTest.values)
  colnames(pTest.values.adj) <- colnames(pTest.values)
  
  pTest.standard_errors <- NULL
  pTest.standard_errors.adj <- NULL
  
  if(!is.null(mStandard_errors1) & !is.null(mStandard_errors2) & !is.null(df1)){
    if(is.null(df2)) df2 <- df1
    var_bar <- (mStandard_errors1^2 + mStandard_errors2^2) / 2
    tVar <- (mStandard_errors1^2) - (mStandard_errors2^2)
    vare <- var_bar^2 * (2 / (df1 - 1) + 2 / (df2 - 1) - 2 * mCorrelationEstimate.standard_errors * sqrt(2 / (df1 - 1)) * sqrt(2 / (df2 - 1)))
    pTest.standard_errors <- 2 * pnorm(q = abs(tVar), mean = 0, sd = sqrt(vare), lower.tail = F)
    pTest.standard_errors.adj <- matrix(data = p.adjust(pTest.standard_errors, method = "fdr", n = effectiveNumberOfTests), nrow = nrow(pTest.standard_errors), ncol = ncol(pTest.standard_errors))
    rownames(pTest.standard_errors.adj) <- rownames(pTest.standard_errors)
    colnames(pTest.standard_errors.adj) <- colnames(pTest.standard_errors)
  }
  
  return(
    list(
      pTest.values = pTest.values,
      pTest.values.adj = pTest.values.adj,
      pTest.standard_errors = pTest.standard_errors,
      pTest.standard_errors.adj = pTest.standard_errors.adj
    )
  )
}

trait_names <- c("Alcohol Use", "Neuroticism", "Risky Behaviour (speed)", "Risky Behaviour (sex)", "Smoking", "Subjective Well-Being",
                 "BMI", "Chronotype", "Sleep Duration", "Social Deprivation", "Tiredness", "Income", "Physical Activity", "Autism",
                 "Bipolar Disorder", "Cannabis Use", "Schizophrenia", "ADHD", "Diabetes type 2", "GAD", "Anorexia Nervosa",
                 "Alcohol Dependence Disorder", "BD I", "BD II", "Migraine", "ALS", "MDD", "Mania")

trait_labels <- setNames(trait_names, c("ALCO04", "NEUR04", "RISK02", "RISK03", "SMOK10", "SUBJ01", "BODY13B",
                                        "CHRO01", "CHRO02", "INCO03", "TIRE01", "INCO02", "PHYS01", "AUTI09",
                                        "BIPO03", "CANU03", "SCHI06", "ADHD06", "DIAB06", "ANXI06", "ANOR21",
                                        "ALCD01", "BIPO04", "BIPO05", "MIGR01", "NPAT12", "DEPR14", "Mania"))


mega_correlation <- readRDS("MegaCorrelation.rds")

rownames(mega_correlation$S_Stand) <- trait_names
colnames(mega_correlation$S_Stand) <- trait_names

# Create correlation matrix
S <- mega_correlation$S_Stand
R <- cor(S)

# Create standard error matrix
SE_R <- matrix(nrow = nrow(S), ncol = ncol(S))
rownames(SE_R) <- trait_names
colnames(SE_R) <- trait_names

SE_R[lower.tri(SE_R, diag = TRUE)] <- sqrt(diag(mega_correlation$V_Stand))
SE_R[upper.tri(SE_R, diag = TRUE)] <- t(SE_R)[upper.tri(SE_R, diag = TRUE)]

# Extract correlation vectors and standard errors for Mania and Bipolar Disorder
R_mania <- R["Mania", ]
R_bd <- R["Bipolar Disorder", ]
seR_mania <- SE_R["Mania", ]
seR_bd <- SE_R["Bipolar Disorder", ]



p_overall <- difftest.matrix(effectiveNumberOfTests = 28,
                               mValues1 = as.matrix(R_mania), 
                               mStandard_errors1 = as.matrix(seR_mania), 
                               mValues2 = as.matrix(R_bd), 
                               mStandard_errors2 = as.matrix(seR_bd), mCorrelationEstimate.values = 0.897847609)



p_adjusted <- p_overall$pTest.values.adj[!is.na(p_overall$pTest.values.adj)]

# Pair adjusted p-values with trait names
paired_pval <- cbind(trait_names, p_adjusted)

# Print the two-column vector
print(paired_pval)



R_mania <- data.frame(Trait = names(R_mania), Correlation = R_mania)
R_bd <- data.frame(Trait = names(R_bd), Correlation = R_bd)
seR_mania <- data.frame(Trait = names(seR_mania), Correlation = seR_mania)
seR_bd <- data.frame(Trait = names(seR_bd), Correlation = seR_bd)

dfs <- list(R_mania, R_bd, seR_mania, seR_bd)



merged_df <- Reduce(function(x, y) merge(x, y, by = "Trait"), dfs)
new_colnames <- c("Trait", "Correlation_Mania", "Correlation_BD", "SE_Mania", "SE_BD")

names(merged_df) <- new_colnames

print(paired_pval)
paired_pval[, "p_adjusted"] <- as.numeric(paired_pval[, "p_adjusted"])

# Merge p_adjusted to merged_df
merged_df$p_adjusted <- paired_pval[match(merged_df$Trait, paired_pval[, "trait_names"]), "p_adjusted"]

custom_order <- c(
  "Bipolar Disorder",
  "BD I",
  "BD II",
  "Schizophrenia",
  "MDD",
  "Neuroticism",
  "Autism",
  "ADHD",
  "GAD",
  "Anorexia Nervosa",
  "Alcohol Use",
  "Alcohol Dependence Disorder",
  "Risky Behaviour (speed)",
  "Risky Behaviour (sex)",
  "Smoking",
  "Cannabis Use",
  "Subjective Well-Being",
  "Income",
  "Social Deprivation",
  "BMI",
  "Chronotype",
  "Sleep Duration",
  "Tiredness",
  "Physical Activity",
  "Diabetes type 2",
  "Migraine",
  "ALS",
  "Mania"
)
merged_df <- merged_df[match(custom_order, merged_df$Trait), ]
merged_df$p_adjusted <- as.numeric(merged_df$p_adjusted)
merged_df$significance <- sapply(merged_df$p_adjusted, add_significance)
merged_df$label_sign <- paste(merged_df$Trait, merged_df$significance)
merged_df <- merged_df[nrow(merged_df):1, ]

merged_df <- merged_df[, -1]
colnames(merged_df)[ncol(merged_df)] <- "Trait"
merged_df <- merged_df[, c(ncol(merged_df), 1:(ncol(merged_df)-1))]
merged_df$p_adjusted <- NULL
merged_df$significance <- NULL
merged_df <- merged_df[-1, ]

print(merged_df)
write.csv(merged_df, "trait_correlation_data.csv", row.names = FALSE)
