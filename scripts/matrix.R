library(corrplot)
library(ggplot2)
library(dplyr)
library(Matrix)
library(shru)

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

difftest.matrix <- function(mValues1, mStandard_errors1, mValues2, mStandard_errors2, mCorrelationEstimate.values = 0.95, effectiveNumberOfTests=NULL, mValueCovariances=NULL, symmetric = T, eigenSumLimit = 0.995) {
  
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
  
  return(list(pTest.values = pTest.values, pTest.values.adj = pTest.values.adj))
}

trait_names <- c(read.csv("corr.csv")$trait, "Mania")
cluster <- c(read.csv("corr.csv")$Cluster)

mega_correlation <- readRDS("MegaCorrelation.rds")

rownames(mega_correlation$S_Stand) <- trait_names
colnames(mega_correlation$S_Stand) <- trait_names

S <- mega_correlation$S_Stand

SE_R <- matrix(nrow = nrow(S), ncol = ncol(S))
rownames(SE_R) <- trait_names
colnames(SE_R) <- trait_names

SE_R[lower.tri(SE_R, diag = TRUE)] <- sqrt(diag(mega_correlation$V_Stand))
SE_R[upper.tri(SE_R, diag = TRUE)] <- t(SE_R)[upper.tri(SE_R, diag = TRUE)]

R_mania <- S["Mania", ]
R_bd <- S["Bipolar Disorder", ]
seR_mania <- SE_R["Mania", ]
seR_bd <- SE_R["Bipolar Disorder", ]

p_overall <- difftest.matrix(effectiveNumberOfTests = 35, mValues1 = as.matrix(R_mania), mStandard_errors1 = as.matrix(seR_mania), mValues2 = as.matrix(R_bd), mStandard_errors2 = as.matrix(seR_bd), mCorrelationEstimate.values = 0.897847609)

p_adjusted <- p_overall$pTest.values.adj[!is.na(p_overall$pTest.values.adj)]
paired_pval <- cbind(trait_names, p_adjusted)

R_mania <- data.frame(Trait = names(R_mania), Correlation = R_mania)
R_bd <- data.frame(Trait = names(R_bd), Correlation = R_bd)
seR_mania <- data.frame(Trait = names(seR_mania), Correlation = seR_mania)
seR_bd <- data.frame(Trait = names(seR_bd), Correlation = seR_bd)

dfs <- list(R_mania, R_bd, seR_mania, seR_bd)

merged_df <- Reduce(function(x, y) merge(x, y, by = "Trait"), dfs)
new_colnames <- c("Trait", "Correlation_Mania", "Correlation_BD", "SE_Mania", "SE_BD")
names(merged_df) <- new_colnames

paired_pval[, "p_adjusted"] <- as.numeric(paired_pval[, "p_adjusted"])
merged_df$p_adjusted <- paired_pval[match(merged_df$Trait, paired_pval[, "trait_names"]), "p_adjusted"]

merged_df_with_pvalues <- merged_df

custom_order <- c("Bipolar Disorder", "BD I", "BD II", "Schizophrenia", "MDD", "Neuroticism", "Autism", "ADHD", "GAD", "Anorexia Nervosa", "Alcohol Use", "Alcohol Dependence Disorder", "Risky Behaviour (speed)", "Risky Behaviour (sex)", "Smoking", "Cannabis Use", "Subjective Well-Being", "Income", "IQ", "Educational Attainment", "Educational Attainment (Non-Cog)", "Educational Attainment (Cog)", "Social Deprivation", "BMI", "Chronotype", "Sleep Duration", "Tiredness", "Physical Activity", "Diabetes type 2", "Migraine", "ALS", "Alzheimer's Disease", "Parkinson's Disease", "Epilepsy", "Mania")
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

print(merged_df_with_pvalues)

merged_df$Cluster <- rev(cluster)

write.csv(merged_df, "trait_correlation_data.csv", row.names = FALSE) 
