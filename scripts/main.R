print("---------------------------------------------------------CodeStart - SubGwas")

args <- commandArgs(trailingOnly = TRUE)

library("GenomicSEM")
library("readr")
library("dplyr")
library("ggplot2")
library("qqman")
library("shru")
library("data.table")
library("utils")
library("reshape2")
library("devtools")
library("corrplot")
library("Matrix")


files_input = c( #########################################################################
  "/scratch/prj/gwas_sumstats/cleaned/SCHI06.gz", 
  "/users/k2473476/proj/NewGwasSub/DEPR14_ori_modified.gz", 
  "/scratch/prj/gwas_sumstats/cleaned/BIPO03.gz")
ref_file = "~/proj/CommonFiles/reference.1000G.maf.0.005.txt.gz" ####################################################
hm3 = "~/proj/CommonFiles/w_hm3.snplist"  ############################################################
paths_corr <- "/scratch/prj/gwas_sumstats/munged/" ############################################################
ld <- "~/proj/CommonFiles/eur_w_ld_chr/" ############################################################
wld <- ld


#traitnames =c("SCZ", "MDD", "BD")
#latentnames =c("Psychosis", "Depression", "Mania")
#output_name <- "Mania_Supermunged"

traitnames =c("MAN1", "MAN2", "MAINMAN")
latentnames =c("LAT1", "LAT2", "MAINLAT")
output_name <- "SYNTHETIC_PHENOTYPE"


infofilter = 0.9
maffilter = 0.05

sample.prev <- c(0.425, 0.346, 0.1013) ############################################################
population.prev <- c(0.01, 0.10, 0.02) ############################################################

se.logit_vector = c(FALSE, FALSE, FALSE) ############################################################
OLS_vector=c(TRUE, TRUE, TRUE) ############################################################
linprob_vector= c(FALSE, FALSE, FALSE) ############################################################
ncores = 32 ############################################################





get_script_directory <- function() {
  # Check if the script is being run interactively
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  
  # Check command-line arguments
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args)
  if (length(file_arg) > 0) {
    script_path <- sub("--file=", "", args[file_arg])
    return(dirname(normalizePath(script_path)))
  }
  
  # Check if running in an RStudio session
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
  }
  
  # Check common environment variables in cluster environments
  if (!is.null(Sys.getenv("PBS_O_WORKDIR", unset = NULL))) {
    return(normalizePath(Sys.getenv("PBS_O_WORKDIR")))
  }
  
  if (!is.null(Sys.getenv("SLURM_SUBMIT_DIR", unset = NULL))) {
    return(normalizePath(Sys.getenv("SLURM_SUBMIT_DIR")))
  }
  
  # If none of the above methods work, return NA
  return(NA)
}

# Get and print the script directory
script_directory <- get_script_directory()
cat(script_directory, "\n")

setwd(script_directory) ###############################################################################

new_dir_name <- args[1]
writeLines(new_dir_name, "dir.txt")
new_dir_path <- file.path(script_directory, new_dir_name) ##############################################################################
if (!dir.exists(new_dir_path)) {
  dir.create(new_dir_path)
}
setwd(new_dir_path)




 corr <- read.csv("../Correlation_input.csv")
 corr$sampleprev <- as.character(corr$sampleprev)
 corr$popprev <- as.character(corr$popprev)

 # Replacing "," with "." in the columns
 corr$sampleprev <- gsub(",", ".", corr$sampleprev)
 corr$popprev <- gsub(",", ".", corr$popprev)

 # Converting the columns back to numeric
 corr$sampleprev <- as.numeric(corr$sampleprev)
 corr$popprev <- as.numeric(corr$popprev)



 output_name_gz <-paste0(output_name, ".gz")

 corr_input <- paste0(paths_corr, corr$code) ############# path to munged sumstats
 corr_input <- c(corr_input)
 corr_input <- c(corr_input, output_name_gz)



print("---------------------------------------------------------PreprocessingFinished - SubGwas")

if("1" %in% args){
print("---------------------------------------------------------Part 1 Started - Preparation")

munge(files = files_input,
      hm3 = hm3,
      trait.names=traitnames,
      info.filter = infofilter,
      maf.filter = maffilter)

LDSCtraits <- paste0(traitnames, ".sumstats.gz")

ldsc_file <- ldsc(traits = LDSCtraits,sample.prev = sample.prev, population.prev = population.prev, ld = ld, wld = wld, trait.names = traitnames)
saveRDS(ldsc_file, file="LDSC_main.rds")

SNP_files <- sumstats(files = files_input,
                      ref = ref_file,
                      trait.names=traitnames, 
                      info.filter = infofilter, 
                      maf.filter = maffilter,
                      se.logit = se.logit_vector,
                      OLS=OLS_vector,
                      linprob=linprob_vector)

saveRDS(SNP_files, file="SNP_files.rds")
  
print("---------------------------------------------------------Part 1 Finished - Preparation")
}


if("2" %in% args){
  print("---------------------------------------------------------Part 2 Started - Model")

LDSCoutput <- readRDS(file="LDSC_main.rds")
SNP_files <- readRDS(file="SNP_files.rds")

model <-
         'Psychosis=~NA*BD + NA*SCZ
         Mania=~NA*BD * start(1)*BD
         Depression=~NA*BD + NA*MDD

         Mania~~1*Mania
         Psychosis~~1*Psychosis
         Depression~~1*Depression
         Mania~~0*Depression
         Psychosis~~NA*Depression
         Mania~~0*Psychosis

         SCZ ~~ 0*BD
         SCZ~~0*SCZ
         BD~~0*BD
         MDD~~0*SCZ
         MDD~~0*BD
         MDD~~0*MDD
         '
model <- gsub("Psychosis", latentnames[1], model)
model <- gsub("Depression", latentnames[2], model)
model <- gsub("Mania", latentnames[3], model)
model <- gsub("SCZ", traitnames[1], model)
model <- gsub("MDD", traitnames[2], model)
model <- gsub("BD", traitnames[3], model)

  
model_LD <- usermodel(LDSCoutput,estimation="DWLS",model=model, CFIcalc = TRUE, std.lv=TRUE)
print(model_LD)

output_file <- "model_LD.csv"
write.csv(model_LD, file = output_file, row.names = FALSE)
cat("Output exported to", output_file)

print("---------------------------------------------------------Part 2 Finished - Model")
}


if("3" %in% args){
  print("---------------------------------------------------------Part 3 Started - GWAS")

LDSCoutput <- readRDS(file="LDSC_main.rds")
SNP_files <- readRDS(file="SNP_files.rds")

if ("r" %in% args) {
  print("Running in Test Mode")
  SNP_files <- SNP_files %>% 
    filter(CHR == 2)
}

sub <- paste(latentnames, "~SNP", sep="")

model <-
         'Psychosis=~NA*BD + NA*SCZ
         Mania=~NA*BD * start(1)*BD
         Depression=~NA*BD + NA*MDD

         Mania~SNP
         Depression~SNP
         Psychosis~SNP

         Mania~~1*Mania
         Psychosis~~1*Psychosis
         Depression~~1*Depression
         Mania~~0*Depression
         Psychosis~~NA*Depression
         Mania~~0*Psychosis

         SCZ ~~ 0*BD
         SCZ~~0*SCZ
         BD~~0*BD
         MDD~~0*SCZ
         MDD~~0*BD
         MDD~~0*MDD

         SNP~~SNP
         '
model <- gsub("Psychosis", latentnames[1], model)
model <- gsub("Depression", latentnames[2], model)
model <- gsub("Mania", latentnames[3], model)
model <- gsub("SCZ", traitnames[1], model)
model <- gsub("MDD", traitnames[2], model)
model <- gsub("BD", traitnames[3], model)

output <- userGWAS(
  covstruc = LDSCoutput, 
  SNPs = SNP_files,
  model = model,
  printwarn = TRUE,
  parallel = T,
  smooth_check = T,
  cores = ncores,
  TWAS = F,
  MPI = F,
  sub = sub,
  estimation = "DWLS", 
  toler = 1e-40,
  GC="standard"
)
saveRDS(output, "output.rds")

data.table::fwrite(x = output[[3]],file = paste0(latentnames[3], "_gwas.gz"), append = F,quote = F,sep = "\t",col.names = T,nThread=20)

print("---------------------------------------------------------Part 3 Finished - GWAS")
}


if("4" %in% args){
  print("---------------------------------------------------------Part 4 Started - Plots and Post-Munge")
  
superm <- supermunge(filePaths = c(paste0(latentnames[3], "_gwas.gz")),
                     refFilePath = ref_file,
                     traitNames = c(output_name), 
                     process=TRUE)

maniadf <- fread(output_name_gz, fill = TRUE)

mahattan <- shru::plot.manhattan.custom(df = as.data.frame(maniadf)) + geom_hline(yintercept = 8, color = "red", linetype = "dashed")
ggsave("outputMANN.png", mahattan, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)
qq <- shru::plot.qq.custom(ps = maniadf$P)
ggsave("outputQQ.png", qq, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)

LDSCoutput<-ldsc(traits = c(paste0(traitnames, ".sumstats.gz"), output_name_gz), 
                 sample.prev = c(sample.prev, NA), 
                 population.prev = c(population.prev, NA), 
                 ld = ld, wld = wld, stand = T)
print(LDSCoutput)

print("---------------------------------------------------------Part 4 Finished - Plots and Post-Munge")
}

if("5" %in% args){
  print("---------------------------------------------------------Part 5 Started - Genetic Correlation")

 Correlation_output<-ldsc(traits = corr_input, 
                 sample.prev = corr$sampleprev, 
                 population.prev = corr$popprev, 
                 ld = ld, wld = wld, stand = T)

 saveRDS(Correlation_output, file="Correlation_output.rds")
 print("---------------------------------------------------------Part 5 Finished - Genetic Correlation")
}

if("6" %in% args){
  print("---------------------------------------------------------Part 6 Started - Matrix")

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

trait_names <- c(corr$trait, latentnames[3])
cluster <- c(corr$Cluster)

mega_correlation <- readRDS("Correlation_output.rds")

rownames(mega_correlation$S_Stand) <- trait_names
colnames(mega_correlation$S_Stand) <- trait_names

S <- mega_correlation$S_Stand

SE_R <- matrix(nrow = nrow(S), ncol = ncol(S))
rownames(SE_R) <- trait_names
colnames(SE_R) <- trait_names

SE_R[lower.tri(SE_R, diag = TRUE)] <- sqrt(diag(mega_correlation$V_Stand))
SE_R[upper.tri(SE_R, diag = TRUE)] <- t(SE_R)[upper.tri(SE_R, diag = TRUE)]

R_mania <- S[latentnames[3], ]
R_bd <- S[traitnames[3], ]
seR_mania <- SE_R[latentnames[3], ]
seR_bd <- SE_R[traitnames[3], ]

p_overall <- difftest.matrix(effectiveNumberOfTests = 35, mValues1 = as.matrix(R_mania), mStandard_errors1 = as.matrix(seR_mania), mValues2 = as.matrix(R_bd), mStandard_errors2 = as.matrix(seR_bd), mCorrelationEstimate.values = 0.897847609)

p_adjusted <- p_overall$pTest.values.adj[!is.na(p_overall$pTest.values.adj)]
paired_pval <- cbind(trait_names, p_adjusted)

R_mania <- data.frame(Trait = names(R_mania), Correlation = R_mania)
R_bd <- data.frame(Trait = names(R_bd), Correlation = R_bd)
seR_mania <- data.frame(Trait = names(seR_mania), Correlation = seR_mania)
seR_bd <- data.frame(Trait = names(seR_bd), Correlation = seR_bd)

dfs <- list(R_mania, R_bd, seR_mania, seR_bd)

merged_df <- Reduce(function(x, y) merge(x, y, by = "Trait"), dfs)
new_colnames <- c("Trait", paste0("Correlation_", latentnames[3]), 
                    paste0("Correlation_", traitnames[3]), 
                    paste0("SE_", latentnames[3]), 
                    paste0("SE_", traitnames[3]))
                    
names(merged_df) <- new_colnames

paired_pval[, "p_adjusted"] <- as.numeric(paired_pval[, "p_adjusted"])
merged_df$p_adjusted <- paired_pval[match(merged_df$Trait, paired_pval[, "trait_names"]), "p_adjusted"]

merged_df_with_pvalues <- merged_df


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
 print("---------------------------------------------------------Part 6 Finished - Matrix")
}



print("---------------------------------------------------------Whole Code Finished - SubGwas")
