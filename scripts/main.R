# Tested with R 4.3.0 

files_input = c( #########################################################################
  "/scratch/prj/gwas_sumstats/cleaned/SCHI06.gz", 
  "/users/k2473476/proj/Recover/DEPR14_ori_modified.gz", 
  "/scratch/prj/gwas_sumstats/cleaned/BIPO03.gz")
ref_file = "~/proj/CommonFiles/reference.1000G.maf.0.005.txt.gz" ####################################################
hm3 = "~/proj/CommonFiles/w_hm3.snplist"  ############################################################
paths_corr <- "/scratch/prj/gwas_sumstats/munged/" ############################################################
ld <- "~/proj/CommonFiles/eur_w_ld_chr/" ############################################################
wld <- ld


#traitnames =c("SCZ", "MDD", "BD")
#latentnames =c("Psychosis", "Depression", "Mania")
#output_name <- "Mania_Supermunged"

traitnames =c("MAN1", "MAN2", "BD") ### 3rd has to appear in the input csv
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

setwd(script_directory)

new_dir_name <- args[1]
writeLines(new_dir_name, "dir.txt")
new_dir_path <- file.path(script_directory, new_dir_name)
if (!dir.exists(new_dir_path)) {
  dir.create(new_dir_path)
}
setwd(new_dir_path)
print(new_dir_path)



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

 corr_input <- paste0(paths_corr, corr$code)
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

  psychosis_on_bd <- round(model_LD$results[1, "STD_Genotype"], 2)
  mania_on_bd <- round(model_LD$results[3, "STD_Genotype"], 2)
  depression_on_bd <- round(model_LD$results[4, "STD_Genotype"], 2)
  psychosis_on_depression <- round(model_LD$results[6, "STD_Genotype"], 2)
  
  square_size <- 1       # Size of the squares
  circle_radius <- 0.5   # Radius of the circles
  x_distance <- 3        # Horizontal distance between the centers of the shapes
  y_distance <- 1        # Vertical distance between the rows
  
  square_labels <- c(traitnames[2], traitnames[3], traitnames[1])
  circle_labels <- c(latentnames[2], latentnames[3], latentnames[1])
  arrow_labels <- c("1.00", mania_on_bd, "1.00", depression_on_bd, psychosis_on_bd, psychosis_on_depression)

  squares <- data.frame(
    x = rep(c(1, 1 + x_distance, 1 + 2 * x_distance), each = 4),
    y = rep(0, 12),
    group = rep(1:3, each = 4)
  )
  
  squares <- squares %>%
    group_by(group) %>%
    mutate(
      x = x + c(-square_size / 2, square_size / 2, square_size / 2, -square_size / 2),
      y = y + c(-square_size / 2, -square_size / 2, square_size / 2, square_size / 2)
    )
  
  circles <- data.frame(
    x = rep(c(1, 1 + x_distance, 1 + 2 * x_distance), each = 100),
    y = rep(square_size + y_distance, 300),
    angle = rep(seq(0, 2 * pi, length.out = 100), 3),
    group = rep(1:3, each = 100)
  )
  
  center_square <- data.frame(
    x = c(1 + x_distance, 1 + x_distance, 1 + x_distance, 1 + x_distance),
    y = c(2 * square_size + 2 * y_distance, 2 * square_size + 2 * y_distance, 2 * square_size + 2 * y_distance, 2 * square_size + 2 * y_distance),
    group = 1
  )
  
  center_square <- center_square %>%
    group_by(group) %>%
    mutate(
      x = x + c(-square_size / 2, square_size / 2, square_size / 2, -square_size / 2),
      y = y + c(-square_size / 2, -square_size / 2, square_size / 2, square_size / 2)
    )
  
  square_labels_df <- data.frame(
    x = c(1, 1 + x_distance, 1 + 2 * x_distance),
    y = rep(-square_size/2 - 0.2, 3),
    label = square_labels
  )
  
  circle_labels_df <- data.frame(
    x = c(1, 1 + x_distance, 1 + 2 * x_distance),
    y = rep(circle_radius + 2*y_distance + 0.2, 3),
    label = circle_labels
  )

  ggplot() +
    geom_polygon(data = squares, aes(x = x, y = y, group = group), fill = NA, color = "black") +
    geom_polygon(data = center_square, aes(x = x, y = y, group = group), fill = NA, color = "black") +
    geom_path(data = circles, aes(x = x + circle_radius * cos(angle), y = y + circle_radius * sin(angle), group = group), color = "black") +
    geom_segment(aes(x = 1 + x_distance, y = 2 * y_distance + 1.5 * square_size, xend = 1, yend = circle_radius + 2*y_distance), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 1 + x_distance, y = 2 * y_distance + 1.5 * square_size, xend = 1 + x_distance, yend = circle_radius + 2*y_distance), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 1 + x_distance, y = 2 * y_distance + 1.5 * square_size, xend = 1 + 2 * x_distance, yend = circle_radius + 2*y_distance), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 1 + 2 * x_distance, y = y_distance+circle_radius, xend = 4.5, yend = square_size/2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 1 + x_distance, y = y_distance+circle_radius, xend = 1 + x_distance, yend = square_size/2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 1, y = y_distance+circle_radius, xend = 4.5-square_size, yend = square_size/2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 1, y = y_distance+circle_radius, xend = 1, yend = square_size/2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 1 + 2 * x_distance, y = y_distance+circle_radius, xend = 1 + 2 * x_distance, yend = square_size/2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = 1, y = circle_radius + 2*y_distance, xend = 1 + 2 * x_distance, yend = circle_radius + 2*y_distance), curvature = -0.15, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = 1-square_size/4, y = -0.5, xend = 1+square_size/4, yend = -0.5), curvature = 0.8, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = 4-square_size/4, y = -0.5, xend = 4+square_size/4, yend = -0.5), curvature = 0.8, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = 7-square_size/4, y = -0.5, xend = 7+square_size/4, yend = -0.5), curvature = 0.8, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = 0.57, y = 2-0.25, xend = 0.57, yend = 2+0.25), curvature = -0.8, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = x_distance+0.57, y = 2-0.25, xend = x_distance+0.57, yend = 2+0.25), curvature = -0.8, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = 2*x_distance+0.43+1, y = 2-0.25, xend = 2*x_distance+0.43+1, yend = 2+0.25), curvature = 0.8, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_curve(aes(x = 4-square_size/4, y = 4.5, xend = 4+square_size/4, yend = 4.5), curvature = -0.8, arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_text(data = square_labels_df, aes(x = x, y = y+0.75, label = label), color = "black") +
    geom_text(data = circle_labels_df, aes(x = x, y = y-0.75, label = label), color = "black") +
    geom_text(aes(x = 4, y = 4, label = "Genetic Effect"), color = "black", size  = 3) +
    geom_text(aes(x = 0.7, y = 0.7, label = arrow_labels[1]), color = "black", vjust = -0.5) +
    geom_text(aes(x = 0.7 + x_distance, y = 1, label = arrow_labels[2]), color = "black", vjust = -0.5) +
    geom_text(aes(x = 0.3 + 2*x_distance + square_size, y = 0.7, label = arrow_labels[3]), color = "black", vjust = -0.5) +
    geom_text(aes(x = 2.4, y = 1, label = arrow_labels[4]), color = "black", vjust = -0.5) +
    geom_text(aes(x = 2.6+x_distance, y = 1, label = arrow_labels[5]), color = "black", vjust = -0.5) +
    geom_text(aes(x = 3.35, y = 2.9, label = arrow_labels[6]), color = "black", vjust = -0.5) +
    geom_text(aes(x = 1, y = -1, label = "0.00"), color = "black", vjust = -0.5) +
    geom_text(aes(x = 1+x_distance, y = -1, label = "0.00"), color = "black", vjust = -0.5) +
    geom_text(aes(x = 1+2*x_distance, y = -1, label = "0.00"), color = "black", vjust = -0.5) +
    geom_text(aes(x = 0.15, y = 1.9, label = "1.00"), color = "black", vjust = -0.5) +
    geom_text(aes(x = 0.15+x_distance, y = 1.9, label = "1.00"), color = "black", vjust = -0.5) +
    geom_text(aes(x = 1.85+2*x_distance, y = 1.9, label = "1.00"), color = "black", vjust = -0.5) +
    geom_text(aes(x = 1+x_distance, y = 4.75, label = "2pq"), color = "black", vjust = -0.5) +
    theme_void() +
    coord_fixed()
     
  ggsave("SEMplot.png", width = 6, height = 6, units = "in")

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

# Function to add significance stars based on p-value
add_significance <- function(p_value) {
  if (p_value < 0.001) "***"
  else if (p_value < 0.01) "**"
  else if (p_value < 0.05) "*"
  else ""
}

# Function to perform differential test between two matrices
difftest.matrix <- function(mValues1, mStandard_errors1, mValues2, mStandard_errors2, mCorrelationEstimate.values = 0.95, effectiveNumberOfTests = NULL, mValueCovariances = NULL, symmetric = TRUE, eigenSumLimit = 0.995) {
  
  mValues1 <- symMatrixConform(mValues1)
  mStandard_errors1 <- symMatrixConform(mStandard_errors1)
  mValues2 <- symMatrixConform(mValues2)
  mStandard_errors2 <- symMatrixConform(mStandard_errors2)
  
  if (is.null(effectiveNumberOfTests) & !is.null(mValueCovariances)) {
    effectiveNumberOfTests <- getEffectiveNumberOfTests(covarianceMatrix = mValueCovariances, symmetric = symmetric, eigenSumLimit = eigenSumLimit)
  }
  
  mTest.values <- abs(mValues1 - mValues2)
  sde <- sqrt((mStandard_errors1^2) + (mStandard_errors2^2) - 2 * mCorrelationEstimate.values * mStandard_errors1 * mStandard_errors2)
  pTest.values <- 2 * pnorm(q = mTest.values, sd = sde, lower.tail = FALSE)
  pTest.values.adj <- matrix(data = p.adjust(pTest.values, method = "fdr", n = effectiveNumberOfTests), nrow = nrow(pTest.values), ncol = ncol(pTest.values))
  rownames(pTest.values.adj) <- rownames(pTest.values)
  colnames(pTest.values.adj) <- colnames(pTest.values)
  
  list(pTest.values = pTest.values, pTest.values.adj = pTest.values.adj)
}

mega_correlation <- readRDS("Correlation_output.rds")

# Set up the trait names
trait_names <- c(corr$trait, latentnames[3])
original_order <- c(latentnames[3], corr$trait)

# Prepare the S matrix
S <- mega_correlation$S_Stand
rownames(S) <- colnames(S) <- trait_names

# Prepare the SE_R matrix
SE_R <- matrix(nrow = nrow(S), ncol = ncol(S))
rownames(SE_R) <- colnames(SE_R) <- trait_names
SE_R[lower.tri(SE_R, diag = TRUE)] <- sqrt(diag(mega_correlation$V_Stand))
SE_R[upper.tri(SE_R, diag = TRUE)] <- t(SE_R)[upper.tri(SE_R, diag = TRUE)]

# Extract relevant correlations and standard errors
R_mania <- S[latentnames[3], ]
R_bd <- S[traitnames[3], ]
seR_mania <- SE_R[latentnames[3], ]
seR_bd <- SE_R[traitnames[3], ]

# Perform differential test
p_overall <- difftest.matrix(effectiveNumberOfTests = length(trait_names), mValues1 = as.matrix(R_mania), mStandard_errors1 = as.matrix(seR_mania), mValues2 = as.matrix(R_bd), mStandard_errors2 = as.matrix(seR_bd))

# Adjust p-values
p_adjusted <- p_overall$pTest.values.adj[!is.na(p_overall$pTest.values.adj)]
paired_pval <- data.frame(trait_names, p_adjusted = as.numeric(p_adjusted))

# Combine data into a single dataframe
dfs <- list(
  data.frame(Trait = names(R_mania), Correlation = R_mania),
  data.frame(Trait = names(R_bd), Correlation = R_bd),
  data.frame(Trait = names(seR_mania), SE = seR_mania),
  data.frame(Trait = names(seR_bd), SE = seR_bd)
)

merged_df <- Reduce(function(x, y) merge(x, y, by = "Trait"), dfs)
merged_df$p_adjusted <- paired_pval[match(merged_df$Trait, paired_pval$trait_names), "p_adjusted"]
merged_df$significance <- sapply(merged_df$p_adjusted, add_significance)
merged_df$label_sign <- paste(merged_df$Trait, merged_df$significance)

# Rename columns dynamically
colnames(merged_df)[2] <- paste0("Correlation_", latentnames[3])
colnames(merged_df)[3] <- paste0("Correlation_", traitnames[3])
colnames(merged_df)[4] <- paste0("SE_", latentnames[3])
colnames(merged_df)[5] <- paste0("SE_", traitnames[3])

# Ensure the order of traits matches original_order
merged_df <- merged_df[match(original_order, merged_df$Trait), ]

# Final modifications and save to CSV
# Remove Trait.1 column if it exists
merged_df <- merged_df[, !grepl("Trait.1", colnames(merged_df))]

# Cluster information
cluster <- c(corr$Cluster)
merged_df$Cluster <- cluster[match(merged_df$Trait, trait_names)]

print(merged_df)
write.csv(merged_df, "trait_correlation_data.csv", row.names = FALSE)

lines <- readLines("../dir.txt")

# Modify the second and third lines
lines[2] <- traitnames[3]
lines[3] <- latentnames[3]

# Write the updated content back to the file
writeLines(lines, "../dir.txt")

 print("---------------------------------------------------------Part 6 Finished - Matrix")
}



print("---------------------------------------------------------Whole Code Finished - SubGwas")
