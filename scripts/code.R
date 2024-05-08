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
library(devtools)
library(corrplot)

setwd(file.path("~", "proj", "NewGwasSub")) ###############################################################################

new_dir_name <- args[1]
writeLines(new_dir_name, "dir.txt")
new_dir_path <- file.path("~", "proj", "NewGwasSub", new_dir_name) ##############################################################################
if (!dir.exists(new_dir_path)) {
  dir.create(new_dir_path)
}
setwd(new_dir_path)

files_input = c( #########################################################################
  "/scratch/prj/gwas_sumstats/cleaned/SCHI06.gz", 
  "/users/k2473476/proj/NewGwasSub/DEPR14_ori_modified.gz", 
  "/scratch/prj/gwas_sumstats/cleaned/BIPO03.gz")
ref_file = "~/proj/CommonFiles/reference.1000G.maf.0.005.txt.gz" ####################################################
hm3 = "~/proj/CommonFiles/w_hm3.snplist"  ############################################################

traitnames =c("SCZ", "MDD", "BD") 
latentnames =c("Psychosis", "Depression", "Mania")

infofilter = 0.9 ############################################################
maffilter = 0.05 ############################################################
ld <- "~/proj/CommonFiles/eur_w_ld_chr/" ############################################################
wld <- ld

sample.prev <- c(0.425, 0.346, 0.1013) ############################################################
population.prev <- c(0.01, 0.10, 0.02) ############################################################
OLS_vector=c(TRUE, TRUE, TRUE) ############################################################
linprob_vector= c(FALSE, FALSE, FALSE) ############################################################
se.logit_vector = c(FALSE, FALSE, FALSE) ############################################################
ncores = 32 ############################################################


print("---------------------------------------------------------PreprocessingFinished - SubGwas")

if("1" %in% args){
print("---------------------------------------------------------Part 1 Started - Preparation")

munge(files = files_input,
      hm3 = hm3,
      trait.names=traitnames, 
      #N = samplesizes,
      info.filter = infofilter, 
      maf.filter = maffilter)

LDSCtraits <- paste0(traitnames, ".sumstats.gz")
LDSCinput<-traitnames

ldsc_file<-ldsc(traits = LDSCtraits,sample.prev = sample.prev, population.prev = population.prev, ld = ld, wld = wld, trait.names = LDSCinput)
saveRDS(ldsc_file, file="LDSC_main.rds")

SNP_files <- sumstats(files = files_input,
                      ref = ref_file,
                      trait.names=traitnames, 
                      #N = samplesizes,
                      info.filter = infofilter, 
                      maf.filter = maffilter,
                      OLS=OLS_vector,
                      linprob=linprob_vector,
                      se.logit = se.logit_vector)

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
output_without_SNPs<-usermodel(LDSCoutput,estimation="DWLS",model=model, CFIcalc = TRUE, std.lv=TRUE)

print(output_without_SNPs)

output_file_without_SNPs <- "model_without_SNPs.csv"
write.csv(output_without_SNPs, file = output_file_without_SNPs, row.names = FALSE)
cat("Output exported to", output_file_without_SNPs)

print("---------------------------------------------------------Part 2 Finished - Model")
}


if("3" %in% args){
  print("---------------------------------------------------------Part 3 Started - GWAS")

LDSCoutput <- readRDS(file="LDSC_main.rds")
SNP_files <- readRDS(file="SNP_files.rds")

if ("r" %in% args) {
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

data.table::fwrite(x = output[[1]],file = paste0(latentnames[1], "_gwas.gz"), append = F,quote = F,sep = "\t",col.names = T,nThread=20)
#data.table::fwrite(x = output[[2]],file = paste0(latentnames[2], "_gwas.gz"), append = F,quote = F,sep = "\t",col.names = T,nThread=20)
#data.table::fwrite(x = output[[3]],file = paste0(latentnames[3], "_gwas.gz"), append = F,quote = F,sep = "\t",col.names = T,nThread=20)
print("---------------------------------------------------------Part 3 Finished - GWAS")
}


if("4" %in% args){
  print("---------------------------------------------------------Part 4 Started - Plots and Post-Munge")
  
superm <- supermunge(filePaths = c("Mania_gwas.gz"),
                     refFilePath = ref_file,
                     traitNames = c("Mania_Supermunged"), 
                     process=TRUE)

maniadf <- fread("Mania_Supermunged.gz", fill = TRUE)

mahattan <- shru::plot.manhattan.custom(df = as.data.frame(maniadf)) + geom_hline(yintercept = 8, color = "red", linetype = "dashed")
ggsave("maniaMANN.png", mahattan, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)
qq <- shru::plot.qq.custom(ps = maniadf$P)
ggsave("maniaQQ.png", qq, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)

LDSCoutput<-ldsc(traits = c(paste0(traitnames, ".sumstats.gz"), "Mania_Supermunged.gz"), 
                 sample.prev = c(sample.prev, NA), 
                 population.prev = c(population.prev, NA), 
                 ld = ld, wld = wld, stand = T)
print(LDSCoutput)

print("---------------------------------------------------------Part 4 Finished - Plots and Post-Munge")
}

if("5" %in% args){
  print("---------------------------------------------------------Part 5 Started - Genetic Correlation")
 corr <- read.csv("/users/k2473476/proj/NewGwasSub/corr.csv") ###############################################################################
 corr$sampleprev <- as.character(corr$sampleprev)
 corr$popprev <- as.character(corr$popprev)

 # Replacing "," with "." in the columns
 corr$sampleprev <- gsub(",", ".", corr$sampleprev)
 corr$popprev <- gsub(",", ".", corr$popprev)

 # Converting the columns back to numeric
 corr$sampleprev <- as.numeric(corr$sampleprev)
 corr$popprev <- as.numeric(corr$popprev)

 files_input <- paste0("/scratch/prj/gwas_sumstats/munged/", corr$code, ".sumstats.gz")
 files_input <- c(files_input)

 files_input <- c(files_input, "Mania_Supermunged.gz")

 print(length(files_input))

 traitnames =corr$code
 samplesizes <- corr$sample_size_discovery
 sample.prev <- corr$sampleprev
 population.prev <- corr$popprev
 
 LDSCoutput<-ldsc(traits = files_input, 
                 sample.prev = corr$sampleprev, 
                 population.prev = corr$popprev, 
                 ld = ld, wld = wld, stand = T)

 saveRDS(LDSCoutput, file="Correlation_output.rds")
 print("---------------------------------------------------------Part 5 Finished - Genetic Correlation")
}

print("---------------------------------------------------------Whole Code Finished - SubGwas")
