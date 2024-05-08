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

install_github("johanzvrskovec/shru")

setwd(file.path("~", "proj", "NewGwasSub"))

new_dir_name <- args[1]
writeLines(new_dir_name, "dir.txt")
new_dir_path <- file.path("~", "proj", "NewGwasSub", new_dir_name)
if (!dir.exists(new_dir_path)) {
  dir.create(new_dir_path)
}
setwd(new_dir_path)

#"/users/k2473476/proj/NewGwasSub/DEPR14_ori_modified.gz", 
#"/scratch/prj/gwas_sumstats/cleaned/DEPR14.gz", 

files_input = c(
  "/scratch/prj/gwas_sumstats/cleaned/SCHI06.gz", 
  "/users/k2473476/proj/NewGwasSub/DEPR14_ori_modified.gz", 
  "/scratch/prj/gwas_sumstats/cleaned/BIPO03.gz")
ref_file = "~/proj/CommonFiles/reference.1000G.maf.0.005.txt.gz"
hm3 = "~/proj/CommonFiles/w_hm3.snplist" 

traitnames =c("SCZ", "MDD", "BD")
latentnames =c("Psychosis", "Depression", "Mania")



#samplesizes <- c(175799,173005,413466)
infofilter = 0.9
maffilter = 0.05
ld <- "~/proj/CommonFiles/eur_w_ld_chr/"
wld <- ld

sample.prev <- c(0.425, 0.346, 0.1013)
#BD https://www.sciencedirect.com/science/article/pii/S0165032720332419
# 0.0048 An evaluation of variation in published estimates of schizophrenia prevalence from 1990─2013: a systematic literature review
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075362%200. prev 25
population.prev <- c(0.01, 0.10, 0.02)
OLS_vector=c(TRUE, TRUE, TRUE)
linprob_vector= c(FALSE, FALSE, FALSE)
se.logit_vector = c(FALSE, FALSE, FALSE)
ncores = 32



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

LDSCnomod<-ldsc(traits = LDSCtraits,sample.prev = sample.prev, population.prev = population.prev, ld = ld, wld = wld, trait.names = LDSCinput)
saveRDS(LDSCnomod, file="LDSC_main.rds")

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

model_locked <-
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
model_locked_new <-
         'NonMania=~NA*BD + NA*SCZ + NA*MDD
         Mania=~NA*BD

         Mania~~1*Mania
         NonMania~~1*NonMania
         Mania~~NA*NonMania

         SCZ ~~ 0*BD
         SCZ~~0*SCZ
         BD~~0*BD
         MDD~~0*SCZ
         MDD~~0*BD
         MDD~~0*MDD
         '



output_without_SNPs<-usermodel(LDSCoutput,estimation="DWLS",model=model_locked, CFIcalc = TRUE, std.lv=TRUE)


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
#latentnames =c("NonMania", "Mania")
latentnames =c("Mania")
sub <- paste(latentnames, "~SNP", sep="")

#print(head(SNP_files))
print(latentnames)
#NA,0,NA,NA,NA,NA,"Psychosis","=~","BD",0.191573078732575,"0.00294251509803735",0.649587525517403,"NaN",0.649587525517403,"< 5e-300"
#NA,0,NA,NA,NA,NA,"Mania","=~","BD",0.367977656013408,"0.00602764740351205",0.708409201906156,"0.0251122599602286",0.708409201906156,"< 5e-300"
#NA,0,NA,NA,NA,NA,"Depression","=~","BD",0.146921984326957,"NaN",0.359313585728609,"0.0338523771015311",0.359313585728609,NA
#NA,0,NA,NA,NA,NA,"Mania","~~","Depression",-0.280447316292466,"NaN",-0.205446284239331,"NaN",-0.205446284239331,NA
#NA,0,NA,NA,NA,NA,"Psychosis","~~","Depression",0.0481095531668951,"0.004540724059861",0.398014814122653,"0.0375680083940393",0.398014814122653,"3.13910863467453e-26"
#NA,0,NA,NA,NA,NA,"Psychosis","~~","Mania",-0.285683296600241,"0.0146285111116016",-0.145730346273058,"0.10974752258409",-0.145730346273058,"6.19883878014931e-85"

#0.6893630 4.63113939756803e-110
#0.5492663              < 5e-300
#0.4723190  3.31685094248736e-32
model_locked <-
         'Psychosis=~NA*BD + NA*SCZ
         Mania=~NA*BD * start(1)*BD
         Depression=~NA*BD + NA*MDD

         Mania~SNP

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
model_locked_new <-
         'NonMania=~NA*BD + NA*SCZ + NA*MDD
         Mania=~NA*BD

         NonMania~SNP
         Mania~SNP

         Mania~~1*Mania
         NonMania~~1*NonMania
         Mania~~NA*NonMania

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
  model = model_locked,
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
  
#psydf <- fread("Psychosis_gwas.gz", fill = TRUE)
#depredf <- fread("Depression_gwas.gz", fill = TRUE)
superm <- supermunge(filePaths = c("Mania_gwas.gz"),
                     refFilePath = ref_file,
                     traitNames = c("Mania_Supermunged"), 
                     process=TRUE)


maniadf <- fread("Mania_Supermunged.gz", fill = TRUE)

#mahattan <- shru::plot.manhattan.custom(df = as.data.frame(psydf)) + geom_hline(yintercept = 8, color = "red", linetype = "dashed")
#ggsave("psyMANN.png", mahattan, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)
#qq <- shru::plot.qq.custom(ps = psydf$Pval_Estimate)
#ggsave("psyQQ.png", qq, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)

#mahattan <- shru::plot.manhattan.custom(df = as.data.frame(depredf)) + geom_hline(yintercept = 8, color = "red", linetype = "dashed")
#ggsave("depreMANN.png", mahattan, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)
#qq <- shru::plot.qq.custom(ps = depredf$Pval_Estimate)
#ggsave("depreQQ.png", qq, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)

mahattan <- shru::plot.manhattan.custom(df = as.data.frame(maniadf)) + geom_hline(yintercept = 8, color = "red", linetype = "dashed")
ggsave("maniaMANN.png", mahattan, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)
qq <- shru::plot.qq.custom(ps = maniadf$P)
ggsave("maniaQQ.png", qq, width = 297, height = 210, units = "mm", dpi = 100, scale = 0.8)


LDSCoutput<-ldsc(traits = "Mania_Supermunged.gz", 
                 sample.prev = NA, 
                 population.prev = NA, 
                 ld = ld, wld = wld, stand = T)
print(LDSCoutput)


LDSCoutput<-ldsc(traits = c(paste0(traitnames, ".sumstats.gz"), "Mania_Supermunged.gz"), 
                 sample.prev = c(sample.prev, NA), 
                 population.prev = c(population.prev, NA), 
                 ld = ld, wld = wld, stand = T)
print(LDSCoutput)


print("---------------------------------------------------------Part 4 Finished - Plots and Post-Munge")
}


if("5" %in% args){
  print("---------------------------------------------------------Part 5 Started - Genetic Correlation")
 corr <- read.csv("/users/k2473476/proj/NewGwasSub/corr.csv")
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

 saveRDS(LDSCoutput, file="MegaCorrelation.rds")


 print("---------------------------------------------------------Part 5 Finished - Genetic Correlation")
}

print("---------------------------------------------------------Whole Code Finished - SubGwas")