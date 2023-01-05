library(cBioPortalData)
library(AnVIL)
library(stringr)
library(tidyr)

cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)

################################################################################
#
# TCGA PanCancer Atlas
#
################################################################################

###---scrape data---------------------------------------------------------------
study_data_ <- studies[grepl("TCGA PanCancer", studies$description), ]
studies_ <- data.frame(studyId = study_data_$studyId, loaded = NA, method = NA, row.names = study_data_$studyId)

for(i in studies_$studyId){
  if(i %in% rownames(studies_)[studies_$loaded != "fail"]){
    cat(paste0("\n", i, " already downloaded\n"))
    next
  }
  try(data_ <- cBioDataPack(i, ask = F))
  if(exists("data_")){
    saveRDS(data_, paste0("~/Documents/BFX_proj/_data_public/TCGA_PanCancer/", i, ".rds"))
    studies_[i, "loaded"] <- date()
    studies_[i, "method"] <- "cBioDataPack()"
    unlink("~/.cache/cBioPortalData/")
    rm(data_)
  } else {
    studies_[i, "loaded"] <- "fail"
    studies_[i, "method"] <- "manual" # download from cBioPortal directly
  }
  write.csv(studies_, "~/Documents/BFX_proj/_data_public/TCGA_PanCancer/checkpoint.csv")
  rm(i)
}

###---extract and merge RSEM, metadata------------------------------------------
checkpoint <- read.csv("Documents/BFX_proj/_data_public/TCGA_PanCancer/checkpoint.csv", row.names = 1)

# scan for universal genes and metadata
## package loaded data
for(i in checkpoint$studyId[checkpoint$loaded != "fail"]){
  data_ <- readRDS(paste0("Documents/BFX_proj/Data/TCGA_PanCancer/", i,".rds"))
  rsem_ <- assay(data_, "mrna_seq_v2_rsem")
  meta_ <- colData(data_)
  if(exists("common_genes")){
    common_genes <- common_genes[common_genes %in% rownames(rsem_)]
  } else {
    common_genes <<- rownames(rsem_)
    common_genes <- common_genes[common_genes != ""]
  }
  if(exists("common_meta")){
    common_meta <- common_meta[common_meta %in% colnames(meta_)]
  } else {
    common_meta <<- colnames(meta_)
  }
  rm(i, list = ls()[grepl("_$", ls())])
}
## manually loaded data
for(i in checkpoint$studyId[checkpoint$loaded == "fail"]){
  data_file_ <- paste0("Documents/BFX_proj/_data_public/TCGA_PanCancer/", i,"/") # point to file
  rsem__ <- read.table(paste0(data_file_, "data_mrna_seq_v2_rsem.txt"), sep = "", header = T, fill = T) # extract rsem
  rsem_ <- as.matrix(rsem__[, 3:ncol(rsem__)])
  rownames(rsem_) <- rsem__$Hugo_Symbol
  patient_data_ <- as.data.frame(read.table(paste0(data_file_, "data_clinical_patient.txt"), sep = "\t", header = T, fill = T)) # extract meta
  sample_data_ <- as.data.frame(read.table(paste0(data_file_, "data_clinical_sample.txt"), sep = "\t", header = T, fill = T)) # extract meta
  meta_ <- left_join(patient_data_, sample_data_, by = "PATIENT_ID")
  
  if(exists("common_genes")){
    common_genes <- common_genes[common_genes %in% rownames(rsem_)]
  } else {
    common_genes <<- rownames(rsem_)
    common_genes <- common_genes[common_genes != ""]
  }
  if(exists("common_meta")){
    common_meta <- common_meta[common_meta %in% colnames(meta_)]
  } else {
    common_meta <<- colnames(meta_)
  }
  rm(i, list = ls()[grepl("_$", ls())])
}

# merge data
## package loaded data
for(i in checkpoint$studyId[checkpoint$loaded != "fail"]){
  data_ <- readRDS(paste0("Documents/BFX_proj/_data_public/TCGA_PanCancer/", i,".rds")) # load data
  rsem_ <- assay(data_, "mrna_seq_v2_rsem") # extract rsem
  rsem_ <- rsem_[common_genes, ]
  meta_ <- as.data.frame(colData(data_), row.names = colData(data_)$SAMPLE_ID) # extract meta
  meta_ <- meta_[, common_meta] # filter and order meta to samples with RNA data; select common metadata
  # call samples with counts and metadata
  samples_ <- rownames(meta_)[rownames(meta_) %in% colnames(rsem_)] 
  rsem_ <- rsem_[, samples_]
  meta_ <- meta_[samples_, ]
  
  if(exists("rsem")){
    rsem <- cbind(rsem, rsem_)
  } else {
    rsem <<- rsem_
  }
  
  if(exists("meta")){
    meta <- rbind(meta, meta_)
  } else {
    meta <<- meta_
  }
  rm(i, list = ls()[grepl("_$", ls())])
}
## manually loaded data
for(i in checkpoint$studyId[checkpoint$loaded == "fail"]){
  data_file_ <- paste0("Documents/BFX_proj/_data_public/TCGA_PanCancer/", i,"/") # point to file
  rsem__ <- read.table(paste0(data_file_, "data_mrna_seq_v2_rsem.txt"), sep = "", header = T, fill = T) # extract rsem
  rsem_ <- as.matrix(rsem__[, 3:ncol(rsem__)])
  rownames(rsem_) <- rsem__$Hugo_Symbol
  colnames(rsem_) <- str_replace_all(colnames(rsem__)[3:ncol(rsem__)], "\\.", "-")
  rsem_ <- rsem_[common_genes, ]
  patient_data_ <- as.data.frame(read.table(paste0(data_file_, "data_clinical_patient.txt"), sep = "\t", header = T, fill = T)) # extract meta
  sample_data_ <- as.data.frame(read.table(paste0(data_file_, "data_clinical_sample.txt"), sep = "\t", header = T, fill = T)) # extract meta
  meta_ <- left_join(sample_data_, patient_data_, by = "PATIENT_ID")
  meta_ <- meta_[, common_meta]
  rownames(meta_) <- meta_$SAMPLE_ID
  # call samples with counts and metadata
  samples_ <- rownames(meta_)[rownames(meta_) %in% colnames(rsem_)] 
  rsem_ <- rsem_[, samples_]
  meta_ <- meta_[samples_, ]
  
  if(exists("rsem")){
    rsem <- cbind(rsem, rsem_)
  } else {
    rsem <<- rsem_
  }
  
  if(exists("meta")){
    meta <- rbind(meta, meta_)
  } else {
    meta <<- meta_
  }
  rm(i, list = ls()[grepl("_$", ls())])
}

saveRDS(rsem, "Documents/BFX_proj/_data_public/TCGA_PanCancer/_processed/tcga_pancancer_rsem_FZ01.rds")
saveRDS(meta, "Documents/BFX_proj/_data_public/TCGA_PanCancer/_processed/tcga_pancancer_metadata_FZ01.rds")
