library(survival)
library(survminer)
library(gtsummary)

################################################################################
#
# Kaplan Meier plotter; hazard ratio calculator
#
################################################################################

# load scraped data
rsem <- readRDS("~/Documents/BFX_proj/_data_public/TCGA_PanCancer/_processed/tcga_pancancer_rsem_FZ01.rds")
metadata <- readRDS("~/Documents/BFX_proj/_data_public/TCGA_PanCancer/_processed/tcga_pancancer_metadata_FZ01.rds")

# analysis
metadata_ <- metadata[metadata$CANCER_TYPE_ACRONYM == "KIRC", ]
rsem_ <- rsem[, metadata_$SAMPLE_ID]
surv_df_ <- data.frame(time_to_event = metadata_[, "OS_MONTHS"],
                       event = as.numeric(substr((metadata_[, "OS_STATUS"]), 1, 1)), 
                       strata = metadata_[, "PATH_T_STAGE"])
ggsurvplot(survfit(Surv(time_to_event, event) ~ strata, data = surv_df_), pval = T)


# # p-value does not calculate correctly when inside a function; do not use
# # base_function
# base_km <- function(data = metadata_,
#                     time_to_event = "OS_MONTHS",
#                     event = "OS_STATUS",
#                     strata = "PATH_T_STAGE"){
#   out <- list()
#   out[["surv_data"]] <- surv_df_ <- data.frame(time_to_event = data[, time_to_event],
#                          event = as.numeric(substr((data[, event]), 1, 1)), 
#                          strata = data[, strata])
#   out[["hazard_ratio"]] <- tbl_regression(coxph(Surv(time_to_event, event) ~ strata, data = surv_df_), exp = T)
#   out[["kaplan_meier"]] <- ggsurvplot(survfit(Surv(time_to_event, event) ~ strata, data = surv_df_), pval = T)
#   return(out)
# }



            