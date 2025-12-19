###########################################################################################
# This code runs the delays distribution with 5 parameters such as: DATE ONSET SYMPTOMS,  # 
# DATE ONSET SYMPTOMS, DATE REPORTING, DATE HOSPITALISATION, DATE RECOVERY.               # 
#                                                                                         #        
#                                                                                         #
# Cite as:                                                                                #
# Quijada M et al, Dengue Epidemiology and Transmission Intensity across Panama           #
# during 2000-2024                                                                        #
# Pre-print available at: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5461361%20  #
#                                                                                         #
# Description:                                                                            #
# See README.md                                                                           #
###########################################################################################

## Load Required Libraries
library(writexl)
library(readxl)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggpubr)
library("fitdistrplus")  # Fitting univariate distributions
library(flexsurv)  # Flexible survival models
library("gridGraphics")
library(boot)
library(lubridate)



# Add data base
# Replace your path here; using "path/to/..." as a placeholder.
# Read the dataset from an RDS file.
DR <- readRDS("path/to/LineList_2000_24.rds") #Delay report

# Convert Date Strings to Date Objects
# Ensure all date columns are formatted properly for calculations
DR$`DATE SAMPLE COLLECTION` <- as.Date(DR$`DATE SAMPLE COLLECTION`, format = "%Y-%m-%d")
DR$`DATE ONSET SYMPTOMS` <- as.Date(DR$`DATE ONSET SYMPTOMS`, format = "%Y-%m-%d")
DR$`DATE REPORTING` <- as.Date(DR$`DATE REPORTING`, format = "%Y-%m-%d")
DR$`DATE HOSPITALISATION` <- as.Date(DR$`DATE HOSPITALISATION`, format = "%Y-%m-%d")
DR$`DATE RECOVERY` <- as.Date(DR$`DATE RECOVERY`, format = "%Y-%m-%d")

# Calculate Delays in Days Between Key Dates
# Subtract dates to compute delays for various stages of medical reporting
DR$OSSC <- as.numeric(DR$`DATE SAMPLE COLLECTION` - DR$`DATE ONSET SYMPTOMS`)  # Onset to Sample Collection
DR$OSR <- as.numeric(DR$`DATE REPORTING` - DR$`DATE ONSET SYMPTOMS`)          # Onset to Reporting
DR$OSDH <- as.numeric(DR$`DATE HOSPITALISATION` - DR$`DATE ONSET SYMPTOMS`)   # Onset to Hospitalization
DR$OSDR <- as.numeric(DR$`DATE RECOVERY` - DR$`DATE ONSET SYMPTOMS`)          # Onset to Recovery
DR$HDR <- as.numeric(DR$`DATE RECOVERY` - DR$`DATE HOSPITALISATION`)          # Hospitalization to Recovery

# Select Relevant Columns
# Filter and organize the dataset for delay analysis
delay <- DR %>%
  dplyr::select(OSR, OSDH, OSSC, OSDR, HDR, Region, Sex, Age)

# Fit Univariate Distributions for Delays
# Subset the data to delays within a reasonable range (0 to 20 days)
delay1 <- subset(delay, OSSC >= 0 & OSSC <= 20)
delay2 <- subset(delay, OSR >= 0 & OSR <= 20)
delay3 <- subset(delay, OSDH >= 0 & OSDH <= 20)
delay4 <- subset(delay, OSDR >= 0 & OSDR <= 20)
delay5 <- subset(delay, HDR >= 0 & HDR <= 20)

# ANALYSIS BY DELAY CATEGORY (National Level Analyses)
# 1. Date Onset Symptoms - Date Sample Collection
data1 <- delay1$OSSC[delay1$OSSC > 0]  # Exclude non-positive delays
delay1$OSSC <- as.numeric(delay1$OSSC)  # Ensure the data is numeric

# Fit Candidate Distributions
descdist(delay1$OSSC, boot = 1000)  # Explore the data distribution
fw1 <- fitdist(data1, "weibull")  # Weibull distribution
fg1 <- fitdist(data1, "gamma")    # Gamma distribution
fln1 <- fitdist(data1, "lnorm")   # Log-normal distribution
fita1 <- fitdist(data1, "llogis") # Log-logistic distribution
cau1 <- fitdist(data1, "cauchy")  # Cauchy distribution

# 2. Date Onset Symptoms - Date Reporting
data2 <- delay2$OSR[delay2$OSR > 0]
descdist(delay2$OSR, boot = 1000)  # Explore the data distribution
fw2 <- fitdist(data2, "weibull")
fg2 <- fitdist(data2, "gamma")
fln2 <- fitdist(data2, "lnorm")
fita2 <- fitdist(data2, "llogis")
cau2 <- fitdist(data2, "cauchy")

# 3. Date Onset Symptoms - Date Hospitalization
data3 <- delay3$OSDH[delay3$OSDH > 0]
descdist(delay3$OSDH, boot = 1000)
fw3 <- fitdist(data3, "weibull")
fg3 <- fitdist(data3, "gamma")
fln3 <- fitdist(data3, "lnorm")
fita3 <- fitdist(data3, "llogis")
cau3 <- fitdist(data3, "cauchy")

# 4. Date Onset Symptoms - Date Recovery
data4 <- delay4$OSDR[delay4$OSDR > 0]
descdist(delay4$OSDR)
fw4 <- fitdist(data4, "weibull")
fg4 <- fitdist(data4, "gamma")
fln4 <- fitdist(data4, "lnorm")
fita4 <- fitdist(data4, "llogis")
cau4 <- fitdist(data4, "cauchy")

# 5. Hospitalization - Date Recovery
data5 <- delay5$HDR[delay5$HDR > 0]
descdist(delay5$HDR)
fw5 <- fitdist(data5, "weibull")
fg5 <- fitdist(data5, "gamma")
fln5 <- fitdist(data5, "lnorm")
fita5 <- fitdist(data5, "llogis")
cau5 <- fitdist(data5, "cauchy")

# Goodness-of-Fit Statistics
# Assess and compare distributions using AIC, BIC, and other criteria
results1 <- gofstat(list(fw1, fg1, fln1, fita1, cau1))
results2 <- gofstat(list(fw2, fg2, fln2, fita2, cau2))
results3 <- gofstat(list(fw3, fg3, fln3, fita3, cau3))
results4 <- gofstat(list(fw4, fg4, fln4, fita4, cau4))
results5 <- gofstat(list(fw5, fg5, fln5, fita5, cau5))

# Compile and Save Results
# Combine results into data frames for each category
combined_C <- rbind(
  data.frame(Category = "DATE ONSET SYMPTOMS - DATE SAMPLE COLLECTION", AIC = results1$aic, BIC = results1$bic),
  data.frame(Category = "DATE ONSET SYMPTOMS - DATE REPORTING", AIC = results2$aic, BIC = results2$bic),
  data.frame(Category = "DATE ONSET SYMPTOMS - DATE HOSPITALISATION", AIC = results3$aic, BIC = results3$bic),
  data.frame(Category = "DATE ONSET SYMPTOMS - DATE RECOVERY", AIC = results4$aic, BIC = results4$bic),
  data.frame(Category = "Hospitalisation - DATE RECOVERY", AIC = results5$aic, BIC = results5$bic)
)

fileC <- "~/path/to/CountryAIC.xlsx"
write_xlsx(combined_C, fileC)  # Save the results
