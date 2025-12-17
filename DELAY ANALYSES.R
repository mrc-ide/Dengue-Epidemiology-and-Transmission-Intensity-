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

#Load packages
library(writexl)
library(readxl)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggpubr)
library("fitdistrplus")  #allows to easily fit distributions.
library(flexsurv)
library("gridGraphics")
library(boot)
library(lubridate)



# Add data base
# Replace your path here; using "path/to/..." as a placeholder.
# Read the dataset from an RDS file.
DR <- readRDS("path/to/LineList_2000_24.rds") #Delay report

# A function to calculate time differences (delays) between two columns.
# Adds the calculated difference as a new column in the dataset.
calculate_delay <- function(data, col1, col2, new_col_name) {
  data[[new_col_name]] <- as.numeric(data[[col1]] - data[[col2]])
  return(data)
}

# Convert date columns to Date objects to enable date-based calculations.
# Ensure column names in the dataset match exactly; otherwise, update column names.
DR$`DATE SAMPLE COLLECTION` <- as.Date(DR$`DATE SAMPLE COLLECTION`, format = "%Y-%m-%d")
DR$`DATE ONSET SYMPTOMS` <- as.Date(DR$`DATE ONSET SYMPTOMS`, format = "%Y-%m-%d")
DR$`DATE REPORTING` <- as.Date(DR$`DATE REPORTING`, format = "%Y-%m-%d")
DR$`DATE HOSPITALISATION` <- as.Date(DR$`DATE HOSPITALISATION`, format = "%Y-%m-%d")
DR$`DATE RECOVERY` <- as.Date(DR$`DATE RECOVERY`, format = "%Y-%m-%d")


# Compute Delay Periods
# Calculate numerical delays in days 
# The result will be new variables appended to the data frame.
DR$OSSC <-  as.numeric(DR$`DATE SAMPLE COLLECTION` - DR$`DATE ONSET SYMPTOMS`) #Delay 1
DR$OSR <-  as.numeric(DR$`DATE REPORTING` - DR$`DATE ONSET SYMPTOMS`) #Delay 2
DR$OSDH <-  as.numeric(DR$`DATE HOSPITALISATION` - DR$`DATE ONSET SYMPTOMS`) #Delay 3
DR$OSDR <-  as.numeric(DR$`DATE RECOVERY` - DR$`DATE ONSET SYMPTOMS`) #Delay 4
DR$HDR <-  as.numeric(DR$`DATE RECOVERY` - DR$`DATE HOSPITALISATION`) #Delay 5


# A function to filter delays within a specific range.
filter_delay <- function(data, delay_col, range_min, range_max) {
  return(subset(data, data[[delay_col]] >= range_min & data[[delay_col]] <= range_max))
}

# Clean and validate date columns
date_columns <- c("DATE SAMPLE COLLECTION", "DATE ONSET SYMPTOMS", "DATE REPORTING", 
                  "DATE HOSPITALISATION", "DATE RECOVERY")

# Create new columns for each delay type
DR <- calculate_delay(DR, "DATE SAMPLE COLLECTION", "DATE ONSET SYMPTOMS", "OSSC") # Delay 1
DR <- calculate_delay(DR, "DATE REPORTING", "DATE ONSET SYMPTOMS", "OSR")          # Delay 2
DR <- calculate_delay(DR, "DATE HOSPITALISATION", "DATE ONSET SYMPTOMS", "OSDH")   # Delay 3
DR <- calculate_delay(DR, "DATE RECOVERY", "DATE ONSET SYMPTOMS", "OSDR")          # Delay 4
DR <- calculate_delay(DR, "DATE RECOVERY", "DATE HOSPITALISATION", "HDR")          # Delay 5


# Filter values between 0 to 20
delay1 <- subset(DR, OSSC >= 0 & OSSC <= 20)
delay2 <- subset(DR, OSR >= 0 & OSR <= 20)
delay3 <- subset(DR, OSDH >= 0 & OSDH <= 20)
delay4 <- subset(DR, OSDR >= 0 & OSDR <= 20)
delay5 <- subset(DR HDR >= 0 & HDR <= 20)



# A function to fit multiple probability distributions to a delay dataset.
# Outputs a list of fitted distributions.
fit_distributions <- function(delay_data, delay_col) {
  delay_data <- delay_data[[delay_col]][delay_data[[delay_col]] > 0]
  descdist(delay_data, boot = 1000) # Summary statistics (fitdistrplus)
  
  # Fit different distributions
  list(
    weibull = fitdist(delay_data, "weibull"),
    gamma = fitdist(delay_data, "gamma"),
    lnorm = fitdist(delay_data, "lnorm"),
    llogis = fitdist(delay_data, "llogis"),
    cauchy = fitdist(delay_data, "cauchy")
  )
}

# A function to extract and structure goodness-of-fit statistics for fitted distributions
summarize_gof <- function(fitted_distributions, category) {
  gof_stats <- gofstat(unlist(fitted_distributions, recursive = FALSE))
  data.frame(
    Category = category,
    Kolmogorov_Smirnov = gof_stats$ks,
    Cramer_von_Mises = gof_stats$cvm,
    Anderson_Darling = gof_stats$ad,
    AIC = gof_stats$aic,
    BIC = gof_stats$bic
  )
}

# Delays to analyze (filtered datasets)
filtered_delays <- list(
  delay1 = filter_delay(DR, "OSSC", 0, 20),
  delay2 = filter_delay(DR, "OSR", 0, 20),
  delay3 = filter_delay(DR, "OSDH", 0, 20),
  delay4 = filter_delay(DR, "OSDR", 0, 20),
  delay5 = filter_delay(DR, "HDR", 0, 20)
)

# Fit distributions and summarize results for each delay
gof_results <- list()
for (i in 1:length(filtered_delays)) {
  delay_name <- names(filtered_delays)[i]
  delay_data <- filtered_delays[[i]]
  fitted_distributions <- fit_distributions(delay_data, substr(delay_name, 7, nchar(delay_name)))
  gof_results[[delay_name]] <- summarize_gof(fitted_distributions, delay_name)
}


# Export Goodness-of-Fit Results
# Combine results into a single data frame
combined_gof <- do.call(rbind, gof_results)

# Save the combined results to an Excel file
output_file <- "~/path/to/DelayAnalysisResults.xlsx"
write_xlsx(combined_gof, output_file)


# Sex-Based Analysis
# Analyze delays by gender (Male/Female)
sex_categories <- unique(DR$Sexo)
sex_gof_results <- list()

for (sex in sex_categories) {
  for (i in 1:length(filtered_delays)) {
    delay_name <- names(filtered_delays)[i]
    delay_data <- filtered_delays[[i]] %>% filter(Sexo == sex)
    
    # Fit distributions and append results
    fitted_distributions <- fit_distributions(delay_data, substr(delay_name, 7, nchar(delay_name)))
    sex_gof_results[[paste0(delay_name, "_", sex)]] <- summarize_gof(fitted_distributions, paste(delay_name, "-", sex))
  }
}

# Combine and export goodness-of-fit results by sex
combined_sex_gof <- do.call(rbind, sex_gof_results)
output_file_sex <- "~/path/to/DelayAnalysisBySex.xlsx"
write_xlsx(combined_sex_gof, output_file_sex)

# Combine all data frames into one
combined_M <- rbind(dfM1, dfM2, dfM3, dfM4, dfM5)
fileAICM <- ("~/Path/to/CountryAICSexM.xlsx")
write_xlsx(combined_M, filAICM)







