#----- Code to simulate & fit dengue constant-in-time FOI (Region) catalytic model
# Produces the analyses and code used for the paper (Dengue Epidemiology and Transmission Intensity 
# across Panama during 2000-2024)

# Cite as:                                                                                
# Quijada M et al, Dengue Epidemiology and Transmission Intensity across Panama           
# during 2000-2024                                                                        
# Pre-print available at: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5461361%20  


#Data format expected (wide): Region, District, Age-group, Sex, 2000,2001,...,2024
#Model: region-level constant FOI (lambda), reporting prob rho and gamma.

# Before running: edit USER CONFIG (file paths, region_to_run, cmdstan dir, stan model path)


#Libraries
library(ggplot2)
library(tidyr)
library("rstan")
library(cmdstanr)
library(readxl)
library(dplyr)
library(cowplot)
library(loo)
library(readr)       
library(purrr)

#Read data from excel
Cases <- readRDS("path/to/DengueCases.xlsx") #Dengue Cases by Region 
Pop <- readRDS("path/to/PopRegion2000_24.xlsx") #Population by Region 


#Cases
# Convert wide Year columns into long format
# Each Year becomes a category under "Year"
# The corresponding case counts are placed in the "Cases" column
Casescolumn <- Cases %>%
  pivot_longer(cols = matches("^\\d{4}$"),   # new column with the Year
               names_to = "Year", values_to = "Cases",  # new column with the case counts
               names_transform = list(Year = as.integer))  # <-- ensures Year is integer

# Aggregate to Year x Age-group x Region
Casescolumn <- Casescolumn %>%
  group_by(`Age-group`, Year, Region) %>%
  summarise(Cases = sum(Cases),.groups = 'drop')


# Average cases across years for each Region × Age-group
#    (Mean over all years after the aggregation above)
CasesT <- Casescolumn  %>%
  group_by(`Age-group`, Region) %>%
  dplyr::summarize(Cases = mean(Cases, na.rm = TRUE))


# Reshape to wide format:
  #Each Age-group becomes its own column
temp <- CasesT %>%
  pivot_wider(names_from = `Age-group`,
              values_from = Cases)

# Remove NOT region categories
temp <- temp %>%
  filter(!Region %in% c("DESCONOCI", "IMPORTADO", "EXTRANJERO", "SIN DEFINIR"))
temp[is.na(temp)] <- 0 #Replace any remaining NA value



#Population
#    Reshape population data from wide (age columns) to long
#    Each age group becomes a row instead of a column
Pop24 <- Pop %>%
  pivot_longer(cols = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "over80"),
               names_to = "Age-group", values_to = "Pop")

#Aggregate population by Year × Age-group × Region
Pop24 <- Pop24 %>%
  group_by(`Age-group`, Year, Region) %>%
  summarise(Pop = sum(Pop)) %>%
  ungroup()

#Compute mean population across years for each Region × Age-group
Popmean <- Pop24 %>%
  group_by(`Age-group`, Region) %>%
  summarise(Pop = as.integer(mean(Pop)),.groups = 'drop')

# Convert long population table into wide format
pop_mat <- Popmean %>%
  pivot_wider(names_from = `Age-group`,
              values_from = Pop) %>%
  ungroup


##DENGUE INTENSITY
#Model Setting
age <- 0:108
amin <- c(1,seq(11,81,10)) #Age minimum
amax <- c(seq(10,80,10), 109) #Age maximum
Location <- unique(temp$Region)


# Prepare Stan data list
  data <- list( nA= 9, # Age structure (9 groups)
               nL = length(Location), 
               max_age = 109,
               cases = as.matrix(temp[, -1]),   # int matrix (regions x age groups)
               pop = as.matrix(pop_mat[, -1]),   # numeric matrix (regions x age groups)
               age=age,
               ageLims=rbind(amin,amax)) # 2 x nA integer matrix of age bounds
  

  #Compile and fit the model
  check_cmdstan_toolchain(fix=T)
  set_cmdstan_path("/Users/marioquijada/.cmdstan/cmdstan-2.33.1")
  setwd("~StanCode")
  mod <- cmdstan_model("Constant-in-time by Region.stan", pedantic=T)
  
  fit <- mod$sample(data=data, 
                    chains=4, 
                    parallel_chains=4, 
                    iter_sampling=10000, 
                    refresh=1000, 
                    iter_warmup=5000)

  #output
  stanfit <- rstan::read_stan_csv(fit$output_files())
  saveRDS(stanfit, "/~path/to/your_file.rds")  #Save CmdStanMCMC object
  
  # Check convergence 
  trace <- traceplot(stanfit, pars = c('rho', 'gamma', 'lam_H'))
  ggsave(
    "~/path/to/your_file.png",
    plot = trace,
    width = 20,
    height = 16,
    units = "in"
  )
  
  #check convergence and quantiles for region-level parameters
  chains <- rstan::extract(stanfit)
  rho <- quantile(chains$rho, c(0.5, 0.025, 0.975)))
  gamma <- quantile(chains$gamma, c(0.5, 0.025, 0.975)))  
  lam <- quantile(chains$lam_H, c(0.5, 0.025, 0.975)))  
  kids <- quantile(chains$kids, c(0.5, 0.025, 0.975))

  # Export parameter 
  gamma <- t(gamma)
  lam <- t(lam)
  rho <- t(rho)
  
  #Arrange and label by region
  rownames(gamma) <- Location
  rownames(lam) <- Location
  rownames(rho) <- Location
  
  # Convert to data frame 
  gamma <- as.data.frame(gamma)
  lam <- as.data.frame(lam)
  rho <- as.data.frame(rho)
  kids <- as.data.frame(kids)

  # Export parameter summaries
  write.csv(rho, "~/path/to/your_file.csv")
  write.csv(gamma, "~/path/to/your_file.csv")
  write.csv(lam_H, "~/path/to/your_file.csv")
  write.csv(kids, "~/path/to/your_file.csv")


  
  #fit
  #age groups 
  ageG <- c('0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', 'over80')
  fit <- data.frame(
    location = rep(Location, each = length(ageG)),
    age = rep(ageG, times = length(Location)),
    cases = as.vector(t(data$cases))
  )
  
  fit[, c('pred', 'ciL', 'ciU')] <- NA
  
  for (l in 1:data$nL) for(a in 1:data$nA) {
    fit[fit$location== Location[l] & fit$age == ageG[a],4:6] <- quantile(chains$Ecases[, l, a], c(0.5, 0.025, 0.975), na.rm = TRUE) 
  }
  
  # Convert 'age' and regions' columns to factors for plotting
  fit$age <- factor(fit$age, levels = unique(fit$age))
  fit$location <- factor(fit$location, levels = Location)
  
  #Plot observed (points) vs predicted (line/ribbon) by region
  fit_plot <- ggplot(fit, aes(x = age, y = cases)) +
    geom_point() +
    geom_line(aes(y = pred, group = location), color = 'blue') +
    geom_ribbon(aes(ymin = ciL, ymax = ciU, group = location), fill = 'blue', alpha = 0.2) +
    ylab('Cases') +
    xlab('Age Group') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~ location, scales = "free_y")
  
  print(fit_plot)
  
  #SavePlot
  ggsave(
    "~/path/to/your_file.png",
    plot = fit_plot,
    width = 12,
    height = 8,
    units = "in"
  )
  
  
  
  


  
  
  
  
  
  
  
  
  
  
  
  
 
