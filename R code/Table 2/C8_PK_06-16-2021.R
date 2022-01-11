############### Approximate Bayesian calibration of exposure estimates combining contact-based exposure estimates with the biological measurement ###############
######################################################################### Date: 06/16/2021 ######################################################################
####################################################################  Last organized on 01/10/2022 ##############################################################

##### Assumptions
### long-term shrinkage water ingestion rate
### Autocorrelation factor = 0.75
### sigma = 6 * Ct in the likelihood in formula (14) 


library(tidyr)
library(tibble)
library(stringr)
library(naniar)
library(lubridate)
library(Hmisc)
library(mvtnorm)
library(readr)
library(dplyr)
library(readxl)

setwd("~/Desktop/My files in Dropbox/PK Bayesian/data")

#demog091006 <- read_csv("VariableInformation/demog091006.csv")
#chem091006 <- read_csv("VariableInformation/chem091006.csv")
BrkmrReshx1 <- read.table("BrkmrReshx_Feb7_updated_savitzonly.csv", header = F, sep = ",")
BrkmrReshx <- read.table("BrkmrReshx_May3_savitzonly.csv", header = F, sep = ",")
BrkmrWorkhx <- read.table("BrkmrWorkhx_savitzonly.csv", header = F, sep = ",")

Tot_Ingestional_dose <- read_csv("Tot_Ingestional_dose.csv")
Tot_Inhalational_dose <- read_csv("Tot_Inhalational_dose.csv")
Ingestional_dose_resid <- read_csv("Ingestional_dose_resid.csv")
Ingestional_dose_work <- read_csv("Ingestional_dose_work.csv")
de_demographicscleaned <- read_csv("de_demographicscleaned.csv") %>% rename(ID = PARTICIPANTID)
Comm_cohort_C8_meas <- read_excel("Comm_cohort_C8_meas.xls") %>% rename(ID = participantid, C8.measure = C8num)

data <- Comm_cohort_C8_meas %>% 
  left_join(de_demographicscleaned %>% dplyr::select(ID, YEAROFBIRTH, GENDER, age, weight), by = "ID") %>% 
  filter(ID %in% unique(as.numeric(BrkmrReshx[, 1]))) %>% 
  mutate(WT.kg = weight * 0.45) %>% 
  mutate(test.year = as.numeric(str_sub(visitdate, start = 7, end = 10))) %>% 
  mutate(WT.kg = replace(WT.kg, is.na(WT.kg) & age >= 16 & age < 21, 61.5)) %>%  # EPA Exposure Handbook Ch8 page 18
  mutate(WT.kg = replace(WT.kg, is.na(WT.kg) & age >= 21 & age < 30, 67.9)) %>% 
  mutate(WT.kg = replace(WT.kg, is.na(WT.kg) & age >= 30 & age < 40, 70.2)) %>% 
  mutate(WT.kg = replace(WT.kg, is.na(WT.kg) & age >= 40 & age < 50, 72.7)) %>% 
  mutate(WT.kg = replace(WT.kg, is.na(WT.kg) & age >= 50 & age < 60, 73.6))


Body.Weight <- data.frame(Age.group = c("0 to <1", "1 to <2", "2 to <3", "3 to <6", "6 to <11", "11 to <16", "16 to <21",
                                        "21 to <30", "30 to <40", "40 to <50", "50 to <60", "60 to <70", "70 to <80", ">80"),
                          Body.weight = c(6.55, 11.1, 13.2, 17.5, 29.0, 53.3, 61.5, 67.9, 70.2, 72.7, 73.6, 73.9, 69.0, 62.8))

water_consumption0 <- read.table("Matlab code for body weight and water consumption rates over time/water_consumption0.csv", header = FALSE, sep = ",")
waterconsumption_2_11_2011 <- read.table("Matlab code for body weight and water consumption rates over time/waterconsumption_2_11_2011.csv", header = TRUE, sep = ",")


waterconsumption <- waterconsumption_2_11_2011 %>% 
  dplyr::select(BrookMarID, Gender, R0) %>% 
  dplyr::rename(ID = BrookMarID) 

data1 <- data %>% left_join(waterconsumption, by = "ID") %>% rename(birth.year = YEAROFBIRTH)

sum(is.na(data1$R0))  # 2462, 40%  
sum(is.na(data1$WT.kg))  # 0


### Shrinkage estimation
Age <- c("0 to <2", "2 to <16", "16 to <21", "21 to <50", "50+")
mean <- 0.001 * c(438, 576, 995, 1552, 1490)
ninetyfive.p <- 0.001 * c(913, 1092, 1809, 2862, 2615)
sd <- (ninetyfive.p - mean)/qnorm(0.95)
logmean <- log(mean^2 / sqrt(mean^2 + sd^2))
logsd <- sqrt(log(1 + (sd^2 / mean^2)))
table <- data.frame(Age, mean, ninetyfive.p, sd, logmean, logsd)
table


# volume of distribution based on Butenhoff et al. (2004) (L/kg)
mu.volume <- 0.198   # arithmetic mean
sigma.volume <- 0.069   # standard deviation
mu_vd <- log( mu.volume^2 / sqrt(mu.volume^2 + sigma.volume^2) )   # log mean = -1.68
sigma_vd <- sqrt( log(1 + (sigma.volume^2 / mu.volume^2)) )   # log sd = 0.34


mu.halflife <- 2.3  # arithmetic mean
sigma.halflife <- (2.4 - 2.1) / (2 * qnorm(0.975)) # standard deviation
mu_hl <- log( mu.halflife^2 / sqrt(mu.halflife^2 + sigma.halflife^2) ) # log mean = 0.83
sigma_hl <- sqrt( log( 1 + (sigma.halflife^2 / mu.halflife^2) ) ) # log sd = 0.03
# ggplot() + stat_function(fun = dlnorm, args = list(meanlog = mu_hl, sdlog = sigma_hl), colour = "red") + xlim(0, 10)




ids <- data.frame(ID = unique(BrkmrReshx[,1])) %>% filter(ID %in% data$ID)  # 6134

WT <- data1$WT.kg  # body weight 
sum(is.na(WT))

### year-by-year Pharmacokinetic weight matrix
W <- array(0, dim = c(6134, 58, 7), 
           dimnames = list(ids$ID, paste0("Year", 1951:2008), c("yr1", "yr2", "yr3", "yr4", "yr5", "yr6", "test")))



################################################################## R code written by Yachen Zhu ########################################################################
############################################# translated from Matlab code: ExConc_water_BR6_savitzonly_geo_2_work.m  ###################################################

BrkmrReshx1[, 14] <- BrkmrReshx1[, 3] * 100 + 1
BrkmrReshx[, 14:16] <- 0

Covar_BR <- matrix(0, nrow = nrow(BrkmrReshx)-1, ncol = 1)  

for(k in 1:nrow(BrkmrReshx)){
  if(BrkmrReshx[k, 2] == 0){    
    Covar_BR[k, 1] <- which(BrkmrReshx[k, 3] == BrkmrReshx1[, 14])  # find row of participant ID
    BrkmrReshx[k, 14:16] <- BrkmrReshx1[Covar_BR[k, 1], 11:13]   # participant ID
  }
}

UTM <- read.table("~/Desktop/My files in Dropbox/PK Bayesian/data/Zip_unique.csv", header = F, sep = ",")

ED_BR <- matrix(0, nrow = nrow(BrkmrReshx), ncol = 13)
Match_utm_BR <- matrix(0, nrow = nrow(BrkmrReshx), ncol = 1)


ED_BR[which(BrkmrReshx[, 10] > BrkmrReshx[, 8]), 1:3] <- 0 
ED_BR[which(!is.nan(BrkmrReshx[, 10]) & BrkmrReshx[, 10] > BrkmrReshx[, 8]), 1:3] <- 0

# ED_BR[which(is.na(BrkmrReshx[, 8])), 1] <- 1    # If the Exposure End Year is missing, then assign exposure duration = 1
# ED_BR[which(is.na(BrkmrReshx[, 8])), 2] <- BrkmrReshx[, 6]  # If the Exposure End Year is missing, assign Exposure Start Year
# ED_BR[which(is.na(BrkmrReshx[, 8])), 3] <- BrkmrReshx[, 6]  # If the Exposure End Year is missing, assign Exposure End Year = Exposure Start Year

# if pipe year is between start year and end year, 
# then private well concentration is assigned from start year to (pipe year - 1)
# and public water concentration is assigned from pipe year through end year
for ( k in 1:nrow(BrkmrReshx) ) {
  if ( is.nan(BrkmrReshx[k, 10]) || BrkmrReshx[k, 10] <= BrkmrReshx[k, 6] ) { # if Pipe Year is earlier than start year, exposure is from start year to end year
    ED_BR[k, 1] <- BrkmrReshx[k, 8] - BrkmrReshx[k, 6] + 1  # Exposure Duration
    ED_BR[k, 2] <- BrkmrReshx[k, 6]   # Exposure Start Year
    ED_BR[k, 3] <- BrkmrReshx[k, 8]   # Exposure End Year
  }
  #  else if ( !is.nan(BrkmrReshx[k, 10]) & BrkmrReshx[k, 10] > BrkmrReshx[k, 8] ) {  # if Pipe Year is later than the Exposure End YEar
  #    ED_BR[k, 1] <- 0  # Exposure Duration
  #    ED_BR[k, 2] <- 0  # Exposure Start Year
  #    ED_BR[k, 3] <- 0  # Exposure End Year
  #  }
  else if (BrkmrReshx[k, 6] < BrkmrReshx[k, 10] & BrkmrReshx[k, 10] <= BrkmrReshx[k, 8]) { # If Pipe Year is between Exposure Start Year and Exposure End Year
    ED_BR[k, 1] <- BrkmrReshx[k, 8] - BrkmrReshx[k, 10] + 1  # Exposure End Year - Pipe Year + 1
    ED_BR[k, 2] <- BrkmrReshx[k, 10]  # Exposure Start Year = Pipe Year
    ED_BR[k, 3] <- BrkmrReshx[k, 8]   # Exposure End Year = Exposure End Year
  }
  
  ED_BR[k, 4] <- BrkmrReshx[k, 12]  # drinking water district 0-9
  ED_BR[k, 5] <- BrkmrReshx[k, 4]  # X coordinate
  ED_BR[k, 6] <- BrkmrReshx[k, 5]  # Y coordinate
  ED_BR[k, 7] <- BrkmrReshx[k, 11] # Zip code
  ED_BR[k, 8] <- BrkmrReshx[k, 1]  # Participant ID
  ED_BR[k, 9] <- BrkmrReshx[k, 7]  # start month
  ED_BR[k, 10] <- BrkmrReshx[k, 9]  # end month
  
  ED_BR[k, 11] <- BrkmrReshx[k, 14]  # water source (1: public water, 2:private water, 3:bottled water, NAN:don't know)
  ED_BR[k, 12] <- BrkmrReshx[k, 15]  # bottled water start month
  ED_BR[k, 13] <- BrkmrReshx[k, 16]  # bottled water start year
}

colnames(ED_BR) <- c("Exposure Duration", "Exposure Start Year", "Exposure End Year", "Drinking Water District",
                     "X Coordinate", "Y Coordinate", "Zip Code", "Participant ID", "Start Month", "End Month", 
                     "Water Source", "Bottled Water Start Month", "Bottled Water Start Year")


for (k in 1:nrow(BrkmrReshx)){
  # Replace zip centroid for missing geocodes
  if (ED_BR[k, 7] >= 0 & is.nan(ED_BR[k, 5]) & length(which(ED_BR[k, 7] == UTM[, 1])) != 0) {
    Match_utm_BR[k, 1] <- which(ED_BR[k, 7] == UTM[, 1])
    ED_BR[k, 5] <- UTM[Match_utm_BR[k, 1], 2]
    ED_BR[k, 6] <- UTM[Match_utm_BR[k, 1], 3]
  }
}



conc_zipcode1 <- matrix(0, nrow = 1, ncol = 58)
conc_zipcode2 <- matrix(0, nrow = 1, ncol = 58)
conc_zipcode3 <- matrix(0, nrow = 1, ncol = 58)
conc_zipcode4 <- matrix(0, nrow = 1, ncol = 58)



library(readxl)
XI <- read_excel("cellcolumn.xls", col_names = F)
YI <- read_excel("cellrow.xls", col_names = F)

X1 = 439641;  X2 = 452270;  X3 = 443586;  X4 = 448538;  
Y1 = 4347249; Y2 = 4345129; Y3 = 4342429; Y4 = 4348942; 

xi <- min(abs(XI - X1))
r1 <- which(abs(XI - X1) == xi)

yi <- min(abs(YI - Y1))
c1 <- which(abs(YI - Y1) == yi)

xi <- min(abs(XI - X2))
r2 <- which(abs(XI - X2) == xi)

yi <- min(abs(YI - Y2))
c2 <- which(abs(YI - Y2) == yi)

xi <- min(abs(XI - X3))
r3 <- which(abs(XI - X3) == xi)

yi <- min(abs(YI - Y3))
c3 <- which(abs(YI - Y3) == yi)

xi <- min(abs(XI - X4))
r4 <- which(abs(XI - X4) == xi)

yi <- min(abs(YI - Y4))
c4 <- which(abs(YI - Y4) == yi)


ED_BR_XI <- matrix(0, nrow = nrow(BrkmrReshx), ncol = nrow(XI))
ED_BR_YI <- matrix(0, nrow = nrow(BrkmrReshx), ncol = nrow(YI))

for(k in 1:nrow(BrkmrReshx)){
  ED_BR_XI[k, ] <- t(XI - ED_BR[k, 5]) 
  ED_BR_YI[k, ] <- t(YI - ED_BR[k, 6])
}

ED_BR_XI_min <- apply(abs(ED_BR_XI), 1, min)
ED_BR_YI_min <- apply(abs(ED_BR_YI), 1, min)

ED_BR_XI_min_which <- as.matrix(apply(abs(ED_BR_XI), 1, which.min))
ED_BR_YI_min_which <- as.matrix(apply(abs(ED_BR_YI), 1, which.min))

# as.matrix(do.call("rbind", lapply(ED_BR_XI_min_which, length)))

for(i in 1:length(ED_BR_XI_min_which)){
  if(length(ED_BR_XI_min_which[[i]]) == 0){
    ED_BR_XI_min_which[[i]] <- NA
  }
}

ED_BR_XI_1 <- data.frame(min = ED_BR_XI_min, which.min = do.call("rbind", ED_BR_XI_min_which))


for(i in 1:length(ED_BR_YI_min_which)){
  if(length(ED_BR_YI_min_which[[i]]) == 0){
    ED_BR_YI_min_which[[i]] <- NA
  }
}

ED_BR_YI_1 <- data.frame(min = ED_BR_YI_min, which.min = do.call("rbind", ED_BR_YI_min_which))



years_run_max <- 58
first_year <- 51  # two-digit year (deposition start date)
first_year1 <- 1951
EF_IH_BR <- matrix(0, nrow = nrow(BrkmrReshx), ncol = 59)  # Exposure fraction
EF_IH_BR1 <- matrix(0, nrow = nrow(BrkmrReshx), ncol = 58)

for ( k in 1:nrow(BrkmrReshx) ) {
  if ( k < nrow(BrkmrReshx) - 1 & k != 1 ) {
    for ( yearindex in 1:years_run_max ) {
      the_year <- yearindex - 1 + first_year1
      if ( ED_BR[k, 2] < the_year & ED_BR[k, 3] > the_year ) { 
        # if the year is between the start year and end year
        EF_IH_BR[k, yearindex] <- 1
      }
      else if ( ED_BR[k, 2] == the_year & ED_BR[k, 3] == the_year & ED_BR[k, 10] > ED_BR[k, 9] ) { 
        # if the year is same as start and end year, and end month > start month, assign fraction of total months/12
        EF_IH_BR[k, yearindex] <- (ED_BR[k, 10] - ED_BR[k, 9])/12
      }
      else if ( ED_BR[k, 2] == the_year & ED_BR[k - 1, 3] == the_year & ED_BR[k, 9] == ED_BR[k - 1, 10] + 1 ) { 
        # if the year is start year for k and also end year for previous residence and start month is the month after last residence end month
        EF_IH_BR[k, yearindex] <- 1 - (ED_BR[k, 9] - 1) / 12
      } 
      else if ( ED_BR[k, 2] == the_year & ED_BR[k, 2] == ED_BR[k - 1, 3] & ED_BR[k, 9] < ED_BR[k - 1, 10] ) {
        # if the year is start year for k and also end year for previous residence and start month is before the end of previous residence end month
        EF_IH_BR[k, yearindex] <- 0.5
      }
      else if ( ED_BR[k, 2] == the_year ) {
        EF_IH_BR[k, yearindex] <- 1 - ED_BR[k, 9] / 12
      }
      else if ( ED_BR[k, 3] == the_year ) {
        EF_IH_BR[k, yearindex] <- ED_BR[k, 10] / 12
      }
      else {
        EF_IH_BR[k, yearindex] <- 0
      }
    }
  }
  else {
    for ( yearindex in 1:years_run_max ) {
      the_year <- yearindex - 1 + first_year1
      if ( ED_BR[k, 2] < the_year & ED_BR[k, 3] > the_year ) {
        EF_IH_BR[k, yearindex] <- 1
      }
      else if (ED_BR[k, 3] == the_year ) {
        EF_IH_BR[k, yearindex] <- ED_BR[k, 10] / 12
      }
    }
  }
  EF_IH_BR[k, yearindex + 1] <- ED_BR[k, 8]
}


index <- 1
index1 <- 1
for ( k in 1:nrow(ED_BR) ) {
  if ( k < nrow(ED_BR) ) {
    if ( EF_IH_BR[k, 59] == EF_IH_BR[k + 1, 59] ) {
      
    }
    else {
      for ( m in 1:58 ) {
        n <- Matrix::nnzero(EF_IH_BR[index1:index, m])
        if ( n >= 2 & sum(EF_IH_BR[index1:index, m]) > 1 ) {
          EF_IH_BR1[index1:index, m] <- EF_IH_BR[index1:index, m] / sum(EF_IH_BR[index1:index, m])
        }
        else {
          EF_IH_BR1[index1:index, m] <- EF_IH_BR[index1:index, m]
        }
      }
      index1 = index + 1
    }
    index = index + 1
  }
}

EF_IH_BR1[nrow(ED_BR) - 4:nrow(ED_BR), 1:58] <- EF_IH_BR[nrow(ED_BR) - 4:nrow(ED_BR), 1:58]
EF_IH_BR1[nrow(ED_BR), 27:55] <- EF_IH_BR[nrow(ED_BR), 27:55] * 0.5
EF_IH_BR1[nrow(ED_BR) - 1, 27:55] <- EF_IH_BR[nrow(ED_BR) - 1, 27:55] * 0.5





Percentpublic <- read.csv("~/Desktop/My files in Dropbox/PK Bayesian/Pipe Networks/Percentpublic_new.csv", header = F)

## public water concentrations by district
conc_Belpre <- read.table("public water concentration files/conc_Belpre1_new.csv", header = F, sep = ",")
conc_Belpre <- conc_Belpre[1:2,]


conc_LittleHocking <- read.table("public water concentration files/conc_LittleHocking1_new.csv", header = F, sep = ",")
conc_LittleHocking <- conc_LittleHocking[1:2,]


conc_Lubeck <- read.table("public water concentration files/conc_Lubeck1_new.csv", header = F, sep = ",")
conc_Lubeck <- conc_Lubeck[1:2,]


conc_TuppersPlains <- read.table("public water concentration files/conc_TuppersPlains1_new.csv", header = F, sep = ",")
conc_TuppersPlains <- conc_TuppersPlains[1:2, ]


conc_Mason <- read.table("public water concentration files/conc_MasonCounty1_new.csv", header = F, sep = ",")
conc_Mason <- conc_Mason[1:2, ]


conc_Pomeroy <- read.table("public water concentration files/conc_Pomeroy_new.csv", header = F, sep = ",") 
conc_Pomeroy <- conc_Pomeroy[1:2, ]



# Multiply calibration coefficient and add one more column for downstream PSDs
conc_Mason[2, ] <- conc_Mason[2, ] * 0.55
conc_Pomeroy[2, ] <- conc_Pomeroy[2, ] * 0.7
conc_Belpre[2, ] <- conc_Belpre[2, ] * 1.35
conc_LittleHocking[2, ] <- conc_LittleHocking[2, ] * 2
conc_Lubeck[2, 1:38] <- conc_Lubeck[2, 1:38] * 2.1 # For old Lubeck
conc_Lubeck[2, 39:58] <- conc_Lubeck[2, 39:58] * 1.65 # For new Lubeck
conc_TuppersPlains[2, ] <- conc_TuppersPlains[2, ] * 2

conc_all <- cbind(t(conc_Belpre), t(conc_Mason)[, 2], t(conc_Pomeroy)[, 2], t(conc_LittleHocking)[, 2], 
                  t(conc_Lubeck)[, 2], t(conc_TuppersPlains)[, 2]) %>% as.data.frame()
colnames(conc_all) <- c("Year", "Belpre", "Mason", "Pomeroy", "Little Hocking", "Lubeck", "Tuppers Plains")


library(tidyr)
conc_all_log <- conc_all %>% 
  mutate(Belpre = log10(Belpre + 0.001), 
         Mason = log10(Mason),
         Pomeroy = log10(Pomeroy),
         `Little Hocking` = log10(`Little Hocking` + 0.001),
         Lubeck = log10(Lubeck), 
         `Tuppers Plains` = log10(`Tuppers Plains`)) %>% 
  gather(key = "Water District", value = "log concentration", -Year)



GAC_BL <- 2006
GAC_LH <- 2007  # November 2
GAC_LB <- 2007  # June 4
GAC_TP <- 2006
GAC_MC <- 2008
GAC_PO <- 2006



L <- 58
corr <- toeplitz(0.75^(0:(L-1)))   # AR(1) correlation structure
sigma <- matrix(0, L, L)
# diag(sigma) <- rep(sqrt(log(1 + 2^2)), L)



vol <- 0.2366  # L/day, drinking water for one 8oz cup

birth.year <- data1$birth.year
test.year <- data1$test.year

FR <- data.frame( matrix(0, nrow = 6134, ncol = 2008 - min(birth.year, na.rm = T) + 1) )
colnames(FR) <- paste0("Year", min(birth.year, na.rm = T):2008)
FR <- cbind(ID = ids$ID, FR)


# fraction of year lived in the birth year and the lab testing year for each participant
FR.end <-
  difftime(
    as.Date(data1$visitdate, format = "%m/%d/%Y"),
    as.Date(paste0("1/1/", data1$test.year), format = "%m/%d/%Y"),
    unit = "days"
  ) / yearDays(as.Date(as.character(data1$test.year), format = "%Y"))
FR.end <- as.numeric(FR.end)

FR.start <- (data1$age - FR.end) - floor(data1$age - FR.end) 
FR.start <- as.numeric(FR.start)

# fraction of year lived in each year for each participant (equal to 1 in any years between the birth year and testing year)
for(i in 1:nrow(FR)){
  col.start <- birth.year[i] - 1946 + 1    # min(birth.year, na.rm = T) == 1947
  col.end <- test.year[i] - 1946 + 1
  
  if(!is.na(col.start)){
    FR[i, col.start] <- FR.start[i]
    FR[i, col.end] <- FR.end[i]
    FR[i, (col.start+1):(col.end-1)] <- 1
  }
}

FR1 <- FR

# calculate the age in each year for each participant, cumsum of FR
FR1[, -1] <- t(apply(FR[, -1], 1, cumsum))
FR1 <- FR1 %>% dplyr::select(ID, Year1951:Year2008)

BW <- matrix(0, nrow = nrow(FR1), ncol = ncol(FR1))
rownames(BW) <- rownames(FR1)
colnames(BW) <- colnames(FR1)
BW[, 1] <- FR1[, 1]

# Age-specific body weight for the participants based on the EPA Exposure Factors Handbook
for(i in 1:nrow(BW)){
  for(j in 2:ncol(BW)){
    BW[i, j] <- Body.Weight[1, 2] * (FR1[i, j] > 0 & FR1[i, j] < 1) +
      Body.Weight[2, 2] * (FR1[i, j] >= 1 & FR1[i, j] < 2) +
      Body.Weight[3, 2] * (FR1[i, j] >= 2 & FR1[i, j] < 3) +
      Body.Weight[4, 2] * (FR1[i, j] >= 3 & FR1[i, j] < 6) +
      Body.Weight[5, 2] * (FR1[i, j] >= 6 & FR1[i, j] < 11) +
      Body.Weight[6, 2] * (FR1[i, j] >= 11 & FR1[i, j] < 16) +
      Body.Weight[7, 2] * (FR1[i, j] >= 16 & FR1[i, j] < 21) +
      Body.Weight[8, 2] * (FR1[i, j] >= 21 & FR1[i, j] < 30) +
      Body.Weight[9, 2] * (FR1[i, j] >= 30 & FR1[i, j] < 40) +
      Body.Weight[10, 2] * (FR1[i, j] >= 40 & FR1[i, j] < 50) +
      Body.Weight[11, 2] * (FR1[i, j] >= 50 & FR1[i, j] < 60) +
      Body.Weight[12, 2] * (FR1[i, j] >= 60 & FR1[i, j] < 70) +
      Body.Weight[13, 2] * (FR1[i, j] >= 70 & FR1[i, j] < 80) +
      Body.Weight[14, 2] * (FR1[i, j] > 80)
  }
}

# Measured body weight for the participants
BW.measured <- data1 %>% 
  dplyr::select(ID, age, weight, WT.kg) %>% 
  mutate(missing.BW = ifelse(is.na(weight), 1, 0))

# Update the body weight in adulthood (>= age 18) for the participants who have reported their body weight in 2005/2006
for(i in 1:nrow(BW)){
  if(BW.measured[i, ]$missing.BW == 0){
    for(j in 2:ncol(BW)){
      if(FR1[i, j] >= 18) BW[i, j] <- BW.measured$WT.kg[i]
    }
  }
}

FR <- FR[, -1] 
FR <- FR %>% dplyr::select(Year1951:Year2008)
FR1 <- FR1[, -1]

Tot_Inhalational_dose$ID <- data.frame(ID = unique(BrkmrReshx[,1]))$ID 
Inhal <- Tot_Inhalational_dose %>% filter(ID %in% ids$ID) 

Ingestional_dose_work$ID <- data.frame(ID = unique(BrkmrReshx[,1]))$ID 
Inges.work <- Ingestional_dose_work %>% filter(ID %in% ids$ID)




y = read.csv("savitzdata_corrected.csv", header=TRUE);
names(y) = c("participantid","C8_estimate","Exposure_yr","MaternalAge","parity", "MaternalEduc", "MaternalSmokeStat", "Preeclampsia", "Preterm", "TermLBW", "BirthDefect")

y %>% rename(ID = participantid) %>% 
  dplyr::select(ID, Exposure_yr) %>% 
  dplyr::group_by(ID) %>% dplyr::summarize(n = n()) %>% 
  View()  # max == 6

pregnant.yrs <- y %>% rename(ID = participantid) %>% 
  dplyr::select(ID, Exposure_yr) %>% 
  reshape2::dcast(ID ~ "Exposure_yr", value.var = "Exposure_yr", fun.aggregate = function(x) paste(x, collapse = ",")) %>% 
  separate(col = "Exposure_yr", into = c("preg1", "preg2", "preg3", "preg4", "preg5", "preg6"), sep = ",") %>% 
  mutate(preg1 = as.numeric(preg1), preg2 = as.numeric(preg2), preg3 = as.numeric(preg3),
         preg4 = as.numeric(preg4), preg5 = as.numeric(preg5), preg6 = as.numeric(preg6),
         test.year = data1$test.year)


################################################################################# Monte Carlo Simulations ##########################################################################################
l <- 1000
C8.est <- array(dim = c(6134, 8, l))
wts <- matrix(NA, nrow = 6134, ncol = l)
wt <- c()


for( q in 1:l ){
  
  # hyper-prior for the ratio of standard deviation to mean of the PFOA water concentration
  tau <- runif(1, min = 1, max = 10) 
  
  diag(sigma) <- rep(sqrt(log(1 + tau^2)), L)
  
  ## Assuming sigma = tau * mu
  log_mu_Mason <- unname( unlist( log( conc_Mason[2, ] / sqrt(1 + tau^2) ) ) )
  log_mu_Pomeroy <- unname( unlist( log( conc_Pomeroy[2, ] / sqrt(1 + tau^2) ) ) )
  log_mu_Belpre <- unname( unlist( log( conc_Belpre[2, ] / sqrt(1 + tau^2) ) ) )
  log_mu_LittleHocking <- unname( unlist( log( conc_LittleHocking[2, ] / sqrt(1 + tau^2) ) ) )
  log_mu_Lubeck <- unname( unlist( log( conc_Lubeck[2, ] / sqrt(1 + tau^2) ) ) )
  log_mu_TuppersPlains <- unname( unlist( log( conc_TuppersPlains[2, ] / sqrt(1 + tau^2) ) ) )
  
  
  conc_Mason[q+2, ] <- exp( rmvnorm(1, mean = log_mu_Mason, sigma = sigma %*% corr %*% sigma) )
  conc_Pomeroy[q+2, ] <- exp( rmvnorm(1, mean = log_mu_Pomeroy, sigma = sigma %*% corr %*% sigma) )
  conc_Belpre[q+2, ] <- exp( rmvnorm(1, mean = log_mu_Belpre, sigma = sigma %*% corr %*% sigma) )
  conc_LittleHocking[q+2, ] <- exp( rmvnorm(1, mean = log_mu_LittleHocking, sigma = sigma %*% corr %*% sigma) )
  conc_Lubeck[q+2, ] <- exp( rmvnorm(1, mean = log_mu_Lubeck, sigma = sigma %*% corr %*% sigma) )
  conc_TuppersPlains[q+2, ] <- exp( rmvnorm(1, mean = log_mu_TuppersPlains, sigma = sigma %*% corr %*% sigma) )
  
  
  conc_Belpre[q+2, ] <- ifelse(conc_Belpre[1, ] > GAC_BL - 1, 0, conc_Belpre[q+2, ])
  
  conc_LittleHocking[q+2, ] = 0 * ( conc_LittleHocking[1, ] > GAC_LH ) + 
    conc_LittleHocking[q+2, ] * 10 / 12 * ( conc_LittleHocking[1, ] == GAC_LH ) +  # accounting for 10 months of not using GAC filter
    conc_LittleHocking[q+2, ] * ( conc_LittleHocking[1, ] < GAC_LH )
  
  conc_Lubeck[q+2, ] = 0 * ( conc_Lubeck[1, ] > GAC_LB ) + 
    conc_Lubeck[q+2, ] * 5/12 * (conc_Lubeck[1, ] == GAC_LB ) +  # accounting for 5 months of not using GAC filter
    conc_Lubeck[q+2, ] * (conc_Lubeck[1, ] < GAC_LB )
  
  conc_TuppersPlains[q+2, ] = ifelse( conc_TuppersPlains[1, ] > GAC_TP - 1, 0, conc_TuppersPlains[q+2, ] )
  
  conc_Mason[q+2, ] = ifelse(conc_Mason[1, ] > GAC_MC - 1, 0, conc_Mason[q+2, ]) 
  
  conc_Pomeroy[q+2, ] = ifelse(conc_Pomeroy[1, ] > GAC_PO - 1, 0, conc_Pomeroy[q+2, ])
  
  
  
  # Assign PSD Exposure Concentrations for public water use participant
  EC_IN_BR <- matrix(0, nrow = nrow(ED_BR), ncol = ncol(conc_Belpre))
  
  
  # Median bottled water drink start year
  for ( k in 1:nrow(ED_BR) ) {
    for( m in 1:ncol(conc_Belpre) ) {
      
      EC_IN_BR[k, m] = 0 * (ED_BR[k, 2] > conc_Belpre[1, m] |
                              ED_BR[k, 3] < conc_Belpre[1, m] |
                              ED_BR[k, 4] == 0 |
                              ED_BR[k, 4] == 7 |
                              ED_BR[k, 4] == 8 |
                              is.nan(ED_BR[k, 4]) ||
                              (!is.nan(ED_BR[k, 11]) & ED_BR[k, 11] == 2) ) +  # water source: private water
        conc_Belpre[q+2, m] * ( ED_BR[k, 4] == 1 ) + 
        conc_TuppersPlains[q+2, m] * ( ED_BR[k, 4] == 2 ) + 
        conc_LittleHocking[q+2, m] * ( ED_BR[k, 4] == 3 ) + 
        conc_Lubeck[q+2, m] * ( ED_BR[k, 4] == 4 ) + 
        conc_Mason[q+2, m] * ( ED_BR[k, 4] == 5 ) + 
        conc_Pomeroy[q+2, m] * ( ED_BR[k, 4] == 6 ) 
      
    }
  }
  
  
  for ( k in 1:nrow(ED_BR) ) {
    for ( m in 1:ncol(conc_Belpre) ) {
      # if bottled water used and if bottled water start year is before the current residential year
      EC_IN_BR[k, m] <- ifelse(!is.nan(ED_BR[k, 11]) & ED_BR[k, 11] == 3 & ED_BR[k, 13] < conc_Belpre[1, m] + 1, 
                               0,  
                               EC_IN_BR[k, m]) 
    }
  }  
  
  
  # start from here: yearly GW layer 1 PFOA concentrations from MT3DMS output
  for(index in 1:years_run_max){
    
    the_year <- index - 1 + first_year
    
    the_year <- (the_year - 100) * (the_year > 99) + the_year * (the_year <= 99)
    
    mt_file <- "mt3d_"
    
    if(the_year < 41){
      special_case <- c('00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15',
                        '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31',
                        '32', '33', '34', '35')
      mt_file <- paste0(mt_file, special_case[the_year + 1])
    } else {
      mt_file <- paste0(mt_file, the_year)
    }
    
    mt_file <- paste0(mt_file, '.csv')
    
    yearindex <- (the_year + 50) * (the_year < 49) + (the_year - 50) * (the_year >= 49)
    yearindex1 <- (the_year + 2000) * (the_year < 49) + (the_year + 1900) * (the_year >= 49)
    
    file_to_open9 <- paste0('~/Desktop/My files in Dropbox/PK Bayesian/data/mt3d10_csv/', mt_file)
    gw_conc <- read.table(file_to_open9, sep = ",", header = F) %>% t() 
    ## Groundwater concentration from the shallowest layer 1. 235*135 grid cells for each year as the index goes up 58 years 
    
    gw_conc_1 <- matrix(0, nrow = nrow(BrkmrReshx), ncol = 1)
    
    for(k in 1:nrow(BrkmrReshx)){
      gw_conc_1[k, 1] <- gw_conc[ED_BR_XI_1[k, 2], ED_BR_YI_1[k, 2]] 
    }
    
    # 6.1.1 Assign drinking water concentrations for missing (WD = blank), out of area (WD = 0), other water district (WD = 7),
    # and private wells (WD = 8) records from the shallowest layer 1 IF these are inside groundwater model domain from their geocode
    # 6.1.2 OR public water district with the condition that pipe year is later than end year or pipe year is between start year and end year
    
    EC_IN_BR[, yearindex] <- ifelse(
      ( !is.nan(ED_BR[, 5]) & !is.nan(ED_BR[, 6]) & !is.nan(BrkmrReshx[, 10]) &
          ((BrkmrReshx[, 10] > BrkmrReshx[, 8] | (BrkmrReshx[, 6] < BrkmrReshx[, 10] & BrkmrReshx[, 10] < BrkmrReshx[, 8])) & 
             ED_BR[, 2] <= yearindex & ED_BR[, 3] >= yearindex & 
             ED_BR[, 6] >= min(XI[, 1]) & ED_BR[, 6] <= max(XI[, 1]) & 
             ED_BR[, 7] >= min(YI[, 1]) & ED_BR[, 7] <= max(YI[, 1]) ) ) | 
        ED_BR[, 4] %in% c(7, 8, 9, 0) & ED_BR[, 3] > 2006,
      
      gw_conc_1[, 1] *
        ( BrkmrReshx[, 10] > BrkmrReshx[, 8] | (BrkmrReshx[, 6] < BrkmrReshx[, 10] & BrkmrReshx[, 10] < BrkmrReshx[, 8]) & 
            ED_BR[, 2] <= yearindex & ED_BR[, 3] >= yearindex & 
            ED_BR[, 6] >= min(XI[, 1]) & ED_BR[, 6] <= max(XI[, 1]) & 
            ED_BR[, 7] >= min(YI[, 1]) & ED_BR[, 7] <= max(YI[, 1]) ) +
        0 * (ED_BR[, 4] %in% c(7, 8, 9, 0) & ED_BR[, 3] > 2006),  # install GAC in private wells in 2006 if PFOA concentration is greater than 0.5 ug/L
      
      EC_IN_BR[, yearindex])
    
    
    # 6.2 Estimate private well concentrations for addresses that could not be geocoded but have a ZIP code
    
    conc_zipcode1[1, yearindex] = gw_conc[r1, c1]  # Zipcode = 45742
    conc_zipcode2[1, yearindex] = gw_conc[r2, c2]  # Zipcode = 26101
    conc_zipcode3[1, yearindex] = gw_conc[r3, c3]  # Zipcode = 26181
    conc_zipcode4[1, yearindex] = gw_conc[r4, c4]  # Zipcode = 45714
    
    # 6.3 Estimate drinking water concentrations for addresses that could not be geocoded with missing water district code or out of area (WD = 0) but have a ZIP code 
    # -->a weighted average of the public and private waters C8 values for each ZIP code is used   
    
    EC_IN_BR[, yearindex] = 
      ifelse(ED_BR[, 4] == 9 & ED_BR[, 2] <= yearindex1 & ED_BR[, 3] >= yearindex1, 
             (Percentpublic[1, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45701) +
               (Percentpublic[2, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45710) +
               (Percentpublic[3, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45735) +
               (Percentpublic[4, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45743) +
               (Percentpublic[5, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45770) +
               (Percentpublic[6, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45772) +
               (Percentpublic[7, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45776) +
               (Percentpublic[8, yearindex + 1] * conc_all[3, yearindex]) * (ED_BR[, 7] == 45724) +
               (Percentpublic[9, yearindex + 1] * conc_all[3, yearindex]) * (ED_BR[, 7] == 45729) +
               (Percentpublic[10, yearindex + 1] * conc_all[3, yearindex] + (1 - Percentpublic[10, yearindex]) * conc_zipcode1[1, yearindex]) * (ED_BR[, 7] == 45742) +
               (Percentpublic[11, yearindex + 1] * conc_all[3, yearindex]) * (ED_BR[, 7] == 45784) +
               (Percentpublic[12, yearindex + 1] * conc_all[3, yearindex]) * (ED_BR[, 7] == 45786) +
               (Percentpublic[13, yearindex + 1] * conc_all[4, yearindex] + (1 - Percentpublic[13, yearindex]) * conc_zipcode2[1, yearindex]) * (ED_BR[, 7] == 26101) +
               (Percentpublic[14, yearindex + 1] * conc_all[4, yearindex]) * (ED_BR[, 7] == 26133) +
               (Percentpublic[15, yearindex + 1] * conc_all[4, yearindex] + (1 - Percentpublic[15, yearindex]) * conc_zipcode3[1, yearindex]) * (ED_BR[, 7] == 26181) +
               (Percentpublic[16, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25082) +
               (Percentpublic[17, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25106) +
               (Percentpublic[18, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25123) +
               (Percentpublic[19, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25187) +
               (Percentpublic[20, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25241) +
               (Percentpublic[21, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25253) +
               (Percentpublic[22, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25260) +
               (Percentpublic[23, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25287) +
               (Percentpublic[24, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25502) +
               (Percentpublic[25, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25503) +
               (Percentpublic[26, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25515) +
               (Percentpublic[27, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25520) +
               (Percentpublic[28, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25550) +
               (Percentpublic[29, yearindex + 1] * conc_all[6, yearindex]) * (ED_BR[, 7] == 45760) +
               (Percentpublic[30, yearindex + 1] * conc_all[1, yearindex] + Percentpublic[31, yearindex + 1] * conc_all[3, yearindex] +
                  (1 - Percentpublic[30, yearindex] - Percentpublic[31, yearindex]) * conc_zipcode4[1, yearindex]) * (ED_BR[, 7] == 45714) +
               (Percentpublic[32, yearindex + 1] * conc_all[3, yearindex] + Percentpublic[33, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45723) +
               (Percentpublic[34, yearindex + 1] * conc_all[2, yearindex] + Percentpublic[35, yearindex + 1] * conc_all[6, yearindex]) * (ED_BR[, 7] == 45769) +
               (Percentpublic[36, yearindex + 1] * conc_all[2, yearindex] + Percentpublic[37, yearindex + 1] * conc_all[6, yearindex]) * (ED_BR[, 7] == 45771) +
               (Percentpublic[38, yearindex + 1] * conc_all[3, yearindex] + Percentpublic[39, yearindex + 1] * conc_all[2, yearindex]) * (ED_BR[, 7] == 45778) +
               (Percentpublic[40, yearindex + 1] * conc_all[3, yearindex]) * (ED_BR[, 7] == 45712) +
               (Percentpublic[41, yearindex + 1] * conc_all[3, yearindex]) * (ED_BR[, 7] == 45713) +
               (Percentpublic[42, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25239) +
               (Percentpublic[43, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25264) +
               (Percentpublic[44, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25265) +
               (Percentpublic[45, yearindex + 1] * conc_all[5, yearindex]) * (ED_BR[, 7] == 25541),
             EC_IN_BR[, yearindex]
      )
    
  }
  
  EC_IN_BR[is.na(EC_IN_BR)] <- 0
  
  EC_IN_BR2 <- EC_IN_BR * EF_IH_BR1   # end of translated version Matlab code ExConc_water_BR6_savitzonly_geo_2_work.m 
  
  
  ############################################################################# R code written by Yachen Zhu, 2020 #######################################################################
  EC_IN_BR2 <- cbind(BrkmrReshx[, 1], EC_IN_BR2)
  colnames(EC_IN_BR2) <- c("ID",  paste0("Year", 1:58))
  
  ID <- EC_IN_BR2[, 1]
  # year-by-year weighted water concentrations for each participant
  EC_IN_BR3 <- aggregate(EC_IN_BR2, data.frame(ID), sum)[, -2]  # get the column sums by Participant ID
  
  
  EC_IN_BR4 <- EC_IN_BR3 %>% 
    inner_join(data1, by = "ID") 
  
  EC_IN_BR5 <- EC_IN_BR4[, 2:59]
  
  Dose <- data.frame(matrix(0, nrow = nrow(EC_IN_BR4), ncol = 58))
  colnames(Dose) <- colnames(Tot_Ingestional_dose)[-1]
  
  p <- runif(nrow(Dose), 0, 1)  # sampled percentile for the participants
  
  # unit: ug/L * L = ug
  Dose <- (FR * EC_IN_BR5 * qlnorm(p, table[1, 5], table[1, 6]) * 365.25) * (FR1 >= 0 & FR1 < 2) +  # 0-2 years old
    (FR * EC_IN_BR5 * qlnorm(p, table[2, 5], table[2, 6]) * 365.25) * (FR1 >= 2 & FR1 < 16) +  # 2-16 years old
    (FR * EC_IN_BR5 * qlnorm(p, table[3, 5], table[3, 6]) * 365.25) * (FR1 >= 16 & FR1 < 21) +  # 16-21 years old
    (FR * EC_IN_BR5 * qlnorm(p, table[4, 5], table[4, 6]) * 365.25) * (FR1 >= 21 & FR1 < 50) +  # 21-50 years old
    (FR * EC_IN_BR5 * qlnorm(p, table[5, 5], table[5, 6]) * 365.25) * (FR1 >= 50)  # 50-60 years old
  
  
  # for individuals with self-reported water consumption rate, update the doses for them 
  for(j in 1:58){
    Dose[, j] <- ifelse(!is.na(EC_IN_BR4[, 69]) & (EC_IN_BR4[, 69] %in% c("1", "2", "3", "4", "5", "6")),
                        EC_IN_BR5[, j] * runif(1, min = 0, max = vol) * 365.25 * (!is.na(EC_IN_BR4[, 69]) & EC_IN_BR4[, 69] == "1") +   # R0 = 1
                          EC_IN_BR5[, j] * runif(1, min = vol, max = 2*vol) * 365.25 * (!is.na(EC_IN_BR4[, 69]) & EC_IN_BR4[, 69] == "2") +   # R0 = 2
                          EC_IN_BR5[, j] * runif(1, min = 3*vol, max = 4*vol) * 365.25 * (!is.na(EC_IN_BR4[, 69]) & EC_IN_BR4[, 69] == "3") +   # R0 = 3
                          EC_IN_BR5[, j] * runif(1, min = 5*vol, max = 6*vol) * 365.25 * (!is.na(EC_IN_BR4[, 69]) & EC_IN_BR4[, 69] == "4") +   # R0 = 4
                          EC_IN_BR5[, j] * runif(1, min = 7*vol, max = 8*vol) * 365.25 * (!is.na(EC_IN_BR4[, 69]) & EC_IN_BR4[, 69] == "5") +   # R0 = 5
                          EC_IN_BR5[, j] * runif(1, min = 8*vol, max = 10*vol) * 365.25 * (!is.na(EC_IN_BR4[, 69]) & EC_IN_BR4[, 69] == "6"),   # R0 = 6
                        Dose[, j])
  }
  
  
  # 70% water ingestion at home vs. 30% water ingestion at work
  Dose <- (Dose * 0.7 + Inhal[, -1] + Inges.work[, -1] * 0.3) * (Inges.work[, -1] > 0) + (Dose + Inhal[, -1]) * (Inges.work[, -1] == 0)
  # Ingest.work from the output of previous studies may have been multiplied by 30% already
  
  V <- rlnorm(nrow(W), mu_vd, sigma_vd)  # L/kg
  HL <- rlnorm(nrow(W), mu_hl, sigma_hl)  
  # HL <- rlnorm(nrow(W), mu_hl, 0)  # 0 variance constant - sensitivity analysis
  K <- log(2) / HL
  
  for (i in 1:nrow(W)) {
    for (r in 1:7){
      if(!is.na(pregnant.yrs[i, r+1])){
        for (j in (max(birth.year[i], 1951) - 1950):(pregnant.yrs[i, r+1] - 1950)) {
          W[i, j, r] <- ((1 - exp(-K[i])) / (K[i] * V[i] * BW[i, (j+1)])) * exp(-K[i] * (pregnant.yrs[i, r+1] - j - 1950))
        }
      }
    }
  }
  
  # candidate, induced prior, ug/L
  C8.est[ , , q] <- cbind( pregnant.yrs$ID, apply( W[, 1:58, 1] * Dose, 1, sum ), apply( W[, 1:58, 2] * Dose, 1, sum ), 
                           apply( W[, 1:58, 3] * Dose, 1, sum ), apply( W[, 1:58, 4] * Dose, 1, sum ), apply( W[, 1:58, 5] * Dose, 1, sum ), 
                           apply( W[, 1:58, 6] * Dose, 1, sum ), apply( W[, 1:58, 7] * Dose, 1, sum ) )
  colnames(C8.est) <- colnames(pregnant.yrs)
  
  # add the background serum concentration
  for (i in 1:nrow(W)) {
    for (r in 1:7){
      if(!is.na(pregnant.yrs[i, r+1])){
        C8.est[i, r+1, q] <- (C8.est[i, r+1, q] + 0.11 * pregnant.yrs[i, r+1, q]) * (pregnant.yrs[i, r+1, q] < 1999) +  # Year < 1999
          (C8.est[i, r+1, q] + 5.2 - 0.33*(pregnant.yrs[i, r+1, q] - 1999)) * (pregnant.yrs[i, r+1, q] >= 1999)  # Year >= 1999
      }
    }
  }
  
  
  wts[, q] <-  -1/2*(EC_IN_BR4$C8.measure - C8.est[, 8, q])^2 / (6*EC_IN_BR4$C8.measure)^2  
  wt[q] <- sum(wts[, q], na.rm = T) 
  
  print( paste("This is iteration", q) )
  print( Sys.time() )
  
}



# wt <- apply(wts, 2, sum)

# turn it from an Array into a List
# C8.est.List <- lapply(seq(dim(C8.est)[3]), function(x) C8.est[ , , x])

# mean(wt)  # mean center of q iterations
# wt1 <- exp(wt-mean(wt))

# SIR resampling
# SIR <- sample(x = C8.est.List, size = 1000, replace = T, prob = exp(wt) )
# length(unique(SIR))   # 2


############################################################### Approximate Bayesian Computation ####################################################################

# reject the sample when the weight is 0, select the sample when the weight is not 0
C8.est.abc.1 <- C8.est[ , , exp(wt) != 0]  # acceptance rate = 44.2%
dim(C8.est.abc.1)
C8.est.abc <- C8.est[ , , exp(wt) > 1e-300]  # acceptance rate = 42.9%
dim(C8.est.abc)

# posterior mean
C8.est.abc.mean <- apply(C8.est.abc, 1:2, mean)   

pregnant.yrs.long <- gather(pregnant.yrs, preg, year, -ID)
C8.est.abc.mean.long <- gather(data.frame(C8.est.abc.mean), preg, C8.est.abc, -ID)
C8.long <- pregnant.yrs.long %>% 
  left_join(C8.est.abc.mean.long, by = c("ID", "preg")) %>% 
  filter(!is.na(year)) %>% 
  filter(preg != "test.year") %>% 
  rename(participantid = ID, Exposure_yr = year)

C8.test <- pregnant.yrs.long %>% 
  left_join(C8.est.abc.mean.long, by = c("ID", "preg")) %>% 
  filter(!is.na(year)) %>% 
  filter(preg == "test.year") %>% 
  rename(Exposure_yr = year) %>% 
  left_join(Comm_cohort_C8_meas, by = "ID") %>% 
  rename(participantid = ID)

cor(C8.test$C8.est.abc, C8.test$C8.measure) # 0.66
cor(C8.test$C8.est.abc, C8.test$C8.measure, method = "spearman") # 0.67
mean((C8.test$C8.measure - C8.test$C8.est.abc)^2)  # 7057.48
mean(abs(C8.test$C8.measure - C8.test$C8.est.abc)) # 27.22

y1 <- y %>% left_join(C8.long %>% select(-preg), by = c("participantid", "Exposure_yr")) %>% distinct()
cor(y1$C8_estimate, y1$C8.est.abc, method = "spearman") # 0.34


############################################################################### GEE models ###########################################################################

library(geepack)
library(splines)
source("http://www.biostat.jhsph.edu/~pmurakam/mysummaries.R")


### First strategy to compute AOR: epidemiological analysis using the average PFOA concentrations across all accepted samples
y1$logexp = log(y1$C8.est.abc)
q = quantile(y1$logexp)
iqr = q[4] - q[2]
gee4 = geese(Preeclampsia ~ I(logexp/iqr) + ns(Exposure_yr, 4) + ns(MaternalAge, 3) + parity + MaternalEduc + MaternalSmokeStat, 
             family = binomial(link = "logit"), id = participantid, corstr = "exchangeable", data = y1)
g = geese.mysummary(gee4, alpha = .05, dig = 3, p.dig = 4)
(r = c(g[2, 3], g[2, 4], g[2, 5]))  
# 1.17, 95% CI: [0.92, 1.49]



### Second strategy to compute AOR: analyses based on all accepted samples and then took an average of the effect estimates
Y <- gee4 <- g <- list()
r <- matrix(NA, nrow = dim(C8.est.abc)[3], ncol = 5)
colnames(r) <- c("estimate", "robust", "Odds.Ratio", "OR.lower", "OR.upper")

for(i in 1:dim(C8.est.abc)[3]){     # 429
  
  C8.est.abc.long <- gather(data.frame(C8.est.abc[,,i]), preg, C8.est.abc, -ID) %>% 
    left_join(pregnant.yrs.long, by = c("ID", "preg")) %>% 
    filter(!is.na(year)) %>% 
    filter(preg != "test.year") %>% 
    select(-preg)
  
  Y[[i]] <- y %>% rename(ID = participantid, year = Exposure_yr) %>% 
    inner_join(C8.est.abc.long, by = c("ID", "year")) %>% distinct()
  
  Y[[i]]$logexp = log(Y[[i]]$C8.est.abc)
  
  q = quantile(Y[[i]]$logexp)
  iqr = q[4] - q[2]
  gee4[[i]] = geese(Preeclampsia ~ I(logexp/iqr) + ns(year, 4) + ns(MaternalAge, 3) + parity + MaternalEduc + MaternalSmokeStat, 
                    family = binomial(link = "logit"), id = ID, corstr = "exchangeable", data = Y[[i]])
  g[[i]] = geese.mysummary(gee4[[i]], alpha = .05, dig = 3, p.dig = 4)
  r[i, ] = c(g[[i]][2, 1:5])
}

exp( mean(r[, 1]) )
exp( mean(r[, 1]) - qnorm(0.975) * sqrt(mean(r[, 2]^2) + var(r[, 1])) )  # computed the variance using the law of total variance
exp( mean(r[, 1]) + qnorm(0.975) * sqrt(mean(r[, 2]^2) + var(r[, 1])) )

# 1.14, 95% CI: [0.86, 1.52]
