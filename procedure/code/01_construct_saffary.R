#------------------------------------------------------------------------#
# Program: 01_construct_saffary                                           #
# Description: Reads in input data and imputes PCP variable to create    #
#              analysis data file for reproducing Saffary et al. paper   #
# Authors: Sarah Bardin and Josh Gilman                                  #
# Date: 3-21-2021                                                        #
#------------------------------------------------------------------------#

##-----------------------##
#----- LOAD LIBRARIES ----#
##-----------------------##

rm(list=ls())
library(tidyverse)
library(readxl)
library(tidycensus)
library(here)
library(sf)

##-----------------------##
#------ READ IN DATA -----#
##-----------------------##

#---County Health Rankings Data---#
chr.pcp <- read_xlsx((here("data", "raw", "public", "2020 County Health Rankings Data - v2.xlsx")), 
                     sheet = "Ranked Measure Data", range =("A2:DF3195"))

chr.pcp <- chr.pcp %>%
  transmute(GEOID = FIPS,
            PCP = `Primary Care Physicians Rate`)
head(chr.pcp)

#---Saffery Data---#
orig.data <- read_xlsx(here("data", "raw", "public", "data_22.05.2020.xlsx"))

#-Clean up variable names-#
names(orig.data) <- substr(names(orig.data),1,regexpr(",",names(orig.data)))
names(orig.data) <- substr(names(orig.data),1,nchar(names(orig.data))-1)
str(orig.data)

#---Boundary Data---#
boundaries <- get_estimates(geography = "county", 
                            product = "population", 
                            year = 2018,                  ## Saffary used 2018 boundary data
                            geometry = TRUE)

boundaries <- boundaries %>% 
  filter(variable == "POP") 

str(boundaries)
str(orig.data)
str(chr.pcp)

##-----------------------##
#------ MERGE DATA  -----#
##-----------------------##
int.data <- left_join(orig.data, chr.pcp, by = "GEOID")
int.data <- int.data %>%
  dplyr::select(GEOID,
                ICUBEDS, 
                PCP,
                DIABETS, 
                OBESITY,
                BLACK, 
                HISPANC, 
                WHITE, 
                UNINSRD, 
                VACCNTN,
                CASS100, 
                DEATH100,
                CASES,
                DEATHS,
                STATENM) %>% 
  dplyr::mutate(PCP_impute_zero = ifelse(is.na(PCP), 0, PCP)) %>% # impute all NA values with 0
  dplyr::mutate(PCP_impute_mean = ifelse(is.na(PCP), 54.5, PCP)) # impute all NA values with mean

str(int.data)
ds_full <- inner_join(boundaries, int.data, by = "GEOID") 
ds_full <- st_transform(ds_full, 5070)
ds_full <- ds_full %>% filter(!is.na(STATENM))  ## Results in 3,142 obs which matches article
                                                ## NOTE: ds_full is used only for descriptive stats in Table 1
summary(ds_full)  

##-----------------------##
#------ EXPORT DATA  -----#
##-----------------------##
st_write(ds_full,here("data", "derived", "ds_full.shp"), delete_dsn=TRUE)

ds_contiguous <- ds_full %>% filter(substr(GEOID,1,2) != "02" & substr(GEOID,1,2) != "15") 
st_write(ds_contiguous,here("data", "derived", "ds_contiguous.shp"), delete_dsn=TRUE)

ds_pcp <- ds_contiguous %>% filter(!is.na(PCP))
st_write(ds_pcp,here("data", "derived", "ds_pcp.shp"), delete_dsn=TRUE)

