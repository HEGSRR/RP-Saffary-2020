#------------------------------------------------------------------------#
# Program: 02_analysis_saffary                                           #
# Description: Reads in analysis data and runs global/local univariate   #
#              and bivariate moran's I analyses and produces report figs #
# Authors: Sarah Bardin and Josh Gilman                                  #
# Date: 3-21-2021                                                        #
#------------------------------------------------------------------------#

##-----------------------##
#----- LOAD LIBRARIES ----#
##-----------------------##
set.seed(1)
rm(list=ls())

library(sf)
library(here)
library(tidyverse)
library(spdep)
library(spatialreg)
library(rgdal)
library(DCluster)
library(maptools)
library(maps)
library(patchwork)
library(ggpubr)
library(fdrtool)
output_maps <- here("results", "maps")
output_figures <- here("results", "figures") #follows R&R template better

##------------------------------------##
#--  DEFINE SPATIAL WEIGHTS MATRICES --#
##------------------------------------##

#--Main Analysis Data--#
ds_contiguous <- st_read(here("data", "derived", "ds_contiguous.shp"))
ds_contiguous <- ds_contiguous[-c(329,1168,1816),] ## remove cases with no neighbors

#--Define Contiguous Weighting Scheme--#
# Create symmetrical, row standardized queen weight matrix w/ first order contiguity
nb <- poly2nb(ds_contiguous, queen = TRUE)
lw <- nb2listw(nb, style = "B", zero.policy = TRUE) ## there are empty neighbor sets, so used zero.policy option
W  <- as(lw, "symmetricMatrix") ## make symmetrical
W <- as.matrix(W/rowSums(W)) ## row standardize
W[which(is.na(W))] <- 0 ## assign NA to zero

#--PCP Analysis Data--#
# NOTE: Because the PCP variable contains missing data, we use the data set which filters out these NA values
#       and perform any PCP analysis using this smaller dataframe.
ds_pcp <- st_read(here("data", "derived", "ds_pcp.shp"))

#--Define PCP Weighting Scheme--#
# Create symmetrical, row standardized queen weight matrix w/ first order contiguity
nb2 <- poly2nb(ds_pcp, queen = TRUE)
lw2 <- nb2listw(nb2, style = "B", zero.policy = TRUE) ## there are still empty neighbor sets
## so used zero.policy option
W2 <- as(lw2, "symmetricMatrix") ## make symmetrical
W2 <- as.matrix(W2/rowSums(W2)) ## row standardize
W2[which(is.na(W2))] <- 0 ## assign NA to zero


##-----------------------------------------##
#- STEP 1) Calculate Table 1 Summary Stats -#
##-----------------------------------------##
# Note: include non-contiguous US counties in this part of analysis because original authors did
ds_full <- st_read(here("data", "derived", "ds_full.shp"))
summary(ds_full)  ## Need to report the mean, median, and IQR for all variables for Table 1

##---------------------------------------------##
#- STEP 2) Perform Global Univariate Moran's I -#
##---------------------------------------------##
#--DEATH RATE--#
moran.test(as.numeric(ds_contiguous$DEATH10), lw, zero.policy = TRUE)

#--CASE RATE--#
moran.test(as.numeric(ds_contiguous$CASS100), lw, zero.policy = TRUE)

##---------------------------------------------##
#- STEP 3a) Perform Local Univariate Moran's I  -#
##---------------------------------------------##

#---------------#
#-- CASE RATE --#
#---------------#
ds_contiguous$localI_case <- localmoran(as.numeric(ds_contiguous$CASS100), lw, zero.policy = TRUE)[,4]

#----Plot Results----#

# set breaks to critical z values
breaks <- c(-100, -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, 100)

# findInterval() assigns ranks to the zvalues based on which bin the z values would fall into 
# where the bins are broken up by the "breaks" variable created above
ds_contiguous$localI_case_sigbreaks <- findInterval(ds_contiguous$localI_case, breaks, all.inside = TRUE)

# Identify high and low clusters using the interaction between the weights matrix and the standardized cases
cases_z <- scale(ds_contiguous$CASS100)[,1]
patterns <- as.character( interaction(cases_z > 0, W%*%cases_z > 0) )
patterns <- patterns %>%
  str_replace_all("TRUE","High") %>%
  str_replace_all("FALSE","Low")

# Replace non-significant I statistics with value of "Not significant"
patterns[ds_contiguous$localI_case_sigbreaks == 4] <- "Not significant"
ds_contiguous$patterns <- patterns

# Create factor version of patterns variable and assign labels for mapping
ds_contiguous$patterns2 <- factor(ds_contiguous$patterns,
                        levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                        labels=c("High  - High", "High - Low", "Low - High","Low - Low", "Not significant"))

theme_set(theme_void())
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) # for highlighting state boundaries as

cases <- ggplot() +
  geom_sf(data=ds_contiguous, aes(fill=patterns2), color="NA") +
  geom_sf(data = states, fill = "NA") +
  scale_fill_manual(values = c("red", "pink", "light blue", "grey80", "blue")) + ## SWITCHE "blue" and "grey08'
  guides(fill = guide_legend(title="Patterns")) + theme(
    legend.position = "bottom"
)
#cases

## ERROR: "Not significant" gets filled in as blue, not grey: we found that the Low-Low are not counted in the 
## levels as there are none (the solution commented line 117) 
# unique(ds_contiguous$patterns2)

#----------------#
#-- DEATH RATE --#
#----------------#
ds_contiguous$localI_death <- localmoran(as.numeric(ds_contiguous$DEATH10), lw, zero.policy = TRUE)[,4]

#----Plot Results----#
ds_contiguous$localI_death_sigbreaks <- findInterval(ds_contiguous$localI_death, breaks, all.inside = TRUE)

# Identify high and low clusters using the interaction between the weights matrix and the standardized cases
deaths_z <- scale(ds_contiguous$DEATH10)[,1]
patterns <- as.character( interaction(deaths_z > 0, W%*%deaths_z > 0) )
patterns <- patterns %>%
  str_replace_all("TRUE","High") %>%
  str_replace_all("FALSE","Low")

patterns[ds_contiguous$localI_death_sigbreaks == 4] <- "Not significant"
ds_contiguous$patterns <- patterns

# Create factor version of patterns variable and assign labels for mapping
ds_contiguous$patterns2 <- factor(ds_contiguous$patterns,
                        levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                        labels=c("High  - High", "High - Low", "Low - High","Low - Low", "Not significant"))

deaths <- ggplot() +
  geom_sf(data=ds_contiguous, aes(fill=patterns2), color="NA") +
  geom_sf(data = states, fill = "NA") +
  scale_fill_manual(values = c("red", "pink", "light blue", "grey80", "blue")) + #Switch "grey80" and "blue"
  guides(fill = guide_legend(title="Patterns")) + theme(
    legend.position = "bottom"
)
#deaths

## ERROR: "Not significant" gets filled in as blue, not grey: we found that the Low-Low are not counted in the 
## levels as there are none (the solution commented line 153) 
# unique(ds_contiguous$patterns2)


#-----------------------------------------------------------------------------------#
#-- STEP 3b) Multiple comparison corrections for Cases/Deaths univariate Moran's I--#
#-----------------------------------------------------------------------------------#
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) # for highlighting state in maps

#----Define function for performming adjustments----#
# NOTE: function takes the variable, dataframe, weight matrices, and multiple comparisons correction (as a string)

univariate_mc <- function(variable, df, WM1, WM2, mc_approach) {
  # quosures https://stackoverflow.com/questions/43438001/how-to-pass-column-names-into-a-function-dplyr
  variable <- enquo(variable)
  var1 <- quo_name(variable)
  
  var_vec <- df %>%
    dplyr::pull(var1)
  var_vec
  
  # Store Moran's I values and associated p-values
  df$localI_case <- localmoran(as.numeric(var_vec), WM1, zero.policy = TRUE)[,4]
  cases_p_values <- localmoran(as.numeric(var_vec), WM1, zero.policy = TRUE)[,5]
  
  if (mc_approach == "bonferroni") {
    #----Bonferroni Adjustment----#
    cases_p <- p.adjust(cases_p_values, method = "bonferroni")
    cases_sig <- (cases_p < 0.05)
    df$cases_sig <- cases_sig
    
  } else {
    #----FDR Adjustments----#
    cases_p <- fdrtool(cases_p_values, statistic = "pvalue", plot = FALSE)
    cases_p <- unlist(cases_p[1])
    cases_sig <- (cases_p < 0.05)
    df$cases_sig <- cases_sig
  }
  
  #----Plot Results----#
  
  # Identify high and low clusters using the interaction between the weights matrix and the standardized cases
  cases_z <- scale(df$CASS100)[,1]
  patterns <- as.character( interaction(cases_z > 0, WM2%*%cases_z > 0) )
  patterns <- patterns %>%
    str_replace_all("TRUE","High") %>%
    str_replace_all("FALSE","Low")
  
  patterns[df$cases_sig == FALSE] <- "Not significant"
  df$patterns <- patterns
  
  df$patterns2 <- factor(df$patterns,
                         levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                         labels=c("High  - High", "High - Low", "Low - High","Low - Low", "Not significant"))
  
  colors <- c("High  - High" = "red", "High - Low" = "pink", "Low - High" = "light blue", "Low - Low" = "dark blue", "Not significant" = "grey80")
  
  theme_set(theme_void())
  g <- ggplot() +
    geom_sf(data=df, aes(fill=patterns2), color="NA") +
    geom_sf(data = states, fill = "NA") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(title="Patterns")) + ggtitle(paste("  ", mc_approach, " ", var1)) + theme(
      legend.position = "bottom"
    )
  return(g)
}
  
#----Perform multiple comparison adjustments for univariate LISA----#
cases_bonferroni <- univariate_mc("CASS100", ds_contiguous, lw, W, "bonferroni")
cases_fdr <- univariate_mc("CASS100", ds_contiguous, lw, W, "fdr")
deaths_bonferroni <- univariate_mc("DEATH10", ds_contiguous, lw, W, "bonferroni")
deaths_fdr <- univariate_mc("DEATH10", ds_contiguous, lw, W, "fdr")
  
##-------------------------------------------------------##
#- STEP 4) Perform Global and Local Bivariate Moran's I  -#
##-------------------------------------------------------##

#-- Takes vectors of predictors/responses to be tested with Bivariate Moran's I --#
variables_func <- function(predictors, responses, df, df_states, WM, mc_approach) {
    ggplots <- list()
    count = 1
    for (i in predictors) {
      for (j in responses) {
        ggplots[[count]] <- individual_tests(i,j,df, df_states, WM, mc_approach)
        count = count + 1
      }
    }
    return(ggplots)
}
  
#-- Performs individual Bivariate Moran's I tests for x,y pairs --#
individual_tests <- function(predictor, response, df, df_states, WM, mc_approach) {
    
    x <- df %>%
      dplyr::pull(predictor)
    
    y <- df %>%
      dplyr::pull(response)
    
    #-- PERFORM BIVARIATE MORAN'S I --#
    m <- moran_I(x, y, WM)
    
    #global bivariate
    (global_moran <- m[[1]][1])
    
    # Local values
    m_i <- m[[2]] # calculated Moran's I value for each county.
    
    # local simulations
    local_sims <- simula_moran(x, y, WM)$local_sims # returns matrix. 999 simulations (columns) for each county (rows).
    
    # global pseudo p-value
    global_sims <- simula_moran(x, y, WM)$global_sims
    
    # Proportion of simulated global values that are higher (in absolute terms) than the actual index
    moran_pvalue <- sum(abs(global_sims) > abs( global_moran )) / length(global_sims)
    
    
    # Logic to deal with no correction, bonferroni, or fdr
    if (mc_approach == "none") {
      # Identifying the significant values
      alpha <- .05  # for a 95% confidence interval; bonferroni adjustment would be dividing this value by 3108.
      probs <- c(alpha/2, 1-alpha/2)
      intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs))) # t() transposes matrix; apply(X, margin, function) applies function to margins of matrix;
      sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] ) # Is the Moran's I value outside the intervals (significant); row number = total # of counties.
    }  else {
      
      n = nrow(WM)
      p_vec <- numeric(n) # allocate p-value vector
      
      # calculate local pseudo p-value
      for (i in 1:n) {
        p_vec[[i]] <- sum(abs(local_sims[i,]) > abs(m_i[i])) / length(local_sims[i,])
      }
      
      if (mc_approach == "bonferroni") {
        # Stats package
        # adj_p_vec <- p.adjust(p_vec, method = "bonferroni")
        # sig <- (adj_p_vec < 0.05)
        
        # Custom correction (alpha divided by number of comparisons (3105))
        alpha <- 0.000016103  # for a 95% confidence interval; bonferroni adjustment would be dividing this value by 3108.
        probs <- c(alpha/2, 1-alpha/2)
        intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs))) # t() transposes matrix; apply(X, margin, function) applies function to margins of matrix;
        sig <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] ) # Is the Moran's I value outside the intervals (significant); row number = total # of counties.
        
      }
      else {
        # Stats package
        # adj_p_vec <- p.adjust(p_vec, method = "fdr")
        
        # fdrtools package
        adj_p_vec <- fdrtool(p_vec, statistic = "pvalue", plot = FALSE)
        adj_p_vec <- unlist(adj_p_vec[1])
        sig <- (adj_p_vec < 0.05)
      }
    }
    
    df$sig <- sig
    
    # Plotting
    g <- vis_tests(x, y, predictor, response, df, df_states, WM, mc_approach)
}
  
#-- LISA visualizations --#
vis_tests <- function(x, y, predictor, response, df, df_states, WM, mc_approach) {
    
    # Identifying the LISA clusters
    xp <- scale(x)[,1]
    yp <- scale(y)[,1]
    
    patterns <- as.character( interaction(xp > 0, WM%*%yp > 0))
    
    patterns <- patterns %>%
      str_replace_all("TRUE","High") %>%
      str_replace_all("FALSE","Low")
    
    patterns[df$sig==0] <- "Not significant" # not sure what this is doing.
    # print(patterns)
    df$patterns <- patterns
    
    # Rename LISA clusters
    df$patterns2 <- factor(df$patterns, levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                           labels=c("High  - High", "High - Low", "Low - High","Low - Low", "Not significant"))
    
    ### PLOT
    colors <- c("High  - High" = "red", "High - Low" = "pink", "Low - High" = "light blue", "Low - Low" = "dark blue", "Not significant" = "grey80")
    
    if (mc_approach == "none") {
      g <- ggplot() +
        geom_sf(data=df, aes(fill=patterns2), color="NA") +
        geom_sf(data = df_states, fill = "NA") +
        scale_fill_manual(values = colors) +
        guides(fill = guide_legend(title = "Patterns")) +
        theme_minimal() + theme(
          legend.position = "bottom"
        )
    } else {
      
      g <- ggplot() +
        geom_sf(data=df, aes(fill=patterns2), color="NA") +
        geom_sf(data = df_states, fill = "NA") +
        scale_fill_manual(values = colors) +
        guides(fill = guide_legend(title = "Patterns")) +
        theme_minimal() + ggtitle(paste(mc_approach, " ", predictor, " vs ", response)) + theme(
          legend.position = "bottom"
        )
      
    }
    
    return(g)
    
     ggsave(path = output_maps, paste0(predictor,response,".png"), height = 4, width = 6, scale = 1.5)
    # print(g)

}
  
  
# CODE BORROWED FROM https://gist.github.com/rafapereirabr/5348193abf779625f5e8c5090776a228
#-- Bivariate Moran's I --#
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  
  xp <- scale(x)[, 1]
  yp <- scale(y)[, 1]
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  
  list(global = global, local  = as.numeric(local))
}

#-- Permutations for the Bivariate Moran's I --#
simula_moran <- function(x, y = NULL, W, nsims = 999){
  
  if(is.null(y)) y = x
  
  n   = nrow(W)
  IDs = 1:n
  
  xp <- scale(x)[, 1]
  W[which(is.na(W))] <- 0
  
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  
  list(global_sims = global_sims,
       local_sims  = local_sims)
}

predictors <- c("ICUBEDS", "OBESITY", "DIABETS", "BLACK", "HISPANC", "WHITE", "UNINSRD", "VACCNTN")
responses <- c("CASS100", "DEATH10")


ggp <- variables_func(responses, predictors, ds_contiguous, states, W, "none") #--- Per the figure foot notes, Saffary looked at responses by predictor as opposed to predictor by response
ggp_PCP_filter <- variables_func(responses, ("PCP"), ds_pcp, states, W2, "none")
ggp_PCP_impute_zero <- variables_func(responses, ("PCP_mpt_z"), ds_contiguous, states, W, "none")
ggp_PCP_impute_mean <- variables_func(responses, ("PCP_mpt_m"), ds_contiguous, states, W, "none")

##-------------------------------------------------------##
#- STEP 5) Perform Adjusted Local Bivariate Moran's I    -#
##-------------------------------------------------------##

##-----------------------------------##
##-------- NOTES ON APPROACH --------##
##-----------------------------------##

# Psuedo P-values were calculated for each county using the same formula as the global Moran's I p-value: 
# Proportion of simulated values (in absolute terms) that are higher than the calculated Moran's I value. 
# I am pretty confident that my method for calculating the pseudo p-values for the supp. figures is correct 
# because I am using the same method that reproduced the global Moran's I p-values in table 2. 
# Although I am calculating pseudo p-values for local rather than global, the method should work for both 
# because I am still comparing calculated Moran's I values against the 999 simulated values.

#-------------------------#
#--BONFERRONI CORRECTION--#
#-------------------------#
# I tried two methods to perform the bonferroni correction. 
# First, I divided the alpha value by the number of comparisons (3105), 
# keeping all other code the same as when we did non-corrected maps (which matched the paper well).  Next, 
# I used the "p.adjust" function from the stats package to perform the bonferroni correction. 
# With this method, I calculated psuedo p-values (See note above). This package produced very similar 
# (if not identical?) to the manually dividing the alpha method. This provides evidence that the psuedo 
# p-value calculation is correct.

#------------------#
#--FDR CORRECTION--#
#------------------#
# I tried two different packages (out of many possible http://www.strimmerlab.org/notes/fdr.html) to calculate 
# the false discover rate. First, I used the "p.adjust" function from the stats package 
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust. The documentation says their 
# fdr algorithm is patterned after this paper: Benjamini, Y., and Hochberg, Y. (1995). 
# Controlling the false discovery rate: a practical and powerful approach to multiple testing. 
# Journal of the Royal Statistical Society Series B, 57, 289--300. http://www.jstor.org/stable/2346101.

# Next, I used the "fdrtool" function from the fdrtool package. The algorithm implemented in the function 
# is described in detail in the documentation on page 5: https://cran.r-project.org/web/packages/fdrtool/fdrtool.pdf 
# on page 5


ggp_supp_fdr <- variables_func(responses, predictors, ds_contiguous, states, W, "fdr")
ggp_supp_bon <- variables_func(responses, predictors, ds_contiguous, states, W, "bonferroni")

ggp_supp_fdr_PCP <- variables_func(responses, ("PCP"), ds_pcp, states, W2, "fdr")
ggp_supp_bon_PCP <- variables_func(responses, ("PCP"), ds_pcp, states, W2, "bonferroni")


# placeholder ggplots
void1 <- ggplot() + theme_void() 
void2 <- ggplot() + theme_void()
void3 <- ggplot() + theme_void()

#--------------------------------#
#---- Figure PCP_impute_zero ----#
#--------------------------------#

ggp_PCP_impute_zero[[1]]
ggsave(path = output_figures, "Cases_PCP_impute_zero.png", height = 4, width = 6, scale = 2.5)

ggp_PCP_impute_zero[[2]]
ggsave(path = output_figures, "Deaths_PCP_impute_zero.png", height = 4, width = 6, scale = 2.5)

#--------------------------------#
#---- Figure PCP_impute_mean ----#
#--------------------------------#

ggp_PCP_impute_mean[[1]]
ggsave(path = output_figures, "Cases_PCP_impute_mean.png", height = 4, width = 6, scale = 2.5)
ggp_PCP_impute_mean[[2]]
ggsave(path = output_figures, "Deaths_PCP_impute_mean.png", height = 4, width = 6, scale = 2.5)

#------------------#
#---- Figure 1 ----#
#------------------#

ggarrange(void1, void2, cases, deaths,
          labels = c("A", "B", "C", "D"), heights = c(0.5,5),
          ncol = 2, nrow = 2)

ggsave(path = output_figures, "fig1_cd.png", height = 4, width = 6, scale = 2.5)

#-------------------------------#
#-- Figure 2 (PCP filter NA) ---#
#-------------------------------#

ggarrange(void1, void2,  ggp[[1]], ggp_PCP_filter[[1]], ggp[[9]], ggp_PCP_filter[[2]],
          labels = c("A", "D", "B", "E", "C", "F"), heights = c(1,5,5),
          ncol = 2, nrow = 3)

ggsave(path = output_figures, "fig2_PCP_filtered_bcef.png", height = 4, width = 6, scale = 2.5)

#---------------------------------#
#-- Figure 2 (PCP impute zero) ---#
#---------------------------------#

ggarrange(void1, void2,  ggp[[1]], ggp_PCP_impute_zero[[1]], ggp[[9]], ggp_PCP_impute_zero[[2]],
          labels = c("A", "D", "B", "E", "C", "F"), heights = c(1,5,5),
          ncol = 2, nrow = 3)

ggsave(path = output_figures, "fig2_PCP_impute_zero_bcef.png", height = 4, width = 6, scale = 2.5)

#---------------------------------#
#-- Figure 2 (PCP impute mean) ---#
#---------------------------------#

ggarrange(void1, void2,  ggp[[1]], ggp_PCP_impute_mean[[1]], ggp[[9]], ggp_PCP_impute_mean[[2]],
          labels = c("A", "D", "B", "E", "C", "F"), heights = c(1,5,5),
          ncol = 2, nrow = 3)

ggsave(path = output_figures, "fig2_PCP_impute_mean_bcef.png", height = 4, width = 6, scale = 2.5)


#------------------#
#---- Figure 3 ----#
#------------------#

ggarrange(void1, void2, ggp[[3]], ggp[[11]],
          labels = c("A", "B", "C", "D"), heights = c(0.5,5),
          ncol = 2, nrow = 2)

ggsave(path = output_figures, "fig3_diabetes_cd.png", height = 4, width = 6, scale = 2.5)


# ------------------#
# ---- Figure 5 ----#
# ------------------#
ggarrange(void1, void2, void3, ggp[[4]], ggp[[5]], ggp[[6]], ggp[[12]], ggp[[13]], ggp[[14]], font.label = list(size = 23),
          labels = c("A", "D", "G", "B", "E", "H", "C", "F", "I"), heights = c(1,5,5),
          ncol = 3, nrow = 3)

ggsave(path = output_figures, "fig5_demographic_bcefhi.png", height = 4, width = 6, scale = 3.8)
