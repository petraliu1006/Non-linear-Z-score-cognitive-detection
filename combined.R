################################################################################
# data loading
################################################################################
rm(list = ls())
setwd('/Users/petra/Desktop/study/UCSF/capstone project/Codes for sharing-selected')

library(fda)
library(splines)
library(splines2)
library(dplyr)
library(magrittr)
library(caret)
library(quadprog)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(mgcv) ## generalized additive models
library(scam) ## shape-constrained additive models, monotonic constraints on the smooth terms.
library(censReg)
library(tidyverse)
library(patchwork)
library(scam)
library(MASS)
library(readxl)
library(pdftools)

source('functions.R')
'%!in%' <- function(x,y)!('%in%'(x,y))
cen = function(x,left=-Inf,right=Inf){
  x[x>=right] = right
  x[x<=left] = left
  return(x)
}
select = dplyr::select

datfull <- read.csv("updated NACC dataset 9.17.21.csv")
### Filtering
datfull %<>% filter(NACCUDSD==1) %<>% # Normal cognition defined by NACC derived variable
  filter(CDRGLOB==0) %<>% # No impairment defined by clinical Dementia Rating (CDR)
  filter(EDUC!=99 & EDUC>=5) %<>% # No missing data in education AND education >= 5 years
  filter(MOCACOMP==1 | CRAFTVRS %!in% c(-4,95:98) | CRAFTURS %!in% c(-4,95:98) | DIGFORCT %!in% c(-4,95:98) | DIGBACCT %!in% c(-4,95:98) | MINTTOTS %!in% c(-4,95:98))     # UDS3 neuropsych test available
datfull = datfull %>% group_by(NACCID) %>% filter(NACCVNUM==min(NACCVNUM)) # Only the first visit where a UDS3 neuropsych test was available
datfull %<>% filter(NPSYLAN==1) %<>% # Language of test administration is English
  filter(DEP==0 | NACCGDS <= 5 | NACCGDS == 88) %<>% # No active depression
  filter(BIPOLDX==0 & SCHIZOP==0 & DELIR==0 & PTSDDX==0 & ALCABUSE!=1) %<>% # No active bipolar disorder, schizophrenia, delirium, PTSD, or current alcohol abuse
  filter(OTHPSY==0 | OTHPSYX %in% c("adj d/o i mixed depression + anxiety","adjustment d/o","adjustment disorder",
                                    "adjustment d/o depression","Agitation; Irritability","ANXIETY","Bereavement","bereavement",
                                    "depression d/o NOS 2o caregiver stress","dysthymia","DYSTHMIA","Dysthymine","Irritability",
                                    "MDD in full remission","mdd in full remission")) %<>% # No active other mood disorders
  filter(NACCFADM==0 & NACCFFTD==0 & NACCOMX %!in% c("?","C9orf72","GSS","GSS 117","GSS 198","GSS F198","GSS F198-IK",
                                                     "GSS Q217P","GSS-102","GSS-217","HDLS","Iron mutation","Jack 2, for blood cancer",
                                                     "Lewy Body", "LRRK","PD LAARK","SCA-3","unknown","Unknown","X-linked recessive ataxia")) # No dominantly inherited AD mutation in the family, family history of FTLD mutation, or potential mutation other than an AD or FTLD mutation

### Data manipulation
datfull %<>% mutate(SEX = 2-SEX) # Sex = 1: male; 0: female
datfull %<>% filter(RACE %in% c(1,2,4,5) & HISPANIC != 9)  #filter white, black, native, asian
datfull %<>% mutate(RACEETHwhitenonhispanic = ifelse(RACE==1 & HISPANIC==0, 1,0),
                    RACEETHblacknonhispanic = ifelse(RACE==2 & HISPANIC==0, 1,0),
                    RACEETHasiannonhispanic = ifelse(RACE %in% c(4,5) & HISPANIC==0, 1,0),
                    # RACEETHindiannonhispanic = ifelse(RACE==3 & HISPANIC==0, 1,0),
                    RACEETHhispanic = ifelse(HISPANIC==1, 1,0))
datfull %<>% mutate(EDUCCAT = cut(EDUC,
                                  breaks = c(0,12.5,16.5,30),
                                  labels=c("possibily high school or lower","possibily some college","possibily higher than college"),
                                  include.lowest=TRUE)) %<>%
  mutate(EDUCCAThighschoollower = ifelse(EDUCCAT=="possibily high school or lower",1,0),
         EDUCCATsomecollege = ifelse(EDUCCAT=="possibily some college",1,0),
         EDUCCAThigherthancollege = ifelse(EDUCCAT=="possibily higher than college",1,0))

# Age truncated at 40 and 95
datfull %<>% mutate(NACCAGE = ifelse(NACCAGE<=40,40,NACCAGE)) %<>% 
  mutate(NACCAGE = ifelse(NACCAGE>=95,95,NACCAGE))

### Creating ANIVEG as sum of ANIMALS and VEG
datfull$ANIVEG <- rep(200,nrow(datfull)) ## 200 is above value where would be cut out
for(i in 1:nrow(datfull)) 
  if(is.na(datfull$ANIMALS[i])==FALSE & is.na(datfull$VEG[i]) == FALSE & datfull$ANIMALS[i] >=(-0.5) & datfull$ANIMALS[i] < 77.5 & datfull$VEG[i] >=(-0.5) & datfull$VEG[i] < 77.5) 
    datfull$ANIVEG[i] <- datfull$ANIMALS[i] + datfull$VEG[i]
### sum of Number of correct F&L-words generated in 1 minute 
datfull$UDSVERTN2 <- rep(200,nrow(datfull)) ## 200 is above value where would be cut out
for(i in 1:nrow(datfull)) 
  if(is.na(datfull$UDSVERFC[i])==FALSE & is.na(datfull$UDSVERLC[i]) == FALSE & datfull$UDSVERFC[i] >=(-0.5) & datfull$UDSVERFC[i] < 40.5 & datfull$UDSVERLC[i] >=(-0.5) & datfull$UDSVERLC[i] < 40.5) 
    datfull$UDSVERTN2[i] <- datfull$UDSVERFC[i] + datfull$UDSVERLC[i]

### Extract subset data of interest 
varset <- c("TRAILA","TRAILB","UDSVERFC","UDSVERLC","UDSVERTN","ANIMALS","VEG","ANIVEG","MINTTOTS","DIGFORSL","DIGBACLS","CRAFTVRS","CRAFTDVR","UDSBENTC","UDSBENTD","MOCATOTS","DIGFORCT","DIGBACCT")
varnames <- c("Trail Making Test A","Trail Making Test B","Letter fluency F","Letter fluency L","Letter fluency F + L total","Category fluency - animals","Category fluency - vegetables","Category fluency - total","Multilingual naming test (MINT) total","Number Span longest digit forward","Number Span longest digit backward","Craft memory - immediate","Craft memory - delay","Benson figure - copy","Benson figure - recall","Montreal Cognitive Assessment (MoCA) total","Number Span forward total correct trials","Number Span backward total correct trials")
varset_transformed <- paste0(varset,'_t')


var = c("NACCID","NACCAGE","SEX","EDUC","RACE",varset)
### Coerce as numeric variable
subdata = datfull %>% dplyr::select(var,starts_with("RACEETH"),starts_with("EDUCCAT")) %>% mutate_at(var[-(1:3)], ~as.numeric(as.character(., na.rm = TRUE))) %>% as.data.frame()

################################################################################

p.list = list()
for(indx in 1:length(varset)) {
  var_name = varset[indx]
  print(paste("Start ",var_name,", n = ",indx,".",sep=""))
  dir.create(paste0("Results0209/",var_name), recursive = TRUE, showWarnings = FALSE)
  source(paste0("Parsets/",var_name,"pars.R"))
  
  #filter
  if(var_name == "TRAILB") {
    dat <- subdata %>% 
      filter(get(var_name) >= mnval & get(var_name) <= mxval, !(get(var_name) == 300))
    datj <- subdata %>% 
      filter(get(var_name) >= mnval & get(var_name) <= mxval)
  } else {
    dat <- subdata %>% 
      filter(get(var_name) >= mnval & get(var_name) <= mxval)
  }
  
  # Box-Cox transformation with error handling
  if(min(dat[[var_name]], na.rm=TRUE) <= 0) {
    constant <- abs(min(dat[[var_name]], na.rm=TRUE)) + 1
    dat[[var_name]] <- dat[[var_name]] + constant
    shift_constant <- constant
  } else {
    shift_constant <- 0
  }
  
  # Perform Box-Cox transformation
  bc <- boxcox(lm(dat[[var_name]] ~ 1), lambda = seq(-2, 2, 0.1))
  best_lambda <- bc$x[which.max(bc$y)]
  transformed_values <- if (best_lambda == 0) {
    log(dat[[var_name]])
  } else {
    (dat[[var_name]]^best_lambda - 1) / best_lambda
  }
  
  # Add transformed values to dataset
  dat[[paste0(var_name, "_transformed")]] <- transformed_values
  print(paste("Transformation completed for", var_name))
  print(paste("Lambda:", best_lambda))
  print(paste("Shift constant:", shift_constant))
  
  if(incage == TRUE) agedrn <- "mpi"
  if(incage == FALSE) agedrn <- "mpd"
  
  sd <- sd(dat[[var_name]], na.rm=T)
  sd_t <- sd(dat[[paste0(var_name, "_transformed")]], na.rm=T)
  dat <- dat %>% mutate(RACE_model = case_when(RACEETHwhitenonhispanic == 1 ~ 0,RACEETHblacknonhispanic == 1 ~ 1,TRUE ~ NA_real_))
  
  # Fit models
  john <- try(eval(substitute(scam(foo ~ s(NACCAGE,bs=agedrn) + EDUC + SEX + NACCAGE:RACE_model + RACE_model,  data=dat), list(foo=as.name(paste0(var_name, "_transformed"))))))
  john_orginal <- try(eval(substitute(scam(foo ~ s(NACCAGE,bs=agedrn) + EDUC + SEX + NACCAGE:RACE_model + RACE_model,  data=dat), list(foo=as.name(var_name)) )))
  age.range <- range(dat$NACCAGE)
  
  # Create prediction data frames
  scam_t_df <- john$model[1:(age.range[2]-age.range[1]+1),]
  scam_t_df_white <- scam_t_df %>% mutate(   NACCAGE = (age.range[1]:age.range[2]), EDUC = 17,RACE_model = 0, SEX = 1)
  scam_t_df_black <- scam_t_df %>% mutate( NACCAGE = (age.range[1]:age.range[2]),EDUC = 17,  RACE_model = 1,  SEX = 1 )
  scam_org_df = john_orginal$model[1:(age.range[2]-age.range[1]+1),] 
  scam_org_df_white = scam_org_df %<>% mutate(NACCAGE = (age.range[1]:age.range[2]), EDUC = 17, RACE_model = 0, SEX = 1)
  scam_org_df_black = scam_org_df %<>% mutate(NACCAGE = (age.range[1]:age.range[2]), EDUC = 17, RACE_model = 1, SEX = 1)
  
  # Generate predictions
  pred_john_t_white <- data.frame( pred = as.numeric(predict.scam(john, newdata=scam_t_df_white, type="response")), age = (age.range[1]:age.range[2]))
  pred_john_t_black <- data.frame( pred = as.numeric(predict.scam(john, newdata=scam_t_df_black, type="response")), age = (age.range[1]:age.range[2]))
  pred_john_org_white = data.frame(pred=as.numeric(predict.scam(john_orginal,newdata=scam_org_df_white,type="response")),age=(age.range[1]:age.range[2]))
  pred_john_org_black = data.frame(pred=as.numeric(predict.scam(john_orginal,newdata=scam_org_df_black,type="response")),age=(age.range[1]:age.range[2]))
  
  # Calculate residuals and model frames
  dat_john_full <- model.frame(john)
  dat_john_t_white <- dat_john_full %>% filter(RACE_model == 0)
  dat_john_t_black <- dat_john_full %>% filter(RACE_model == 1)
  dat_john_t_white$residuals <- john$y[dat_john_full$RACE_model == 0] - john$fitted.values[dat_john_full$RACE_model == 0]
  dat_john_t_white$yVec <- john$y[dat_john_full$RACE_model == 0]
  dat_john_t_white$XB <- john$fitted.values[dat_john_full$RACE_model == 0]
  dat_john_t_black$residuals <- john$y[dat_john_full$RACE_model == 1] - john$fitted.values[dat_john_full$RACE_model == 1]
  dat_john_t_black$yVec <- john$y[dat_john_full$RACE_model == 1]
  dat_john_t_black$XB <- john$fitted.values[dat_john_full$RACE_model == 1]
  
  dat_john_org_full = model.frame(john_orginal)
  dat_john_org_white = dat_john_org_full %>% filter(RACE_model == 0)
  dat_john_org_black = dat_john_org_full %>% filter(RACE_model == 1)
  dat_john_org_white$residuals <- john_orginal$y[dat_john_org_full$RACE_model == 0] - john_orginal$fitted.values[dat_john_org_full$RACE_model == 0]
  dat_john_org_white$yVec <- john_orginal$y[dat_john_org_full$RACE_model == 0]
  dat_john_org_white$XB <- john_orginal$fitted.values[dat_john_org_full$RACE_model == 0]
  dat_john_org_black$residuals <- john_orginal$y[dat_john_org_full$RACE_model == 1] - john_orginal$fitted.values[dat_john_org_full$RACE_model == 1]
  dat_john_org_black$yVec <- john_orginal$y[dat_john_org_full$RACE_model == 1]
  dat_john_org_black$XB <- john_orginal$fitted.values[dat_john_org_full$RACE_model == 1]
  
  # Set up for SD calculations
  age <- (age.range[1]:age.range[2])
  wdth <- 5
  
  sdResidout1t <- Nout1t <- vector(length=length(age))
  sdResidout3t <- Nout3t <- vector(length=length(age))
  sdResidout1org <- Nout1org <- vector(length=length(age))
  sdResidout3org <- Nout3org <- vector(length=length(age))
  
  # Calculate SDs using moving window
  for(i in age) {
    datsub1t <- dat_john_t_white %>%
      filter(NACCAGE >= (i-wdth) & NACCAGE <= (i+wdth))
    sdResidout1t[i-age[1]+1] <- sd(datsub1t$residuals)
    Nout1t[i - age[1] + 1] <- nrow(datsub1t)
    
    datsub3t <- dat_john_t_black %>%
      filter(NACCAGE >= (i-wdth) & NACCAGE <= (i+wdth))
    sdResidout3t[i-age[1]+1] <- sd(datsub3t$residuals)
    Nout3t[i - age[1] + 1] <- nrow(datsub3t)
  }
  for(i in age){
    datsub1org <- dat_john_org_white %>% filter(NACCAGE >= (i-wdth) & NACCAGE <= (i+wdth)) 
    sdResidout1org[i-age[1]+1] <- sd(datsub1org$residuals) 
    Nout1org[i - age[1] + 1] <- nrow(datsub1org)
    
    datsub3org <- dat_john_org_black %>% filter(NACCAGE >= (i-wdth) & NACCAGE <= (i+wdth))
    sdResidout3org[i-age[1]+1] <- sd(datsub3org$residuals)
    Nout3org[i - age[1] + 1] <- nrow(datsub3org)
  }
  
  print(datsub1t)
  print(datsub3t)
  
  # Create data frames for SD modeling
  resultDFwithoutliers1t <- data.frame(age = age,  sdResidout = sdResidout1t, Nout = Nout1t)
  resultDFwithoutliers3t <- data.frame( age = age,  sdResidout = sdResidout3t, Nout = Nout3t)
  resultDFwithoutliers1org <- data.frame(age,sdResidout=sdResidout1org,Nout=Nout1org)
  resultDFwithoutliers3org <- data.frame(age,sdResidout=sdResidout3org,Nout=Nout3org)
  
  # Fit SD curves
  res_mod1t <- scam(sdResidout ~ s(age, bs = "mpi"),  weights = Nout,data = resultDFwithoutliers1t)
  res_mod3t <- scam(sdResidout ~ s(age, bs = "mpi"), weights = Nout,  data = resultDFwithoutliers3t)
  res_mod1org = scam(sdResidout ~ s(age, bs = "mpi"), weights = Nout, data = resultDFwithoutliers1org) #weight:ensures groups with more data have a greater influence on the model fit
  res_mod3org = scam(sdResidout ~ s(age, bs = "mpi"), weights = Nout, data = resultDFwithoutliers3org)
  # Generate SD predictions
  sd_john_t_white <- data.frame( sdpred = predict.scam(res_mod1t, newdata = resultDFwithoutliers1t, type = "response"),age = age, sd = resultDFwithoutliers1t$sdResidout)
  sd_john_t_black <- data.frame(sdpred = predict.scam(res_mod3t, newdata = resultDFwithoutliers3t, type = "response"), age = age, sd = resultDFwithoutliers3t$sdResidout)
  sd_john_org_white = data.frame(sdpred = predict.scam(res_mod1org, newdata = resultDFwithoutliers1org, type = "response"),age = age,sd = resultDFwithoutliers1org$sdResidout)
  sd_john_org_black = data.frame(sdpred = predict.scam(res_mod3org, newdata = resultDFwithoutliers3org, type = "response"),age = age,sd = resultDFwithoutliers3org$sdResidout)
  
  print(range(sd_john_t_white$sdpred))
  print(range(sd_john_t_black$sdpred))
  print(mean(sd_john_t_white$sdpred))
  print(median(sd_john_t_white$sdpred))
  print(mean(sd_john_t_black$sdpred))
  print(median(sd_john_t_black$sdpred))
  
  print(1 - ((1 - summary(john)$r.sq) * (length(john$y) - 1) / (length(john$y) - length(john$coefficients) - 1)))
  print(1 - ((1 - summary(john_orginal)$r.sq) * (length(john$y) - 1) / (length(john$y) - length(john$coefficients) - 1)))
  
  # Inverse Box-Cox transformation
  inverse_boxcox <- function(x, lambda, constant = 0) {
    if (lambda == 0) {
      result <- exp(x)
    } else {
      result <- (x * lambda + 1)^(1/lambda)
    }
    return(result - constant)
  }

    # Apply inverse transformation to predictions
  pred_john_t_white$pred_original <- inverse_boxcox(pred_john_t_white$pred, best_lambda)
  pred_john_t_black$pred_original <- inverse_boxcox(pred_john_t_black$pred, best_lambda)
  
  # Create modified prediction dataframes for plotting
  pred_john_t_white_mod <- pred_john_t_white %>% select(pred, age) %>% mutate(pred_original = inverse_boxcox(pred, best_lambda)) %>% select(-pred) %>% rename(pred = pred_original) %>% mutate(Population = "White", Model = "Inverse-transform Mean/SD")
  pred_john_t_black_mod <- pred_john_t_black %>%  select(pred, age) %>%  mutate(pred_original = inverse_boxcox(pred, best_lambda)) %>%  select(-pred) %>%  rename(pred = pred_original) %>% mutate(Population = "Black", Model = "Inverse-transform Mean/SD")
  pred_john_org_white_mod <- pred_john_org_white %>% select(pred, age) %>% mutate(Population = "White", Model = "Original Mean/SD")
  pred_john_org_black_mod <- pred_john_org_black %>% select(pred, age) %>% mutate(Population = "Black", Model = "Original Mean/SD")
  
  pred_df1 <- rbind.data.frame(
    pred_john_t_white %>% bind_cols(Population="White"),
    pred_john_t_black %>% bind_cols(Population="Black")) %>% 
    mutate(Population = factor(Population, levels = c("White", "Black")))
  
  p1 <- dat %>%
    filter(RACEETHwhitenonhispanic == 1 | RACEETHblacknonhispanic == 1) %>%
    rename(y = !!sym(paste0(var_name, "_transformed"))) %>%
    mutate(Sample = ifelse(RACEETHwhitenonhispanic == 1, "White Mean/SD", "Black Mean/SD")) %>%
    ggplot(aes(x = NACCAGE, y = y, colour = Sample)) +
    geom_point(position = position_jitter(h = 0.3, w = 0.3), alpha = 0.3,size = 0.5, show.legend = TRUE) +
    geom_line(data = pred_df1, aes(x = age, y = pred, colour = Population), linewidth = 1.2, show.legend = TRUE) +
    scale_colour_manual(
      name = NULL,
      values = c(
        "White Mean/SD" = "#9D7FCC",     
        "Black Mean/SD" = "grey10",     
        "White" = "#5E3C99", 
        "Black" = "grey10"),
      guide = guide_legend(override.aes = list(shape = c(19, 19, NA, NA), linetype = c("blank", "blank", "solid", "solid"),size = c(2, 2, 1, 1),alpha = c(1, 1, 1, 1)))) + 
    xlab("Age") + 
    ylab(var_name) + 
    ggtitle("Fitted Mean (Box-transformed)") + 
    theme_bw()  
  
  p2 <- ggplot() +
    geom_point(data = sd_john_t_white, aes(x = age, y = sd), colour = "#5E3C99", alpha = 0.2,size = 0.7) +
    geom_line(data = sd_john_t_white, aes(x = age, y = sdpred), colour = "#5E3C99", linewidth = 1.2) +
    geom_point(data = sd_john_t_black, aes(x = age, y = sd), colour = "black", alpha = 0.2,size = 0.7) +
    geom_line(data = sd_john_t_black, aes(x = age, y = sdpred), colour = "black", linewidth = 1.2) +
    xlab("Age") + ylab(var_name) + ggtitle("Fitted SD  (Box-transformed)") + theme_bw() + theme(legend.position = "right")
  
  pred_df <- rbind( pred_john_t_white_mod, pred_john_t_black_mod,pred_john_org_white_mod,pred_john_org_black_mod) %>%
    mutate(  Model = factor(Model, levels = c("Original Mean/SD", "Inverse-transform Mean/SD")),  
             Population = factor(Population, levels = c("White", "Black")))

  p3 <- dat %>%
    filter(RACEETHwhitenonhispanic == 1 | RACEETHblacknonhispanic == 1) %>%
    rename(y = !!sym(var_name)) %>%
    mutate(Sample = ifelse(RACEETHwhitenonhispanic == 1, "White Mean/SD", "Black Mean/SD")) %>%
    ggplot(aes(x = NACCAGE, y = y, colour = Sample)) +
    geom_point(position = position_jitter(h = 0.3, w = 0.3), alpha = 0.2, size = 0.5,show.legend = TRUE) +
    geom_line(data = pred_df, aes(x = age, y = pred, colour = Population,linetype = Model), linewidth = 1.2,show.legend = TRUE) +
    scale_colour_manual(
      name = NULL,
      values = c(
        "White Mean/SD" = "#9D7FCC",
        "Black Mean/SD" = "grey10",
        "White" = "#5E3C99",
        "Black" = "grey10")) +
    scale_linetype_manual(name = NULL,values = c("Original Mean/SD" = "longdash","Inverse-transform Mean/SD" = "dotdash")) +
    xlab("Age") +
    ylab(var_name) +
    ggtitle("Fitted Mean (Original vs Inverse-transform)") +
    theme_bw() +
    guides(colour = guide_legend(order = 1,override.aes = list(shape = c(19, 19, NA, NA),linetype = c("blank", "blank", "solid", "solid"),size = c(2, 2, 1, 1),alpha = c(1, 1, 1, 1))),
      linetype = guide_legend(order = 2, override.aes = list(shape = NA)))
  
  sd_john_t_white <- sd_john_t_white %>% mutate(sdpred_original = sdpred * abs(inverse_boxcox(pred_john_t_white$pred, best_lambda))^(1-best_lambda), sd_original = sd * abs(inverse_boxcox(pred_john_t_white$pred, best_lambda))^(1-best_lambda))
  sd_john_t_black <- sd_john_t_black %>% mutate(sdpred_original = sdpred * abs(inverse_boxcox(pred_john_t_black$pred, best_lambda))^(1-best_lambda), sd_original = sd * abs(inverse_boxcox(pred_john_t_black$pred, best_lambda))^(1-best_lambda))
  
  p4 <- ggplot() +
    geom_point(data = sd_john_t_white, aes(x = age, y = sd_original), colour = "#5E3C99", alpha = 0.2,size = 0.7) +
    geom_line(data = sd_john_t_white, aes(x = age, y = sdpred_original), colour = "#5E3C99", linewidth = 1.2, linetype = "dotdash") +
    geom_point(data = sd_john_t_black, aes(x = age, y = sd_original), colour = "black", alpha = 0.2,size = 0.7) +
    geom_line(data = sd_john_t_black, aes(x = age, y = sdpred_original), colour = "black", linewidth = 1.2, linetype = "dotdash") +
    geom_point(data = sd_john_org_white, aes(x = age, y = sd), colour = "#5E3C99", alpha = 0.2,size = 0.7) +
    geom_line(data = sd_john_org_white, aes(x = age, y = sdpred), colour = "#5E3C99", linewidth = 1, linetype = "longdash") +
    geom_point(data = sd_john_org_black, aes(x = age, y = sd), colour = "black", alpha = 0.2,size = 0.7) +
    geom_line(data = sd_john_org_black, aes(x = age, y = sdpred), colour = "black", linewidth = 1, linewidth = 1, linetype = "longdash") +
    xlab("Age") +
    ylab(var_name) +
    ggtitle("Fitted SD (Original vs Inverse Transformed)") +
    theme_bw() +
    theme(legend.position = "right")
  
  # Modify the combined plot to include all 4 plots
  combined_plot <- p1+p2 + p3+p4+
    plot_layout(ncol = 2) +
    plot_annotation(
      title = paste("Results for", var_name),
      theme = theme(plot.title = element_text(hjust = 0.5)))
  
  # Save combined plot
  ggsave(
    filename = paste0("Results0209/", var_name, "/combined_plots.pdf"),
    plot = combined_plot,
    width = 16,
    height = 15  )
  
  # Store plot in the list for later use
  p.list[[var_name]] <- combined_plot
  
  
  ################################################################################
    # Generate lookup table
  age = 40:95
  educ = 10:20
  sex = 0:1
  race = c("White", "Black")
  
  datout <- expand.grid( NACCAGE = age, EDUC = educ, SEX = sex, RACE = race)
  
  # Add model information
  datout <- datout %>% mutate(model = "transformed", RACE_model = case_when(RACE == "White" ~ 0, RACE == "Black" ~ 1))
  
  # Prepare data for predictions
  datout_white <- datout %>% filter(RACE == "White") %>% mutate(age = NACCAGE)
  datout_black <- datout %>% filter(RACE == "Black") %>% mutate(age = NACCAGE)
  
  # Generate predictions for lookup table
  datout$mean.adj <- NA
  datout$sd.adj <- NA
  
  # Add predictions
  datout$mean.adj[datout$RACE == "White"] <- predict(john, newdata = datout_white, type = "response")
  datout$mean.adj[datout$RACE == "Black"] <- predict(john, newdata = datout_black, type = "response")
  datout$sd.adj[datout$RACE == "White"] <- predict(res_mod1t, newdata = datout_white, type = "response")
  datout$sd.adj[datout$RACE == "Black"] <- predict(res_mod3t, newdata = datout_black, type = "response")
  
  datout$mean.adj.ori[datout$RACE == "White"] <- predict(john_orginal, newdata = datout_white, type = "response")
  datout$mean.adj.ori[datout$RACE == "Black"] <- predict(john_orginal, newdata = datout_black, type = "response")
  datout$sd.adj.ori[datout$RACE == "White"] <- predict(res_mod1org, newdata = datout_white, type = "response")
  datout$sd.adj.ori[datout$RACE == "Black"] <- predict(res_mod3org, newdata = datout_black, type = "response")
  
  # Save lookup table
  write.csv(datout, paste0("Results_t/", var_name, "/Zscore_lookup_transformed.csv"), row.names = FALSE)
}

# After the loop, combine all plots into one PDF if needed
if(length(p.list) > 0) {
  pdf_path <- "Results0209/all_variables_combined.pdf"
  pdf(pdf_path, width = 13, height = 8)
  for(p in p.list) {
    print(p)
  }
  dev.off()
  print(paste("Combined PDF saved to:", pdf_path))
}


################################################################################
# Table for Julie TRAIL A from box-transformed model
################################################################################
library(readxl)
Z.TMT.A = read.csv("Results/TRAILA/Zscore_lookup_transformed.csv")
tab = read_xlsx("Trails Norms Comparisons_Diff Dem Adj_jw_JAF.xlsx", sheet = 1)
tab = tab %>% 
  slice(2:nrow(tab)) %>%
  select("ALZ", "Age", "Ed", "Sex", "TMTARaw", "TMTBRaw")

# Clean and transform the data
tab = tab %>%  mutate(   NACCAGE = case_when(     Age <= 40 ~ 40,     Age >= 95 ~ 95,     TRUE ~ Age   ),
                         EDUC = case_when(   Ed <= 10 ~ 10,   Ed >= 20 ~ 20,   TRUE ~ Ed ), SEX = ifelse(Sex == "M", 1, 0)  ) %>%
  mutate(  NACCAGE = as.integer(NACCAGE),  EDUC = as.integer(EDUC),  SEX = as.integer(SEX))

# Apply Box-Cox transformation
tab$TMTARawtransformed <-  (tab$TMTARaw^(-0.343434343434343) - 1) / (-0.343434343434343)

# Merge with Z-score lookup tables and calculate Z-scores
tab = tab %>%  left_join(   Z.TMT.A %>% filter(RACE == "White") %>% rename(TMT_A_mean_white = mean.adj, TMT_A_sd_white = sd.adj) %>%  select(NACCAGE, EDUC, SEX, TMT_A_mean_white, TMT_A_sd_white), by = c("NACCAGE", "EDUC", "SEX") ) %>%
  left_join(  Z.TMT.A %>%  filter(RACE == "Black") %>%   rename(TMT_A_mean_black = mean.adj,   TMT_A_sd_black = sd.adj) %>%  select(NACCAGE, EDUC, SEX, TMT_A_mean_black, TMT_A_sd_black), by = c("NACCAGE", "EDUC", "SEX") ) %>%
  mutate(  Z_TMT_A_white = (TMT_A_mean_white - TMTARawtransformed) / TMT_A_sd_white,Z_TMT_A_black = (TMT_A_mean_black - TMTARawtransformed) / TMT_A_sd_black)

# Write final output
write.csv(tab, "TMTA_Zscores_Output_new.csv", row.names = FALSE)

################################################################################
# Table for Julie TRAIL B from box-transformed model
################################################################################
library(readxl)
Z.TMT.B = read.csv("Results/TRAILB/Zscore_lookup_transformed.csv")
tab = read_xlsx("Trails Norms Comparisons_Diff Dem Adj_jw_JAF.xlsx", sheet = 1)
tab = tab %>% 
  slice(2:nrow(tab)) %>%
  select("ALZ", "Age", "Ed", "Sex", "TMTARaw", "TMTBRaw")

# Clean and transform the data
tab = tab %>%  mutate(   NACCAGE = case_when(     Age <= 40 ~ 40,     Age >= 95 ~ 95,     TRUE ~ Age   ),
                         EDUC = case_when(   Ed <= 10 ~ 10,   Ed >= 20 ~ 20,   TRUE ~ Ed ), SEX = ifelse(Sex == "M", 1, 0)  ) %>%
  mutate(  NACCAGE = as.integer(NACCAGE),  EDUC = as.integer(EDUC),  SEX = as.integer(SEX))

# Apply Box-Cox transformation
tab$TMTBRawtransformed <- (tab$TMTBRaw^(-0.303030303030303) - 1) / (-0.303030303030303)

# Merge with Z-score lookup tables and calculate Z-scores
tab = tab %>%  left_join(   Z.TMT.B %>% filter(RACE == "White") %>% rename(TMT_B_mean_white = mean.adj, TMT_B_sd_white = sd.adj) %>%  select(NACCAGE, EDUC, SEX, TMT_B_mean_white, TMT_B_sd_white), by = c("NACCAGE", "EDUC", "SEX") ) %>%
  left_join(  Z.TMT.B %>%  filter(RACE == "Black") %>%   rename(TMT_B_mean_black = mean.adj,   TMT_B_sd_black = sd.adj) %>%  select(NACCAGE, EDUC, SEX, TMT_B_mean_black, TMT_B_sd_black), by = c("NACCAGE", "EDUC", "SEX") ) %>%
  mutate(  Z_TMT_B_white = (TMT_B_mean_white - TMTBRawtransformed) / TMT_B_sd_white,Z_TMT_B_black = (TMT_B_mean_black - TMTBRawtransformed) / TMT_B_sd_black)

# Write final output
write.csv(tab, "TMTB_Zscores_Output_new.csv", row.names = FALSE)




################################################################################
##lambda
################################################################################

for(indx in 1:length(varset)) {
  var_name = varset[indx]
  print(paste("Start ",var_name,", n = ",indx,".",sep=""))
  dir.create(paste0("Results0209/",var_name), recursive = TRUE, showWarnings = FALSE)
  source(paste0("Parsets/",var_name,"pars.R"))
  
  #filter
  if(var_name == "TRAILB") {
    dat <- subdata %>% 
      filter(get(var_name) >= mnval & get(var_name) <= mxval, !(get(var_name) == 300))
    datj <- subdata %>% 
      filter(get(var_name) >= mnval & get(var_name) <= mxval)
  } else {
    dat <- subdata %>% 
      filter(get(var_name) >= mnval & get(var_name) <= mxval)
  }
  
  # Box-Cox transformation with error handling
  if(min(dat[[var_name]], na.rm=TRUE) <= 0) {
    constant <- abs(min(dat[[var_name]], na.rm=TRUE)) + 1
    dat[[var_name]] <- dat[[var_name]] + constant
    shift_constant <- constant
  } else {
    shift_constant <- 0
  }
  
  # Perform Box-Cox transformation
  bc <- boxcox(lm(dat[[var_name]] ~ 1), lambda = seq(-2, 2, 0.1))
  best_lambda <- bc$x[which.max(bc$y)]
  transformed_values <- if (best_lambda == 0) {
    log(dat[[var_name]])
  } else {
    (dat[[var_name]]^best_lambda - 1) / best_lambda
  }
  
  # Add transformed values to dataset
  dat[[paste0(var_name, "_transformed")]] <- transformed_values
  print(paste("Transformation completed for", var_name))
  print(paste("Lambda:", best_lambda))
  print(paste("Shift constant:", shift_constant))
  
  if(incage == TRUE) agedrn <- "mpi"
  if(incage == FALSE) agedrn <- "mpd"
  
  sd <- sd(dat[[var_name]], na.rm=T)
  sd_t <- sd(dat[[paste0(var_name, "_transformed")]], na.rm=T)
  dat <- dat %>% mutate(RACE_model = case_when(RACEETHwhitenonhispanic == 1 ~ 0,RACEETHblacknonhispanic == 1 ~ 1,TRUE ~ NA_real_))
  
  # Fit models
  john <- try(eval(substitute(scam(foo ~ s(NACCAGE,bs=agedrn) + EDUC + SEX + NACCAGE:RACE_model + RACE_model,  data=dat), list(foo=as.name(paste0(var_name, "_transformed"))))))
  john_orginal <- try(eval(substitute(scam(foo ~ s(NACCAGE,bs=agedrn) + EDUC + SEX + NACCAGE:RACE_model + RACE_model,  data=dat), list(foo=as.name(var_name)) )))
  age.range <- range(dat$NACCAGE)
}
  