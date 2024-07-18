## Script name: R codes for baseline PS instruction.R
## Purpose of script: The codes in the manuscript introcudtion PS methods
## Author: Franco, FAN Min
## Date Created: 2024-07-16
## Copyright (c) FAN Min, 2024
## Email: minfan@connect.hku.hk
##
## Notes:
##  1. FM include all dx, rx, into the ps model
##  2.dx.hepb is all with 0, so i delete
##  3. the factors makes the model faster: 14.3 mins vs 7.7 s when the age is continuous and other covariates are factor
##

# temp function

time_consum <- function(gogo){
    start.time <- Sys.time()
    gogo
    end.time <- Sys.time()
    round(end.time - start.time,2)
}

get_coef <- function(model){
    temp <- as.data.table(cbind(exp(coef(model)),exp(confint.default(model))),keep.rownames = T)
    colnames(temp) <- c("var","est","lwr","upr")
    temp <- temp[var=="COVID"]
    return(temp)
}

# Package loading and environment setting ----------------------------------
options(scipen = 6, digits = 6)
memory.limit(30000000)
library(data.table) # data manipulation
library(tableone) # creating table one and obtaining SMD
library(WeightIt) # example with package

# Data preparation
dt <- readRDS("data/ACESO-renal-cohort.sample-20240128.RDS")
dt$age <- 2024-dt$dob_y
## PS manually estimate ----------------------------------------------------
cova <- c("sex","age",grep("dx|rx",colnames(dt),value = T))
cova <- setdiff(cova,"dx.hepb")
dt[,(cova) := lapply(.SD, as.factor),.SDcols = cova]
dt[,age:=as.numeric(age)]

ps_formula <- as.formula(paste0("COVID~",paste(cova,collapse = "+")))
ps_model <- glm(formula = ps_formula,data = dt,family = "binomial")
#saveRDS(ps_model,file = "data/psmodel_result.RDS")
dt$psvalue <- ps_model$fitted.values

tableone_ps <- print(CreateTableOne(vars = c(cova,"psvalue"),strata = "COVID",data = dt,test = FALSE),smd =T)
as.data.table(tableone_ps,keep.rownames = T)[as.numeric(SMD)>0.1]

# regression adjustment ---------------------------------------------------

## outcome model
outcome_model_regadj <- glm(outcome.MAKE~COVID+psvalue,data = dt,family="binomial")
get_coef(outcome_model_regadj)

# matching ----------------------------------------------------------------


# stratification ----------------------------------------------------------
### obtain the max and min of the ps value group by exposure group; get the max of min, min of the max
temp_range <- dt[,.(max=max(psvalue),min=min(psvalue)),COVID][,.(min=max(min),max=min(max))]
dt_stratification<-dt[psvalue>=temp_range$min&psvalue<=temp_range$max]
### make strata accordiing to the ps
dt_stratification<- dt_stratification[order(psvalue),][,rank_ps:=.I][,strata:=floor(rank_ps*5/.N +1)][]

temp_stratification <-
merge(dt_stratification[COVID == 1, .(ct = .N), by = strata][order(strata, -ct)][!duplicated(strata)][, .(strata, exp = 1, ct_exp = ct)],
dt_stratification[COVID == 0, .(ct = .N), by = strata][order(strata, -ct)][!duplicated(strata)][, .(strata, nonexp = 0, ct_nonexp = ct)],by="strata")[
    !is.na(ct_exp) & !is.na(ct_nonexp) & ct_exp != 0 & ct_nonexp != 0
    ][,`:=`(w_unexp=(ct_exp/sum(ct_exp))/(ct_nonexp/sum(ct_nonexp)),w_exp=1)
      ]
temp_stratification<-rbind(temp_stratification[,.(strata,w=w_exp,COVID=exp)],
                           temp_stratification[,.(strata,w=w_unexp,COVID=nonexp)])[]


dt_stratification <- merge(dt_stratification,temp_stratification,by=c("strata","COVID"))

unloadNamespace("tableone")
library(survey)

trim_weight <- svydesign(ids = ~ patient_pssn, strata = ~ COVID, weights = ~ w,
                         data = dt_stratification)
library(tableone)
tab1 <- svyCreateTableOne(vars = varsp,
                          strata = "dm_sas", data = trim_weight,
                          factorVars = catvarsp)
baseline<-as.data.frame(as.data.table(print(tab1,test = F,noSpaces = T,smd = TRUE),
                                      keep.rownames = T))


# weighting ---------------------------------------------------------------

## IPTW
dt[,ipw:=COVID/psvalue+(1-COVID)/(1-psvalue)]
ipw_pkg <- weightit(ps_formula,
                    data = dt,
                    estimand = "ATE",
                    method = "ps")
dt[,ipw_pkg:=ipw_pkg$weights]
summary(dt$ipw_pkg)
summary(dt$ipw)

outcome_model_iptw <- glm(outcome.MAKE ~COVID, data = dt, weights = ipw,family = "binomial")
get_coef(outcome_model_iptw)
outcome_model_iptw_pkg <- glm(outcome.MAKE ~COVID, data = dt, weights = ipw_pkg,family = "binomial")
get_coef(outcome_model_iptw_pkg)

## SMR
dt[,smrw:=ifelse(COVID==1,1,psvalue/(1-psvalue))]
summary(dt$smrw)
outcome_model_smrw <- glm(outcome.MAKE ~COVID, data = dt, weights = smrw,family = "binomial")
get_coef(outcome_model_smrw)

