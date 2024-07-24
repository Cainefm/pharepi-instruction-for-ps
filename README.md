# Propensity Score techniques in a static treatment pharmacoepidemiology question: a Step-by-Step instruction

![](https://img.shields.io/badge/Language-R-blue?style=flat) ![](https://img.shields.io/badge/progress-45%25-yellowgreen?style=flat)

Contact: [Joe Blais](%5Bhttps://www.bryer.org/%5D(https://www.pharma.hku.hk/en/Our-People/Professoriate/Assistant-Professor/Professor-Joseph-Edgar-BLAIS/Professor-Joseph-Edgar-BLAIS)) ([joeblais\@hku.hk](mailto:joeblais@hku.hk){.email})

# Overview

Propensity score methods have been widely utilized in pharmacoepidemiology studies over the past few decades. The codes introduce the methods in the manuscript for research questions involving static treatments or exposures, excluding time-varying or dynamic treatments.

The project includes code for various approaches: regression adjustment, matching, stratification, and weighting. These techniques help address confounding in observational studies, enhancing the validity of causal inferences.

# Propensity Score Approaches

``` r
head(dt)
```

## Manually calculating the propensity score

``` r
ps_formula <- as.formula(paste0("COVID~",paste(cova,collapse = "+")))
ps_model <- glm(formula = ps_formula,data = dt,family = "binomial")
dt$psvalue <- ps_model$fitted.values

tableone_ps <- print(CreateTableOne(vars = c(cova,"psvalue"),strata = "COVID",data = dt,test = FALSE),smd =T)
as.data.table(tableone_ps,keep.rownames = T)[as.numeric(SMD)>0.1]
```

## Calculating the propensity score using the package

# 1. Regression Adjustment

# 2. Matching

# 3. Stratification

# 4. Weighting
