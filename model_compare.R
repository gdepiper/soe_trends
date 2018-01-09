library(dplyr)
library(nlme)
library(MuMIn)

## Function to choose best linear model -------------------
fit_lm <- function(dat) {
  # Remove missing values first so that all models
  # use the same number of observations (important for AIC)
  dat <- dat %>% dplyr::filter(complete.cases(.))
  
  # Constant model (null model used to calculate 
  # overall p-value)
  constant_norm <-
    nlme::gls(series ~ 1, 
              data = dat)
  
  constant_ar1 <-
    nlme::gls(series ~ 1,
              data = dat,
              correlation = nlme::corAR1(form = ~time))
  
  
  
  # Linear model with normal error
  linear_norm <- 
    nlme::gls(series ~ time, 
              data = dat)
  
  # Linear model with AR1 error
  linear_ar1 <- 
    nlme::gls(series ~ time, 
              data = dat,
              correlation = nlme::corAR1(form = ~time))
  
  # Polynomial model with normal error
  dat$time2 <- dat$time^2
  poly_norm <- 
    nlme::gls(series ~ time + time2, 
              data = dat)
  
  # Polynomial model with AR1 error
  poly_ar1 <- 
    nlme::gls(series ~ time + time2, 
              data = dat,
              correlation = nlme::corAR1(form = ~time))
  
  # Calculate AICs for all models
  df_aicc <- 
    data.frame(model = c("poly_norm", 
                         "poly_ar1",
                         "linear_norm", 
                         "linear_ar1"),
               aicc  = c(AICc(poly_norm),
                         AICc(poly_ar1),
                         AICc(linear_norm),
                         AICc(linear_ar1)),
               coefs = rbind(coef(poly_norm),
                             coef(poly_ar1),
                             c(coef(linear_norm), NA),
                             c(coef(linear_ar1),  NA)),
               # Calculate overall signifiance (need to use
               # ML not REML for this)
               pval = c(anova(update(constant_norm, method = "ML"), 
                              update(poly_norm, method = "ML"))$`p-value`[2],
                        anova(update(constant_ar1, method = "ML"), 
                              update(poly_ar1, method = "ML"))$`p-value`[2],
                        anova(update(constant_norm, method = "ML"), 
                              update(linear_norm, method = "ML"))$`p-value`[2],
                        anova(update(constant_ar1, method = "ML"), 
                              update(linear_ar1, method = "ML"))$`p-value`[2]))
  
  best_lm <-
    df_aicc %>%
    dplyr::filter(aicc == min(aicc))
  
  return(best_lm)
}

## Function to choose best GAM ----------------------------
fit_gam <- function(dat) {
    
  # Remove missing values first so that all models
  # use the same number of observations (important for AIC)
  dat <- dat %>% dplyr::filter(complete.cases(.))
  
  # Break out if less than 15 data points (GAMM runs into numerical problems)
  if (nrow(dat) < 15) {
    best_gam <- 
      data.frame(model = NA,
                 aicc  = NA,
                 coefs = NA,
                 pval  = NA) 
    return(best_gam)
  }
  
  # GAM with normal error
  gam_norm <- 
    mgcv::gam(series ~ s(time), 
              data = dat)
  
  # GAM with AR1 error
  gam_ar1 <- 
    mgcv::gamm(series ~ s(time), 
               data = dat,
               correlation = corAR1(form = ~time))
  
  # Calculate AICs for all models
  df_aicc <- 
    data.frame(model = c("gam_norm",
                         "gam_ar1"),
               aicc  = c(AICc(gam_norm),
                         AICc(gam_ar1)),
               coefs = rbind(coef(gam_norm),
                             coef(gam_ar1$gam)),
               pval  = c(summary(gam_norm)$s.table[, "p-value"],
                         summary(gam_ar1$gam)$s.table[, "p-value"]))
  
  best_gam <-
    df_aicc %>%
    dplyr::filter(aicc == min(aicc))
  
  
  return(best_gam)
}


## Demo the functions -------------------------------------
# dat <- data.frame(series = 0:49 + rnorm(50),
#                   time = 0:49)
# 
# plot(series ~ time, dat, type = "b")
# 
# fit_lm(dat = dat)
# 
# fit_gam(dat = dat)
