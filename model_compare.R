library(dplyr)
library(nlme)
library(MuMIn)

## Function to choose best linear model -------------------
fit_lm <- function(dat) {
  # Linear model with normal error
  linear_norm <- 
    nlme::gls(y ~ time, 
              data = dat, 
              na.action = na.omit)
  
  # Linear model with AR1 error
  linear_ar1 <- 
    nlme::gls(y ~ time, 
              data = dat, 
              na.action = na.omit,
              correlation = nlme::corAR1(form = ~time))
  
  # Polynomial model with normal error
  dat$time2 <- dat$time^2
  poly_norm <- 
    nlme::gls(y ~ time + time2, 
              data = dat, 
              na.action = na.omit)
  
  # Polynomial model with AR1 error
  poly_ar1 <- 
    nlme::gls(y ~ time + time2, 
              data = dat, 
              na.action = na.omit,
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
                             c(coef(linear_norm),NA),
                             c(coef(linear_ar1),NA)))
  
  best_lm <-
    df_aicc %>%
    dplyr::filter(aicc == min(aicc))
  
  return(best_lm)
}

## Function to choose best GAM ----------------------------
fit_gam <- function(dat) {
  # GAM with normal error
  gam_norm <- 
    mgcv::gam(y ~ s(time), 
              data = dat, 
              na.action = na.omit)
  
  # GAM with AR1 error
  gam_ar1 <- 
    mgcv::gamm(y ~ s(time), 
               data = dat, 
               na.action = na.omit,
               correlation = corAR1(form = ~time))
  
  # Calculate AICs for all models
  df_aicc <- 
    data.frame(model = c("gam_norm",
                         "gam_ar1"),
               aicc  = c(AICc(gam_norm),
                         AICc(gam_ar1)),
               coefs = rbind(coef(gam_norm),
                             coef(gam_ar1$lme)))
  
  best_gam <-
    df_aicc %>%
    dplyr::filter(aicc == min(aicc))
  
  return(best_gam)
}


## Demo the functions -------------------------------------
dat <- data.frame(y = 0:49 + rnorm(50),
                  time = 0:49)

plot(y ~ time, dat, type = "b")

fit_lm(dat = dat)

fit_gam(dat = dat)
