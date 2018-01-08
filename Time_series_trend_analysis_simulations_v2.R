#Geret DePiper
#Simulations to test robusntess of Mann-Kendall, prewhitening, linear GLM, and GAM to assessing trends in time series
#December 14, 2017
#10, 20, 30 years
#AR(1) ARIMA(1,0,1), linear trend, quadratic trend, no trend
PKG <-c("Kendall",'zyp')

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p)
    require(p,charcter.only=TRUE)  }
}

set.seed(436)

n = 10 #number of simulations
x <- 30 #number of periods

LTREND1 <- 0.4232083+(0.004941981*c(1:x)) #Linear trend
LTREND2 <- 0.3503154+(0.001059701*c(1:x))
LTREND3 <- 0.4961013+(0.00882426*c(1:x))
AR1 <- list(ar = 0.4354565)
AR2<- list(ar = 0.8)
NOAR <- list()

#Adding placeholders for simulated data
LINEAR1_AR1 <- NULL 
LINEAR1_AR2 <- NULL 
LINEAR2_AR1 <- NULL 
LINEAR2_AR2 <- NULL
LINEAR3_AR1 <- NULL 
LINEAR3_AR2 <- NULL 

NOTREND_AR1 <- NULL
NOTREND_AR2 <- NULL
NOTREND_NOAR <- NULL

LINEAR1_NOAR <- NULL
LINEAR2_NOAR <- NULL
LINEAR3_NOAR <- NULL

NOTREND_AR1_RESULTS <- NULL
NOTREND_NOAR_RESULTS <- NULL
NOTREND_AR2_RESULTS <- NULL
LINEAR1_AR1_RESULTS <- NULL
LINEAR1_AR2_RESULTS <- NULL
LINEAR2_AR1_RESULTS <- NULL
LINEAR2_AR2_RESULTS <- NULL
LINEAR3_AR1_RESULTS <- NULL
LINEAR3_AR2_RESULTS <- NULL
LINEAR1_NOAR_RESULTS <- NULL
LINEAR2_NOAR_RESULTS <- NULL
LINEAR3_NOAR_RESULTS <- NULL

#initializing simulations
for (i in 1:n) {
#Generating ar(1) simulations
  for (k in c('AR1','AR2','NOAR')){
TEMP1 <- arima.sim(get(k),n=x,rand.gen=rnorm)
LTEMP1 <- TEMP1+LTREND1
LTEMP2 <- TEMP1+LTREND2
LTEMP3 <- TEMP1+LTREND3
assign(paste0('LINEAR1_',k,sep=""),rbind(get(paste0('LINEAR1_',k,sep="")),LTEMP1))
assign(paste0('LINEAR2_',k,sep=""),rbind(get(paste0('LINEAR2_',k,sep="")),LTEMP2))
assign(paste0('LINEAR3_',k,sep=""),rbind(get(paste0('LINEAR3_',k,sep="")),LTEMP3))

assign(paste0('NOTREND_',k, sep=""),rbind(get(paste0('NOTREND_',k, sep="")), TEMP1))

for (y in c('TEMP1','LTEMP1','LTEMP2','LTEMP3')){
  #Testing 30 year series with prewhitening Mann-Kendall technique
TEMP_TEST1 <- zyp.trend.vector(get(y),method='yuepilon')
  TEMP_P1 <- unlist(TEMP_TEST1[6])
  TEMP_TREND1 <- unlist(TEMP_TEST1[3])
  T_1 <- cbind(TEMP_P1,TEMP_TREND1)
  #30 year standard Mann-Kendall
TEMP_TEST12 <- MannKendall(get(y))
T_2 <- as.double(unlist(TEMP_TEST12[2]))

for (j in c(10,20)) {
  #Testind 20 & 10 year series with prewhitening Mann-Kendall technique
  TEMP_TEST2 <- zyp.trend.vector(get(y)[j:x],method='yuepilon')
  TEMP_P2 <- unlist(TEMP_TEST2[6])
  TEMP_TREND2 <- unlist(TEMP_TEST2[3])
  T_1 <- cbind(T_1,TEMP_P2,TEMP_TREND2)
  #Now standard Mann_Kendall
  TEMP_TEST22 <- MannKendall(get(y)[j:x])
  TEMP_P22 <- as.double(unlist(TEMP_TEST22[2]))
  T_2 <- cbind(T_2,TEMP_P22)
}
  colnames(T_1) <- c('p_Val30_pw','Slope30_pw','p_Val20_pw','Slope20_pw','p_Val10_pw','Slope_10_pw')
  colnames(T_2) <- c('p_Val30_mk','p_Val20_mk','p_Val10_mk')
  
if (y=='TEMP1' & k!='NOAR') {assign(paste0('NOTREND_',k,'_RESULTS',sep=""),rbind(get(paste0('NOTREND_',k,'_RESULTS',sep="")),cbind(T_1,T_2)))}
  else if (y=='TEMP1' & k=='NOAR') {assign(paste0('NOTREND_',k,'_RESULTS',sep=""),rbind(get(paste0('NOTREND_',k,'_RESULTS',sep="")),cbind(T_1,T_2)))}
  else if (y=='LTEMP1') {assign(paste0('LINEAR1_',k,'_RESULTS',sep=""),rbind(get(paste0('LINEAR1_',k,'_RESULTS',sep="")),cbind(T_1,T_2)))}
  else if (y=='LTEMP2') {assign(paste0('LINEAR2_',k,'_RESULTS',sep=""),rbind(get(paste0('LINEAR2_',k,'_RESULTS',sep="")),cbind(T_1,T_2)))}
  else if (y=='LTEMP3') {assign(paste0('LINEAR3_',k,'_RESULTS',sep=""),rbind(get(paste0('LINEAR3_',k,'_RESULTS',sep="")),cbind(T_1,T_2)))}
  }
#rm(TEMP1,LTEMP1,QTEMP1)
}
}