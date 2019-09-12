library("dplyr")
library("forecast")
library("TSA")
library("tseries")
library("quantmod")

df = read.table("~/Documents/DogProject_Jaz/LocalRscripts/AKC/AKC_breedPopularity_1926thru2005_proportional.txt", header = T)
IBDNe = read.table("~/Documents/DogProject_Jaz/LocalRscripts/IBDNe/AllSitesMerge/boxer_IBDNE_usingIBDSeqIBDSegs_May18ver.ne", header = T)

#Remake time axes
df$TimeToPresent = abs(max(df$Time)-df$Time)
IBDNe$TimeToPresent = IBDNe$GEN*4
  
#new df
comboDF = df %>% select("TimeToPresent", "boxer") %>% mutate(NE = IBDNe$NE[match(TimeToPresent, IBDNe$TimeToPresent)]) %>% arrange(TimeToPresent) #arrange df reverse order so time series is correct

#Make time series data
ts_popularity = ts(comboDF$boxer, start=min(comboDF$TimeToPresent), end=max(comboDF$TimeToPresent)) #yearly data
ts_IBDNe = ts(comboDF$NE, start = min(comboDF$TimeToPresent), end=max(comboDF$TimeToPresent))

#Check whether data is stationary
plot(na.remove(ts_popularity)) #data is not stationary
plot(na.remove(ts_IBDNe)) #data is not stationary

#Make data sets comparable
comparableDF = na.remove(ts.intersect(ts_IBDNe,ts_popularity))

#Check for correlation
ccf(comparableDF[,"ts_IBDNe"],comparableDF[,"ts_popularity"]) #it's there but it could be do to noise 

#Make data stationary 
IBDNe_optimumModel = auto.arima(comparableDF[,"ts_IBDNe"])

#pre-whiten and obtain cross-correlation
prewhiten( comparableDF[,"ts_IBDNe"], comparableDF[,"ts_popularity"], x.model = IBDNe_optimumModel)

#Alternatively
ccf(comparableDF[,"ts_popularity"],comparableDF[,"ts_IBDNe"]) #it's there but it could be do to noise
popularity_optimumModel = auto.arima(comparableDF[,"ts_popularity"])
prewhiten(comparableDF[,"ts_popularity"], comparableDF[,"ts_IBDNe"],  x.model = popularity_optimumModel)


#Multi-variate linear regression
#install.packages("vars") #If not already installed
#install.packages("astsa") #If not already installed
library(vars)
library(astsa)

fitvar1=VAR(comparableDF, p=1, type="both")
summary(fitvar1)


fitvar2=VAR(comparableDF, p=2, type="both")
summary(fitvar2)


#Granger test of causality 
#need to examine lags from VAR and choose appropriatley
grangertest(comparableDF[,"ts_popularity"] ~ comparableDF[,"ts_IBDNe"], order = 6)
grangertest( comparableDF[,"ts_IBDNe"] ~ comparableDF[,"ts_popularity"] , order = 6)
