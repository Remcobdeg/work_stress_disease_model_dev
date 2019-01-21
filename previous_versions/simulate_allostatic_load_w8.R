
##############################################
# This script ..............................
# # ----CREATING CORTISOL PROFILE BY SIMULATING SUDDEN CORTISOL PRODUCTION BURST USING SIMPLE FORMULA----
# ..............................
# ..............................
##############################################

# The circardian cortisol profile will be constructed as a periodically varying cortisol production of the HPA-axis,
# minus decay of cortisol using:

load_from_file <- F 
save_settings <- F

#setwd("/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simsettings")
setwd("D:/Surfdrive/BS28A - Major project/simulations/simsettings")

####### LOAD SIMULATION SETTINGS FROM FILE ########

if (load_from_file == T){
  settings <- read.csv(file.choose())}

# define work vs. non-work days per person (not include in loaded file)
workday <- rep( c( rep( c(rep(T,5),rep(F,2)),12 ), rep(F,7) ), 20) #12 working weeks with 2 days weekend, then a week of holidays, repeated 20 times

######## DEFINE SIMULATION SETTINGS MANUALLY ########

if (load_from_file == F){

  #settings <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL))
  settings <- data.frame(dt = NA) #initialize dataframe
  
  ## simulation characteristics
  settings$dt <- 1/10 #every 6 minutes
  settings$days <- 1
  settings$N_people <- 250
  settings$N <- with(settings, 24/dt*days)
  settings$integration <- "RK4"
  
  ## define population characteristics
  settings$halftime_mean <- 80/60 #half time of cortisol is 80 minutes or 4/3hour
  settings$halftime_sd <- 8/60
  settings$buffer_capacity_mean <- 15
  settings$buffer_capacity_sd <- 0
  settings$r_recovery_mean <- .3 #rate at which the buffer restores
  settings$r_recovery_sd <- 0
  settings$tolerance_mean <- 1000
  settings$tolerance_sd <- 0
  
  #circadian pulse
  settings$circ_peaktime_mean <- 6 
  settings$circ_intensity_mean <- 0 #1
  settings$circ_peaktime_sd <- 1
  settings$circ_intensity_sd <- .2
  settings$multiplier_circ_intensity_sd_pp <- 0
  
  #awakening impulse
  settings$wake_delay_mean <- 1/2; settings$wake_delay_sd <- 0 #time after onset of the awakening impulse that people (report) waking up 
  settings$wake_intensity_mean <- 12; settings$wake_intensity_sd <- 6
  settings$wake_duration_mean <- 1; settings$wake_duration_sd <- 0 #assume little variation
  settings$wake_time_mean <- 7; settings$wake_time_sd <-0 #people wake up at 7AM on average with a sd of 20 minutes
  #these settings define the variation of means between individuals. However, they say nothing about the intra-individual
  #variation. Therefore, we need also define the intra-individual variation:
  #We will let the intra-individual variation depend on the individual mean
  settings$multiplier_wake_delay_sd_pp <- 0
  settings$multiplier_wake_intensity_sd_pp <- 0
  settings$multiplier_wake_duration_sd_pp <- 0
  settings$wake_time_sd_pp <- 0 #not use a multiplier here
  # wake_delay_sd_pp <- rep(0,N_people)
  # wake_intensity_sd_pp <- rep(0,N_people)
  # wake_duration_sd_pp <- rep(0,N_people)
  # wake_time_sd_pp <- rep(0,N_people) 
  
  #acute stressors
  settings$acute_f_mean <- 10; settings$acute_f_sd <- 5 #frequency of acute stressors
  settings$acute_intensity_mean <- 10; settings$acute_intensity_sd <- 0 #intensity of acute stressors
  settings$acute_duration_mean <- 1/40; settings$acute_duration_sd <- 0 #duration of acute stressors
  settings$work_start <- 8
  settings$work_end <- 17
  #these settings define the variation of means between individuals. However, they say nothing about the intra-individual
  #variation. Therefore, we need also define the intra-individual variation:
  #We will let the intra-individual variation depend on the individual mean
  settings$multiplier_acute_f_sd_pp <- 0
  settings$multiplier_acute_intensity_sd_pp <- 0
  settings$multiplier_acute_duration_sd_pp <- 0
}

## define numbers for set.seed

#v <- data.frame(v = sample(1:999, 100, replace = F)); write.csv(v, 'seeds.csv'); v <- v$v
v <- read.csv('seeds.csv')[,2]

## load packages #####

library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")
#library("reshape")


##### define functions #####


#BLOCK SHAPED BURST
blockburst <- function(duration = 2, start = 6, N = 1440, dt = 1/60, intensity = 1){
  Y <- rep(0,N)
  Y[ round(start/dt) : round((start+duration)/dt)] <- intensity
  return(Y)
}

#SINUS SHAPED CIRCADIAN FLUCTUATION WITH PEAK AT 6AM
circ <- function(N=24/dt, intensity = 1, peak_time = 6){
  Y <- rep(NA,N) #initialize
  for (t in 1:N){
    Y[t] <- sin( 2*pi * (t/N - (peak_time - 6)/24) ) +1 # create the sinus shape (fluctuates from 0 to 2)
  }
  Y <- Y * intensity #scale the signal to the requested intensity
  return(Y)
}

#CORTISOL & ALLOSTATIC LOAD FUNCTION
ALT <- function(HPA, dt, halftime, buffer_capacity, r_recovery, integration = c("Euler", "RK4")[1]){
  #HPA is the vector that describes cortisol production by the HPA-axis
  #dt is the time interval between samples
  #halftime is the cortisol half-time
  #!! note: make sure that dt and halftime are given in the same unit, e.g. h
  
  #calculate decay rate
  #standard decay function: dY/dt = Y*r (decay is a function of the current concentration)
  #one can rewrite this as: dY/Y = r*dt. Integrating over time gives ln(Y(t))-ln(Y(0)) = r*t or Y(t)/Y(0) = exp(r*t) or Y(t) = Y(0)*exp(r*t)
  #filling in: Y(t=1/2) = 1/2*Y(0) = Y(0)*exp(r*t_halftime) or -ln(2) = r*t_halftime or r = -ln(2)/t_halftime
  r <- -log(2)/halftime
  
  #calculate Y[t+1] from dY/dt = Y*r and additional input (cortisol produced by the HPA axis)
  #dY/dt = Y*r + input or dY = (Y*r + input)*dt.
  #this can be integrated numerically using Euler: Y[t+1] = (Y[t]*r + input[t])*dt + Y[t], but also using RK4 (see below)
  
  #initiate vector C. Note: C(0) = 0 + HPA[1]
  N <- length(HPA)
  C <- as.numeric(c(HPA[1], rep(NA,N-2)))
  Dload <- as.numeric(c(0, rep(NA,N-2))) #initialize
  Aload <- as.numeric(c(0, rep(NA,N-2))) #initialize
  
  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           for (t in 1:(N-1)){
             
             C[t+1] <- Euler(Y1 = C[t], func = lingrowth, param = c(r, HPA[t]), h = dt) 
             
             Dload[t+1] <- Euler(Y1 = Dload[t], func = lingrowth, param = c(-r_recovery, C[t]), h = dt) 
             #Dload[t+1] <- Dload[t] - Dload[t] * r_recovery * dt + C[t] * dt 
             
             if (Dload[t] > buffer_capacity){
               Aload[t+1] <- Aload[t] + (Dload[t] - buffer_capacity)*dt}
             else {Aload[t+1] <- Aload[t]}
           },
         
         "RK4" =   # Runge-Kutta 4th Order Integration
         for (t in 1:(N-1)){
           
           C[t+1] <- RK4(Y1 = C[t], func = lingrowth, param = c(r, HPA[t]), h = dt) #RK4 numerical integration
           
           Dload[t+1] <- RK4(Y1 = Dload[t], func = lingrowth, param = c(-r_recovery, C[t]), h = dt) 
           
           if (Dload[t] > buffer_capacity){Aload[t+1] <- Aload[t] + (Dload[t] - buffer_capacity)*dt}
           else {Aload[t+1] <- Aload[t]}
         }
         )

  return(data.frame(C, Dload, Aload))  
  
}

Euler <- function(Y1, func, param, h){
  ## This function performs one iteration using the Euler Integration

  # Iterative process
  Y2 <- Y1 + func(Y1, param)*h
  
  return(Y2)
}

RK4 <- function(Y1, func, param, h){
  ## This function performs one iteration using the Runge-Kutta 4th Order Integration

  k1 <- func(Y1, param)
  k2 <- func(Y1 + 1/2*h*k1, param)
  k3 <- func(Y1 + 1/2*h*k2, param)
  k4 <- func(Y1 + h*k3, param)
  
  # Iterative process
  Y2 <- Y1 + (1/6)*h*(k1+2*k2+2*k3+k4)
  
  return(Y2)
}

## This function performs one iteration of linear growth
lingrowth <- function(Y1,param){
  r <- param[1]; input <- param[2]
  Y2 <- Y1*r + input
  return(Y2)
}



#--------SIMULATE---------

## time stamp the start of the script execution, so we can determine how long the whole thing took
start.time <- Sys.time()



### SOME TRY OUT, IGNORE THIS ###
# #considering a positively skewed distribution of f
# f_median <- 10; f_dev <- 1/5 #a positively skewed distribution of f is assumed
# f_median_pp <- rnorm( N_people, mean = f_median**1/3, sd = (f_median**1/3)*f_dev ) **3  #assuming thas the distirbution of stressor frequencies looks like a strongly positively skewed normal distribution
# densityplot(f_median_pp)
### END OF TRY OUT ###


##pack it in a dataframe that summarizes information for each person per day 
data.day <- data.frame()

#create a variable with person codes
simperson <- c()
for (j in 1:settings$N_people){simperson <- c(simperson, rep(as.character(paste0('SimPerson_',j)), settings$days))} #person variable
data.day <- data.frame(simperson)
rm(simperson)


#add day information
data.day$day <- with(settings, rep(1:days,N_people))

#add holiday information

data.day$workday <- rep(workday[1:settings$days], settings$N_people)
# workday <- rep(T, settings$days)
# workday[ seq(6,settings$days,7) ] <- F #6th day of the week is not a working day
# workday[ seq(7,settings$days,7) ] <- F #7th day of the week is not a working day
# for(i in seq(12,settings$days/7,13)){workday[(i*7+1):(i*7+7)]<-F} #every 13th week is a holiday
# workday <- workday[1:settings$days] #the previous line may have elongated the vector. Make sure it has the appropriate length
rm(workday)


## create people from population characteristics

#create data frame to store all information about individuals
simpopulation <- data.frame(simperson = as.character(paste0('SimPerson_',1:settings$N_people)) )
# simperson <- c()
# for (j in 1:settings$N_people){simperson <- c(simperson, rep(as.character(paste0('SimPerson_',j)), settings$days))} #person variable
# simpopulation <- as.data.frame(simperson)
# rm(simperson)

set.seed(v[1]); simpopulation$halftime_pp <- abs( rnorm(settings$N_people, mean = settings$halftime_mean, sd = settings$halftime_sd) )
set.seed(v[2]); simpopulation$buffer_capacity_pp <- abs( rnorm(settings$N_people, mean = settings$buffer_capacity_mean, sd = settings$buffer_capacity_sd) ) 
set.seed(v[3]); simpopulation$r_recovery_pp <- abs( rnorm(settings$N_people, mean = settings$r_recovery_mean, sd = settings$r_recovery_sd) ) 


##awakening impulse
#indivual mean & within subject variation (same as 'intra-individual variation')
set.seed(v[4]); simpopulation$wake_delay_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_delay_mean, sd = settings$wake_delay_sd) ) #can have only positive values
set.seed(v[5]); simpopulation$wake_delay_sd_pp <- simpopulation$wake_delay_mean_pp * settings$multiplier_wake_delay_sd_pp 
set.seed(v[6]); simpopulation$wake_intensity_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_intensity_mean, sd = settings$wake_intensity_sd) ) #can have only positive values
set.seed(v[7]); simpopulation$wake_intensity_sd_pp <- simpopulation$wake_intensity_mean_pp * settings$multiplier_wake_intensity_sd_pp 
set.seed(v[8]); simpopulation$wake_duration_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_duration_mean, sd = settings$wake_duration_sd) ) #can have only positive values
set.seed(v[9]); simpopulation$wake_duration_sd_pp <- simpopulation$wake_duration_mean_pp * settings$multiplier_wake_duration_sd_pp 
set.seed(v[10]); simpopulation$wake_time_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_time_mean, sd = settings$wake_time_sd) ) #can have only positive values
set.seed(v[11]); simpopulation$wake_time_sd_pp <- settings$wake_time_sd_pp #already defined above, without using a multiplier

#circadian pulse
set.seed(v[12]); simpopulation$circ_peaktime_mean_pp <-rnorm(settings$N_people, mean = settings$circ_peaktime_mean, sd = settings$circ_peaktime_sd)  
set.seed(v[13]); simpopulation$circ_intensity_mean_pp <- abs( rnorm(settings$N_people, mean = settings$circ_intensity_mean, sd = settings$circ_intensity_sd) )  
set.seed(v[14]); simpopulation$circ_intensity_sd_pp <- simpopulation$circ_intensity_mean_pp * settings$multiplier_circ_intensity_sd_pp 

##acute stressors
#indivual mean
set.seed(v[15]); simpopulation$acute_f_mean_pp <- abs( rnorm(settings$N_people, mean = settings$acute_f_mean, sd = settings$acute_f_sd) ) #can have only positive values
set.seed(v[16]); simpopulation$acute_f_sd_pp <- simpopulation$acute_f_mean_pp * settings$multiplier_acute_f_sd_pp
set.seed(v[17]); simpopulation$acute_intensity_mean_pp <-  abs( rnorm(settings$N_people, mean = settings$acute_intensity_mean, sd = settings$acute_intensity_sd) ) #can have only positive values
set.seed(v[18]); simpopulation$acute_intensity_sd_pp <- simpopulation$acute_intensity_mean_pp * settings$multiplier_acute_intensity_sd_pp
set.seed(v[19]); simpopulation$acute_duration_mean_pp <-  abs( rnorm(settings$N_people, mean = settings$acute_duration_mean, sd = settings$acute_duration_sd) ) #can have only positive values
set.seed(v[20]); simpopulation$acute_duration_sd_pp <- simpopulation$acute_duration_mean_pp * settings$multiplier_acute_duration_sd_pp
#within subject variation (same as 'intra-individual variation')


## Create HPAA stimulation based on people characteristics

##intialize
data <- with(settings, data.frame(
  simperson=as.character(c(rep(NA,N*N_people))),
  time=as.numeric(c(rep(NA,N*N_people))),
  HPA=as.numeric(c(rep(NA,N*N_people))),
  cortisol=as.numeric(c(rep(NA,N*N_people))),
  Dload=as.numeric(c(rep(NA,N*N_people))),
  Aload=as.numeric(c(rep(NA,N*N_people))),
  stringsAsFactors=FALSE))

#we want to store the wake time for later use
data.day$wake_time <- NA #with(settings, matrix(NA, days, N_people))

for (j in seq(settings$N_people)){
  
  interval <- (1+(j-1)*settings$N):(j*settings$N) #set the sample range that belongs to the simulated person that is concerned
  
  #create acute stress responses
  HPA_s <- c() #initialize
  for (k in seq(settings$days)){
    f <- round( abs( rnorm(1, mean = simpopulation$acute_f_mean_pp[j], sd = simpopulation$acute_f_sd_pp[j]) )) #can have only positive values
    if (f>0 & data.day$workday[k] ==T){
      HPA_s_day <- matrix(NA,24/settings$dt,f) #initialize
      #create each individual stressor
      for (i in seq(f)){
        
        set.seed(v[21]); duration <- abs(rnorm(1, mean = simpopulation$acute_duration_mean_pp[j], sd = simpopulation$acute_duration_sd_pp[j]))
        set.seed(v[22]); start <- runif(1,settings$work_start,settings$work_end)
        set.seed(v[23]); intensity <- abs(rnorm(1, mean = simpopulation$acute_intensity_mean_pp[j], sd = simpopulation$acute_intensity_sd_pp[j]))
        #note that the absolute is taken to prevent against occasional negative values
        
        HPA_s_day[,i] <- blockburst(duration = duration, start = start, N = 24/settings$dt, dt = settings$dt, intensity = intensity)
      }
      HPA_s_day <- apply(HPA_s_day,1,max) #if stress impulses overlap, the strongest impulse is used. They are not summed
    } else{HPA_s_day <- rep(0,24/settings$dt)} #if f <= 0 it produces a row of zeros
    HPA_s <- c(HPA_s, HPA_s_day) #pasting the days together
  }
  
  #create awakening inpuls
  HPA_a <- c() #initialize
  for (k in seq(settings$days)){
    
    set.seed(v[24]); duration <- abs(rnorm(1, mean = simpopulation$wake_duration_mean_pp[j], sd = simpopulation$wake_duration_sd_pp[j]))
    set.seed(v[25]); wake <- abs( rnorm(1, mean = simpopulation$wake_time_mean_pp[j], sd = simpopulation$wake_time_sd_pp[j]))
    set.seed(v[26]); start <- abs( wake - rnorm(1, mean = simpopulation$wake_delay_mean_pp[j], sd = simpopulation$wake_delay_sd_pp[j]) )
    set.seed(v[27]); intensity <- abs(rnorm(1, mean = simpopulation$wake_intensity_mean_pp[j], sd = simpopulation$wake_intensity_sd_pp[j]))
    #note that the absolute is taken to prevent against occasional negative values
    
    HPA_a_day <- blockburst(duration = duration, start = start, N = 24/settings$dt, dt = settings$dt, intensity = intensity)
    
    HPA_a <- c(HPA_a, HPA_a_day) #pasting the days together
    
    #store the wake time for later use
    data.day$wake_time[k+(j-1)*settings$days] <- wake
  }
  
  #create a circadian sinus fluctuation
  HPA_c <- c() #initialize
  for (k in seq(settings$days)){
    
    set.seed(v[28]); intensity <- abs(rnorm(1, mean = simpopulation$circ_intensity_mean_pp[j], sd = simpopulation$circ_intensity_sd_pp[j]))
    #note that the absolute is taken to prevent against occasional negative values
    
    HPA_c_day <- circ(N = 24/settings$dt, intensity = intensity, peak_time = simpopulation$circ_peaktime_mean_pp[j])
    
    HPA_c <- c(HPA_c, HPA_c_day) #pasting the days together
  }
  
  
  #define total cortisol release
  data$time[interval] <- settings$dt*(1:settings$N)
  data$HPA[interval] <- (HPA_a + HPA_s + HPA_c) #; rm(HPA_a_day, HPA_s_day, HPA_c_day, HPA_a, HPA_s, HPA_c)
  data$simperson[interval] <- as.character(paste0('SimPerson_',j))
  rm(interval)
  
  # ########### TEMP ################
  # plot(cbind(
  # ts(HPA_s+HPA_a+HPA_c), ts(HPA_s), ts(HPA_a), ts(HPA_c)
  # ))
  # #################################
}

data$simperson <- as.factor(data$simperson) #change to a factor for later use


# #create dataframe that describes cortisol in the bloodstream for each person
# Cortisol <- data.frame(matrix(NA, N, N_people)) #initialize
# for (j in seq(N_people)){
#   Cortisol[1:N,j] <- cortisol(HPA = HPA[,j], dt = dt, halftime = halftime_pp)
# }


#create dataframe that describes cortisol in the bloodstream for each person
for (j in seq(settings$N_people)){
  interval <- (1+(j-1)*settings$N):(j*settings$N) #set the sample range that belongs to the simulated person that is concerned
  data$cortisol[interval] <- ALT(HPA = data$HPA[interval], dt = settings$dt, halftime = simpopulation$halftime_pp[j], buffer_capacity = simpopulation$buffer_capacity_pp[j], r_recovery = simpopulation$r_recovery_pp[j], integration = settings$integration)$C
  data$Dload[interval] <- ALT(HPA = data$HPA[interval], dt = settings$dt, halftime = simpopulation$halftime_pp[j], buffer_capacity = simpopulation$buffer_capacity_pp[j], r_recovery = simpopulation$r_recovery_pp[j], integration = settings$integration)$Dload
  data$Aload[interval] <- ALT(HPA = data$HPA[interval], dt = settings$dt, halftime = simpopulation$halftime_pp[j], buffer_capacity = simpopulation$buffer_capacity_pp[j], r_recovery = simpopulation$r_recovery_pp[j], integration = settings$integration)$Aload

  ########### TEMP ################
  plot(ts(data$cortisol[interval]))
  plot(ts(data$Aload[interval]))
  #################################
}

########### TEMP ################

#plot curves of 5 random people, first 2 days
set.seed(v[29])
if (settings$N_people >= 5){
  psample <- sample(simpopulation$simperson, 5, replace = FALSE, prob = NULL) #5 random people
} else {psample <- simpopulation$simperson}

# #HPA output
# ggplot(subset(data, ((simperson %in% psample) & (time <= 48))), aes(x=time, y=HPA)) +
#   geom_line() + facet_wrap( ~ simperson, ncol=1) +
#   ggtitle("Cortisol production (= HPA-axis activity) of 5 randomly selected people (first 2 days)") +
#   ylab("Cortisol production [mmol/L]") + xlab("time [h]")
# 
# #cortisol profile
# ggplot(subset(data, ((simperson %in% psample) & (time <= 48))), aes(x=time, y=cortisol)) +
#   geom_line() + facet_wrap( ~ simperson, ncol=1) +
#   ggtitle("Bloodstream cortisol of 5 randomly selected people (first 2 days)") +
#   ylab("Cortisol [mmol/L]") + xlab("time [h]")

#dynamic load profile
ggplot(subset(data, ((simperson %in% psample) & (time <= 48))), aes(x=time, y=Dload)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Dynamic load of 5 randomly selected people (first 2 days)") +
  ylab("Dynamic load") + xlab("time [h]")

#allostatic load profile
ggplot(subset(data, ((simperson %in% psample) & (time <= 48))), aes(x=time, y=Aload)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Allostatic load of 5 randomly selected people (first 2 days)") +
  ylab("Allostatic load") + xlab("time [h]")
# 
# 
# 
# # plot quantiles
# 
# df.qt <- data %>% group_by(time) %>% summarise(
#   Q05=quantile(cortisol,probs=.05),
#   Q25=quantile(cortisol,probs=.25),
#   Q50=quantile(cortisol,probs=.50),
#   Q75=quantile(cortisol,probs=.75),
#   Q95=quantile(cortisol,probs=.95))
# 
# with(df.qt, ts.plot( 
#   ts(Q05, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
#   ts(Q25, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
#   ts(Q50, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
#   ts(Q75, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
#   ts(Q95, start = 0, end = 24-settings$dt, frequency = 1/settings$dt)
# ))
# 
# 
# plot quantiles since awakening

#make a new dataframe that includes only the first day
cortisol.sinceawakening <- c(NULL)
time.sinceawakening <- c(NULL)
for (j in seq(settings$N_people)){
  #per person we want the data from the first day from awakening until end of the day
  #wake time for person one is data.day$wake_time[1], for person two this is data.day$wake_time[1+settings$days]
  interval <- which( (data$time > data.day$wake_time[1+(j-1)*settings$days]) & (data$simperson == data$simperson[j*settings$N])
                     & data$time < 24)
  cortisol.sinceawakening <- c(cortisol.sinceawakening, data$cortisol[interval])
  time.sinceawakening <- c(time.sinceawakening , seq_along(interval) * settings$dt )
}

data.sinceawakening <- data.frame(cortisol.sinceawakening, time.sinceawakening)

df.qt <- data.sinceawakening %>% group_by(time.sinceawakening) %>% summarise(
  Q05=quantile(cortisol.sinceawakening,probs=.05),
  Q25=quantile(cortisol.sinceawakening,probs=.25),
  Q50=quantile(cortisol.sinceawakening,probs=.50),
  Q75=quantile(cortisol.sinceawakening,probs=.75),
  Q95=quantile(cortisol.sinceawakening,probs=.95))

######### TRYOUT TO SHAPE TO LONG FORMAT ############
# df.qt$dummy <- paste0('dummy',df.qt$time.sinceawakening)
# data2 <- melt(df.qt, id = names(df.qt[1]), measured = names(df.qt[2:6]))
# #reshape(df.qt, idvar = "dummy", varying = c("Q05", "Q25", "Q50", "Q75", "Q95"), timevar = "Quantile" , v.names = "Quantile", direction = "long")

with(df.qt, ts.plot(
  ts(Q05, start = 0, end = (length(Q05)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q25, start = 0, end = (length(Q25)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q50, start = 0, end = (length(Q50)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q75, start = 0, end = (length(Q75)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q95, start = 0, end = (length(Q95)-1)*settings$dt, frequency = 1/settings$dt)
))


#### ------EXPLORING PARAMETER SETTINGS FOR DEVELOPING ALLOSTATIC LOAD------ ######

# #Note: this exploration requires only 1 day to be simulated per person
#
# AUC <- colSums(Cortisol)*dt
# hist(AUC)
# quantile(AUC,probs=.10, na.rm = T)
# quantile(AUC,probs=.50, na.rm = T) #~27.5 --> set buffer capacity to 27.5 to start with
#
# data <- allostatic_load(Cortisol[,1], dt = dt, a = .2, k = 5)
#
# input <- Cortisol[,1]
# a <- .005 #use order of .2
# k <- 5 #use ~10


## visualize distribution of buffer and allostatic load at the end of each day

# for (i in 1:days){
#   hist(as.matrix(Buffer[24/dt*i,]), main = paste0('Distribution of remain buffer for ', N_people, ' people at end of day ', i, ' (capacity = ', buffer_capacity, ')'), xlab = 'Remaining buffer' )
#   hist(as.matrix(Aload[24/dt*i,]), main = paste0('Distribution of allostatic load for ', N_people, ' people at end of day ', i), xlab = 'Allostatic load' )
# }

#hist(cbind( as.matrix(Buffer[24/dt*i,]) , as.matrix(Aload[24/dt*i,])) )

## determine how long the script took
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# ## visualize for each person how their allostatic load increases
# 
# ggplot(data,aes(x=time, y=buffer, colour=simperson)) +
#   geom_line() +
#   theme(legend.position="none") + #remove legend
#   ggtitle(paste("Load buffer variation - N =",N_people)) +
#   ylab("available buffer")
# 
# ## visualize for eacht person how their allostatic load increases
# 
# ggplot(data,aes(x=time, y=Aload, colour=simperson)) +
#   geom_line() +
#   #theme_bw() +
#   theme(legend.position="none") + #remove legend
#   ggtitle(paste("Allostatic load development - N =",N_people)) +
#   ylab("Allostatic load")

# combine data for displaying in data.day dataframe

#add personal mean to data.day
for (i in 1:settings$N_people){
  data.day$acute_f_mean_pp[(i-1)*settings$days + 1:settings$days] <- simpopulation$acute_f_mean_pp[i]
}


########### TEMP ################
# 
# #add allostatic load at noon each day to data.day
# data.day$allostatic_load <- data$Aload[with(settings, seq(12/dt,N*N_people,24/dt))]
# 
# ## plot stressor - allostatic load relationship
# 
# #choose subset of days to plot
# if(settings$days<=10){plotdays <- 1:settings$days
# } else {plotdays <- seq( as.integer(settings$days/10), settings$days, as.integer(settings$days/10) )}
# 
# ggplot(subset(data.day, day %in% plotdays), aes(x=acute_f_mean_pp, y=allostatic_load)) +
#   geom_point() + facet_wrap( ~ day, ncol=5) + geom_smooth(method=lm, se=F) + 
#   #theme(legend.position="none") + #remove legend
#   ggtitle(paste("Allostatic load vs. work stress, N =",settings$N_people)) +
#   ylab("Allostatic load") + xlab("average frequency of stressors")
# 

if (save_settings == T){
  file = paste("Settings", format(Sys.time(),"%Y-%m-%d_%H-%M"),".csv", sep="_")
  write.csv(settings, file)
}




#loop in script die variablen doorgeeft naar de markdown. Of dataset laten wegschrijven die door markdwon
#worden ingelezen. Zoeken naar tutorials generating reports. Ook cheat sheets voor codes

