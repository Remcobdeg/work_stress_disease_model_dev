##############################################
# This script ..............................
# ..............................
# ..............................
# ..............................
##############################################

# The circardian cortisol profile will be constructed as a periodically varying cortisol production of the HPA-axis,
# minus decay of cortisol using:


library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")

# ----CREATING CORTISOL PROFILE BY SIMULATING SUDDEN CORTISOL PRODUCTION BURST USING SIMPLE FORMULA----


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
ALT <- function(HPA, dt, halftime, buffer_capacity, r_recovery){
  #HPA is the vector that describes cortisol production by the HPA-axis
  #dt is the time interval between samples
  #halftime is the cortisol half-time
  #!! note: make sure that dt and halftime are given in the same unit, e.g. h
  
  #removal of cortisol:
  #Assuming that cortisol is removed at constant rate --> C[t+1] = a*C. With a < 1, this means that C decays exponentially.
  #The value of a depends on the sample interval (otherwise, if we would keep 'a' constant, a shorter sampling interval
  #would mean faster decay). Knowing 'a' for a specific sampling interval, we can determine the value of 'a' for different
  #sampling intervals --> a(r) = exp(r * new_interval/known_interval). r can be calculated from combinations between 'a'
  #and 'known_interval'. We know the halftime for cortisol: ~60min. So --> 1/2 = exp(r * known_interval/known_interval) =
  # = exp(r) --> r = ln(2). With dt our new_interval --> a = exp( ln(2) * dt / 60 )
  # C[t+1] = decay + HPA_production --> C[t+1] = a(r)*C[t] + HPA[t] =  C[t+1] = exp( -ln(2) * dt/halftime ) * C[t] + HPA[t]
  
  N <- length(HPA)
  
  #initiate vector C. Note: C(0) = 0 + HPA[1]
  C <- as.numeric(c(HPA[1], rep(NA,N-1)))
  buffer <- as.numeric(c(buffer_capacity, rep(NA,N-2))) #initialize
  a_load <- as.numeric(c(0, rep(NA,N-2))) #initialize
  
  #calculate cortisol values
  for (t in 1:(N-1)){
    C[t+1] <- C[t]*exp( -log(2) * dt/halftime ) + HPA[t]*dt #Euler integration
    
    # ########################### RK4 numerical solution ############################
    # #first rewrite the euler equation to the form Y[n+1] = Y[n] + f(Y[n]) -->
    # #C[t+1] <- C[t] + C[t] * ( exp( -log(2) * dt/halftime ) -1 ) + HPA[t]*dt
    #
    # k1 <- C[t] * ( exp( -log(2) * dt/halftime ) -1 ) #+ HPA[t]*dt
    # 
    # k2 <- (C[t] + 1/2*dt*k1 ) * ( exp( -log(2) * (dt*3/2)/halftime ) -1 ) #+ HPA[t]*dt/2
    # 
    # k3 <- (C[t] + 1/2*dt*k2 ) * ( exp( -log(2) * (dt*3/2)/halftime ) -1 ) #+ HPA[t]*dt/2
    # 
    # k4 <- (C[t] + dt*k3) * ( exp( -log(2) * (dt*2)/halftime ) -1 )# + HPA[t]*dt
    # 
    # # Iterative process
    # C[t+1] <- C[t] + (1/6)*dt*(k1+2*k2+2*k3+k4) + HPA[t]*dt
    # ###############################################################################
    
    buffer[t+1] <- buffer[t] - C[t] * dt + r_recovery * (buffer_capacity - buffer[t]) * dt
    
    if (buffer[t] < 0){a_load[t+1] <- a_load[t] + abs(buffer[t])*dt}
    else {a_load[t+1] <- a_load[t]}
  }
  return(data.frame(C, buffer, a_load))  
  
}





#--------SIMULATE---------

## time stamp the start of the script execution, so we can determine how long the whole thing took
start.time <- Sys.time()

#settings <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL))
settings <- data.frame(dt = NA) #initialize dataframe

## simulation characteristics

settings$dt <- 1/60 #every 1 minute
settings$days <- 1000
settings$N_people <- 10000
settings$N <- with(settings, 24/dt*days)

## define population characteristics

settings$halftime_mean <- 80/60 #half time of cortisol is 80 minutes or 4/3hour
settings$halftime_sd <- 8/60
settings$buffer_capacity_mean <- 15
settings$buffer_capacity_sd <- 1.5
settings$r_recovery_mean <- .3 #rate at which the buffer restores
settings$r_recovery_sd <- .03
settings$tolerance_mean <- 1000
settings$tolerance_sd <- 100

#circadian pulse
settings$circ_peaktime_mean <- 6
settings$circ_intensity_mean <- 0 #1
settings$circ_peaktime_sd <- .5
settings$circ_intensity_sd <- 0
settings$multiplier_circ_intensity_sd_pp <- 0

#awakening impulse
settings$wake_delay_mean <- 1/2; settings$wake_delay_sd <- 1/20 #time after onset of the awakening impulse that people (report) waking up 
settings$wake_intensity_mean <- 15; settings$wake_intensity_sd <- 6
settings$wake_duration_mean <- 1; settings$wake_duration_sd <- .1 #assume little variation
settings$wake_time_mean <- 7; settings$wake_time_sd <-1/3 #people wake up at 7AM on average with a sd of 20 minutes
#these settings define the variation of means between individuals. However, they say nothing about the intra-individual
#variation. Therefore, we need also define the intra-individual variation:
#We will let the intra-individual variation depend on the individual mean
settings$multiplier_wake_delay_sd_pp <- 1/4
settings$multiplier_wake_intensity_sd_pp <- 1/4
settings$multiplier_wake_duration_sd_pp <- 1/4
settings$wake_time_sd_pp <- 1/3 #not use a multiplier here
# wake_delay_sd_pp <- rep(0,N_people)
# wake_intensity_sd_pp <- rep(0,N_people)
# wake_duration_sd_pp <- rep(0,N_people)
# wake_time_sd_pp <- rep(0,N_people) 

#acute stressors
settings$acute_f_mean <- 5; settings$acute_f_sd <- 1 #frequency of acute stressors
settings$acute_intensity_mean <- 6; settings$acute_intensity_sd <- 1 #intensity of acute stressors
settings$acute_duration_mean <- 1/20; settings$acute_duration_sd <- 1/60 #duration of acute stressors
settings$work_start <- 8
settings$work_end <- 17
#these settings define the variation of means between individuals. However, they say nothing about the intra-individual
#variation. Therefore, we need also define the intra-individual variation:
#We will let the intra-individual variation depend on the individual mean
settings$multiplier_acute_f_sd_pp <- 1/2
settings$multiplier_acute_intensity_sd_pp <- 1/2
settings$multiplier_acute_duration_sd_pp <- 1


### SOME TRY OUT, IGNORE THIS ###
# #considering a positively skewed distribution of f
# f_median <- 10; f_dev <- 1/5 #a positively skewed distribution of f is assumed
# f_median_pp <- rnorm( N_people, mean = f_median**1/3, sd = (f_median**1/3)*f_dev ) **3  #assuming thas the distirbution of stressor frequencies looks like a strongly positively skewed normal distribution
# densityplot(f_median_pp)
### END OF TRY OUT ###

## define work vs. non-work days per person
workday <- rep(T, settings$days)
workday[ seq(6,settings$days,7) ] <- F #6th day of the week is not a working day
workday[ seq(7,settings$days,7) ] <- F #7th day of the week is not a working day
for(i in seq(12,settings$days/7,13)){workday[(i*7+1):(i*7+7)]<-F} #every 13th week is a holiday
workday <- workday[1:settings$days] #the previous line may have elongated the vector. Make sure it has the appropriate length

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
data.day$workday <- rep(workday, settings$N_people)
rm(workday)
                     

## create people from population characteristics

#create data frame to store all information about individuals
simpopulation <- data.frame(simperson = as.character(paste0('SimPerson_',1:settings$N_people)) )
# simperson <- c()
# for (j in 1:settings$N_people){simperson <- c(simperson, rep(as.character(paste0('SimPerson_',j)), settings$days))} #person variable
# simpopulation <- as.data.frame(simperson)
# rm(simperson)

simpopulation$halftime_pp <- abs( rnorm(settings$N_people, mean = settings$halftime_mean, sd = settings$halftime_sd) )
simpopulation$buffer_capacity_pp <- abs( rnorm(settings$N_people, mean = settings$buffer_capacity_mean, sd = settings$buffer_capacity_sd) ) 
simpopulation$r_recovery_pp <- abs( rnorm(settings$N_people, mean = settings$r_recovery_mean, sd = settings$r_recovery_sd) ) 


##awakening impulse
#indivual mean & within subject variation (same as 'intra-individual variation')
simpopulation$wake_delay_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_delay_mean, sd = settings$wake_delay_sd) ) #can have only positive values
simpopulation$wake_delay_sd_pp <- simpopulation$wake_delay_mean_pp * settings$multiplier_wake_delay_sd_pp 
simpopulation$wake_intensity_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_intensity_mean, sd = settings$wake_intensity_sd) ) #can have only positive values
simpopulation$wake_intensity_sd_pp <- simpopulation$wake_intensity_mean_pp * settings$multiplier_wake_intensity_sd_pp 
simpopulation$wake_duration_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_duration_mean, sd = settings$wake_duration_sd) ) #can have only positive values
simpopulation$wake_duration_sd_pp <- simpopulation$wake_duration_mean_pp * settings$multiplier_wake_duration_sd_pp 
simpopulation$wake_time_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_time_mean, sd = settings$wake_time_sd) ) #can have only positive values
simpopulation$wake_time_sd_pp <- settings$wake_time_sd_pp #already defined above, without using a multiplier

#circadian pulse
simpopulation$circ_peaktime_mean_pp <-rnorm(settings$N_people, mean = settings$circ_peaktime_mean, sd = settings$circ_peaktime_sd)  
simpopulation$circ_intensity_mean_pp <- abs( rnorm(settings$N_people, mean = settings$circ_intensity_mean, sd = settings$circ_intensity_sd) )  
simpopulation$circ_intensity_sd_pp <- simpopulation$circ_intensity_mean_pp * settings$multiplier_circ_intensity_sd_pp 

##acute stressors
#indivual mean
simpopulation$acute_f_mean_pp <- abs( rnorm(settings$N_people, mean = settings$acute_f_mean, sd = settings$acute_f_sd) ) #can have only positive values
simpopulation$acute_f_sd_pp <- simpopulation$acute_f_mean_pp * settings$multiplier_acute_f_sd_pp
simpopulation$acute_intensity_mean_pp <-  abs( rnorm(settings$N_people, mean = settings$acute_intensity_mean, sd = settings$acute_intensity_sd) ) #can have only positive values
simpopulation$acute_intensity_sd_pp <- simpopulation$acute_intensity_mean_pp * settings$multiplier_acute_intensity_sd_pp
simpopulation$acute_duration_mean_pp <-  abs( rnorm(settings$N_people, mean = settings$acute_duration_mean, sd = settings$acute_duration_sd) ) #can have only positive values
simpopulation$acute_duration_sd_pp <- simpopulation$acute_duration_mean_pp * settings$multiplier_acute_duration_sd_pp
#within subject variation (same as 'intra-individual variation')


## Create HPAA stimulation based on people characteristics

##intialize
#dataframe of HPAA stimulation per person (persons in columns, sample moments in rows)
data <- with(settings, data.frame(sim_Person=as.character(c(rep(NA,N*N_people))),stringsAsFactors=FALSE))
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

        duration <- abs(rnorm(1, mean = simpopulation$acute_duration_mean_pp[j], sd = simpopulation$acute_duration_sd_pp[j]))
        start <- runif(1,settings$work_start,settings$work_end)
        intensity <- abs(rnorm(1, mean = simpopulation$acute_intensity_mean_pp[j], sd = simpopulation$acute_intensity_sd_pp[j]))
        #note that the absolute is taken to prevent against occasional negative values

        HPA_s_day[,i] <- blockburst(duration = duration, start = start, N = 24/settings$dt, dt = 1/60, intensity = intensity)
      }
      HPA_s_day <- apply(HPA_s_day,1,max) #if stress impulses overlap, the strongest impulse is used. They are not summed
    } else{HPA_s_day <- rep(0,24/settings$dt)} #if f <= 0 it produces a row of zeros
    HPA_s <- c(HPA_s, HPA_s_day) #pasting the days together
  }

  #create awakening inpuls
  HPA_a <- c() #initialize
  for (k in seq(settings$days)){

    duration <- abs(rnorm(1, mean = simpopulation$wake_duration_mean_pp[j], sd = simpopulation$wake_duration_sd_pp[j]))
    wake <- abs( rnorm(1, mean = simpopulation$wake_time_mean_pp[j], sd = simpopulation$wake_time_sd_pp[j]))
    start <- abs( wake - rnorm(1, mean = simpopulation$wake_delay_mean_pp[j], sd = simpopulation$wake_delay_sd_pp[j]) )
    intensity <- abs(rnorm(1, mean = simpopulation$wake_intensity_mean_pp[j], sd = simpopulation$wake_intensity_sd_pp[j]))
    #note that the absolute is taken to prevent against occasional negative values

    HPA_a_day <- blockburst(duration = duration, start = start, N = 24/settings$dt, dt = 1/60, intensity = intensity)

    HPA_a <- c(HPA_a, HPA_a_day) #pasting the days together

    #store the wake time for later use
    print(k+(j-1)*settings$days)
    data.day$wake_time[k+(j-1)*settings$days] <- wake
    #################################################
  }
  
  #create a circadian sinus fluctuation
  HPA_c <- c() #initialize
  for (k in seq(settings$days)){
    
    intensity <- abs(rnorm(1, mean = simpopulation$circ_intensity_mean_pp[j], sd = simpopulation$circ_intensity_sd_pp[j]))
    #note that the absolute is taken to prevent against occasional negative values
    
    HPA_c_day <- circ(N = 24/settings$dt, intensity = intensity, peak_time = simpopulation$circ_peaktime_mean_pp[j])
    
    HPA_c <- c(HPA_c, HPA_c_day) #pasting the days together
  }

  
  
  #plot(ts( (HPA_a + HPA_s)), main = paste0('person = ',j))
  
  #define total cortisol release
  data$time[interval] <- settings$dt*(1:settings$N)
  data$HPA[interval] <- (HPA_a + HPA_s + HPA_c); rm(HPA_a_day, HPA_s_day, HPA_c_day, HPA_a, HPA_s, HPA_c)
  data$sim_Person[interval] <- as.character(paste0('SimPerson_',j))
  rm(interval)
}

data$sim_Person <- as.factor(data$sim_Person) #change to a factor for later use


# #create dataframe that describes cortisol in the bloodstream for each person
# Cortisol <- data.frame(matrix(NA, N, N_people)) #initialize
# for (j in seq(N_people)){
#   Cortisol[1:N,j] <- cortisol(HPA = HPA[,j], dt = dt, halftime = halftime_pp)
# }


#create dataframe that describes cortisol in the bloodstream for each person
for (j in seq(settings$N_people)){
  interval <- (1+(j-1)*settings$N):(j*settings$N) #set the sample range that belongs to the simulated person that is concerned
  data$cortisol[interval] <- ALT(HPA = data$HPA[interval], dt = settings$dt, halftime = simpopulation$halftime_pp[j], buffer_capacity = simpopulation$buffer_capacity_pp[j], r_recovery = simpopulation$r_recovery_pp[j])$C
  data$buffer[interval] <- ALT(HPA = data$HPA[interval], dt = settings$dt, halftime = simpopulation$halftime_pp[j], buffer_capacity = simpopulation$buffer_capacity_pp[j], r_recovery = simpopulation$r_recovery_pp[j])$buffer
  data$a_load[interval] <- ALT(HPA = data$HPA[interval], dt = settings$dt, halftime = simpopulation$halftime_pp[j], buffer_capacity = simpopulation$buffer_capacity_pp[j], r_recovery = simpopulation$r_recovery_pp[j])$a_load
}


# plot curves
for (j in seq(settings$N_people)){
  interval <- (1+(j-1)*settings$N):(j*settings$N)
  plot(cbind(
    ts(data$HPA[interval], start = settings$dt, end = 24*settings$days, frequency = 1/settings$dt),
    ts(data$cortisol[interval], start = settings$dt, end = 24*settings$days, frequency = 1/settings$dt)
  ),
  yax.flip = T, ann = F, col = "blue", frame.plot = T)
  title(main =  paste0("Cortisol profile person ",j),
        ylab = "|   Output (cortisol level)   |   Input (cortisol release)   |",
        xlab = "time (hours)")
}


# plot quantiles

df.qt <- data %>% group_by(time) %>% summarise(
  Q05=quantile(cortisol,probs=.05),
  Q25=quantile(cortisol,probs=.25),
  Q50=quantile(cortisol,probs=.50),
  Q75=quantile(cortisol,probs=.75),
  Q95=quantile(cortisol,probs=.95))

with(df.qt, ts.plot( 
  ts(Q05, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
  ts(Q25, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
  ts(Q50, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
  ts(Q75, start = 0, end = 24-settings$dt, frequency = 1/settings$dt),
  ts(Q95, start = 0, end = 24-settings$dt, frequency = 1/settings$dt)
))


# plot quantiles since awakening

#make a new dataframe that includes only the first day
cortisol.sinceawakening <- c(NULL)
time.sinceawakening <- c(NULL)
for (j in seq(settings$N_people)){
  #per person we want the data from the first day from awakening until end of the day
  #wake time for person one is data.day$wake_time[1], for person two this is data.day$wake_time[1+settings$days]
  interval <- which( (data$time > data.day$wake_time[1+(j-1)*settings$days]) & (data$sim_Person == data$sim_Person[j*settings$N])
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
#   hist(as.matrix(A_load[24/dt*i,]), main = paste0('Distribution of allostatic load for ', N_people, ' people at end of day ', i), xlab = 'Allostatic load' )
# }

#hist(cbind( as.matrix(Buffer[24/dt*i,]) , as.matrix(A_load[24/dt*i,])) )

## determine how long the script took
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# ## visualize for eacht person how their allostatic load increases
# 
# ggplot(data,aes(x=time, y=buffer, colour=sim_Person)) +
#   geom_line() +
#   theme(legend.position="none") + #remove legend
#   ggtitle(paste("Load buffer variation - N =",N_people)) +
#   ylab("available buffer")
# 
# ## visualize for eacht person how their allostatic load increases
# 
# ggplot(data,aes(x=time, y=a_load, colour=sim_Person)) +
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

#add allostatic load at noon each day to data.day
data.day$allostatic_load <- data$a_load[with(settings, seq(12/dt,N*N_people,24/dt))]

## plot stressor - allostatic load relationship

#choose subset of days to plot
plotdays <- seq( as.integer(settings$days/10), settings$days, as.integer(settings$days/10) )

ggplot(subset(data.day, day %in% plotdays), aes(x=acute_f_mean_pp, y=allostatic_load)) +
  geom_point() + facet_wrap( ~ day, ncol=5) + geom_smooth(method=lm, se=F) + 
  #theme(legend.position="none") + #remove legend
  ggtitle(paste("Allostatic load vs. work stress, N =",settings$N_people)) +
  ylab("Allostatic load") + xlab("average frequency of stressors")

