
##############################################
# This script ..............................
# # ----CREATING CORTISOL PROFILE BY SIMULATING SUDDEN CORTISOL PRODUCTION BURST USING SIMPLE FORMULA----
# # The circardian cortisol profile will be constructed as a periodically varying cortisol production of the HPA-axis,
# minus decay of cortisol using:
# ..............................
##############################################

rm(list = ls()) #clear the environment

load_from_file <- F 
save_settings <- F

winDir <- "D:/Surfdrive/BS28A - Major project/simulations/simulations"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simulations"
if (file.exists(winDir)){mainDir <- winDir} else if (file.exists(macDir)){mainDir <- macDir} else {mainDir <- choose.dir()}
setwd(mainDir)


######## 1 DEFINE PARAMETERS AND SIMULATION SETTINGS ########

#### 1.1 define work vs. non-work days per person  ####
#Done seperate at it is (not include in the loaded file)
workday <- rep( c( rep( c(rep(T,5),rep(F,2)),12 ), rep(F,7) ), 20) #12 working weeks with 2 days weekend, then a week of holidays, repeated 20 times

#### 1.2 Load parameters from file ####

if (load_from_file == T){
  settings <- read.csv(file.choose())}

#### 1.3 Set parameters manually ####

if (load_from_file == F){

  #settings <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL))
  settings <- data.frame(dt = NA) #initialize dataframe
  
  ## simulation settings
  settings$dt <- 1/2 #sampling time, in hours
  settings$days <- 100
  settings$N_people <- 1000
  settings$integration <- "RK4"
  
  ## parameters that determine the cortisol profile
  
  #cortisol decay
  settings$halftime_mean <- 80/60 #half time of cortisol is 80 minutes or 4/3hour
  settings$halftime_sd <- 8/60

  #circadian pulse
  settings$circ_peaktime_mean <- 6
  settings$circ_peaktime_sd <- 1
  settings$circ_intensity_mean <- 0 #1
  settings$circ_intensity_sd <- 0 #.2
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
  
  ## parameters that determine how a cortisol profile translates into allostatic load and sickness
  settings$buffer_capacity_mean <- 15
  settings$buffer_capacity_sd <- 0
  settings$r_recovery_mean <- .3 #rate at which the buffer restores
  settings$r_recovery_sd <- 0
  settings$tolerance_mean <- 1000
  settings$tolerance_sd <- 0
}

#### 2 LOAD PACKAGES #####

#function to load a library and install it when it is not yet installed
package <- function(pack) { #function to load a library and install it when it is not yet installed
  if (!pack %in% installed.packages()) {install.packages(pack) }
  library(pack,character.only=TRUE)}

package("lattice")
package("plyr")
package("dplyr")
package("ggplot2")
package("car")


##### 3 DEFINE FUNCTIONS #####


#### 3.1 BLOCK SHAPED BURST ####
blockburst <- function(intensity = 1, duration = 2, start = 6, N = 1440, dt = 1/60){
  Y <- rep(0,N)
  Y[ round(start/dt) : round((start+duration)/dt)] <- intensity
  return(Y)
}

#### 3.2 SINUS SHAPED CIRCADIAN FLUCTUATION ####
circ <- function(intensity, peak_time, N){
  Y <- rep(NA,N) #initialize
  for (t in 1:N){
    Y[t] <- sin( 2*pi * (t/N - (peak_time - 6)/24) ) +1 # create the sinus shape (fluctuates from 0 to 2)
  }
  Y <- Y * intensity #scale the signal to the requested intensity
  return(Y)
}


#### 3.3 CALCULATE CORTISOL & ALLOSTATIC LOAD ####
CTA <- function(HPA_0, cortisol_0, tension_0, Aload_0, HPA, dt, 
                halftime, buffer_capacity, r_recovery, integration = c("Euler", "RK4")[1]){
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
  
  #initiate and add values start values (last values of the previous day)
  N <- length(HPA)
  HPA <- c(HPA_0, HPA)
  cortisol <- as.numeric(c(cortisol_0, rep(NA,N)))
  tension <- as.numeric(c(tension_0, rep(NA,N))) #initialize
  Aload <- as.numeric(c(Aload_0, rep(NA,N))) #initialize
  
  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           for (t in 1:N){
             
             cortisol[t+1] <- Euler(Y1 = cortisol[t], func = lingrowth, param = c(r, HPA[t]), h = dt) 
             
             tension[t+1] <- Euler(Y1 = tension[t], func = lingrowth, param = c(-r_recovery, cortisol[t]), h = dt) 
             #tension[t+1] <- tension[t] - tension[t] * r_recovery * dt + cortisol[t] * dt 
             
             if (tension[t] > buffer_capacity){
               Aload[t+1] <- Aload[t] + (tension[t] - buffer_capacity)*dt}
             else {Aload[t+1] <- Aload[t]}
           },
         
         "RK4" =   # Runge-Kutta 4th Order Integration
         for (t in 1:N){
           
           cortisol[t+1] <- RK4(Y1 = cortisol[t], func = lingrowth, param = c(r, HPA[t]), h = dt) #RK4 numerical integration
           
           tension[t+1] <- RK4(Y1 = tension[t], func = lingrowth, param = c(-r_recovery, cortisol[t]), h = dt) 
           
           if (tension[t] > buffer_capacity){Aload[t+1] <- Aload[t] + (tension[t] - buffer_capacity)*dt}
           else {Aload[t+1] <- Aload[t]}
         }
         )

  return(data.frame(cortisol[2:(N+1)], tension[2:(N+1)], Aload[2:(N+1)])) #remove the start value 
  
}

#### 3.4 Euler integration ####
Euler <- function(Y1, func, param, h){
  ## This function performs one iteration using the Euler Integration

  # Iterative process
  Y2 <- Y1 + func(Y1, param)*h
  
  return(Y2)
}

#### 3.5 Runge-Kutta 4th Order Integration ####
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

#### 3.6 Linear growth (used in Euler and RK4 functions) ####
#This function performs one iteration of linear growth
lingrowth <- function(Y1,param){
  r <- param[1]; input <- param[2]
  Y2 <- Y1*r + input
  return(Y2)
}


################################### 4 EXECUTION PART ###################################

#### 4.1 Initialize ####

## time stamp the start of the script execution, so we can determine how long the whole thing took
start.time <- Sys.time()

N <- 24/settings$dt #number of samples in one day; defined for convenience purposes as this is often used

n <- 1 #seed counter; everytime a random number is drawn, the next random number is drawn using a new seed, but every execution of this script causes the same random numbers to be used
daycount <- 1 #is increased every time a new day is simulated


## Create a dataframe to store information for each person per day 
data.day <- data.frame( simperson = rep(as.character(paste0('SimPerson_',1:settings$N_people)), settings$days) )

#add day information
data.day$day <- with(settings, rep(1:days, each = N_people))

#add holiday information
data.day$workday <- with(settings, rep(workday[1:days], each = N_people))
rm(workday)


#### 4.2 Create people from population characteristics ####

#create data frame to store all information about individuals
simpopulation <- data.frame(simperson = as.character(paste0('SimPerson_',1:settings$N_people)) )

## parameters that determine the cortisol profile

##Cortisol decay
set.seed(n); n=n+1; simpopulation$halftime_pp <- abs( rnorm(settings$N_people, mean = settings$halftime_mean, sd = settings$halftime_sd) )

##awakening impulse
#indivual mean & within subject variation (same as 'intra-individual variation')
set.seed(n); n=n+1; simpopulation$wake_delay_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_delay_mean, sd = settings$wake_delay_sd) ) #can have only positive values
set.seed(n); n=n+1; simpopulation$wake_delay_sd_pp <- simpopulation$wake_delay_mean_pp * settings$multiplier_wake_delay_sd_pp 
set.seed(n); n=n+1; simpopulation$wake_intensity_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_intensity_mean, sd = settings$wake_intensity_sd) ) #can have only positive values
set.seed(n); n=n+1; simpopulation$wake_intensity_sd_pp <- simpopulation$wake_intensity_mean_pp * settings$multiplier_wake_intensity_sd_pp 
set.seed(n); n=n+1; simpopulation$wake_duration_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_duration_mean, sd = settings$wake_duration_sd) ) #can have only positive values
set.seed(n); n=n+1; simpopulation$wake_duration_sd_pp <- simpopulation$wake_duration_mean_pp * settings$multiplier_wake_duration_sd_pp 
set.seed(n); n=n+1; simpopulation$wake_time_mean_pp <- abs( rnorm(settings$N_people, mean = settings$wake_time_mean, sd = settings$wake_time_sd) ) #can have only positive values
set.seed(n); n=n+1; simpopulation$wake_time_sd_pp <- settings$wake_time_sd_pp #already defined above, without using a multiplier

##circadian pulse
set.seed(n); n=n+1; simpopulation$circ_peaktime_pp <-rnorm(settings$N_people, mean = settings$circ_peaktime_mean, sd = settings$circ_peaktime_sd)  
set.seed(n); n=n+1; simpopulation$circ_intensity_mean_pp <- abs( rnorm(settings$N_people, mean = settings$circ_intensity_mean, sd = settings$circ_intensity_sd) )  
set.seed(n); n=n+1; simpopulation$circ_intensity_sd_pp <- simpopulation$circ_intensity_mean_pp * settings$multiplier_circ_intensity_sd_pp 

##acute stressors
#indivual mean
set.seed(n); n=n+1; simpopulation$acute_f_mean_pp <- abs( rnorm(settings$N_people, mean = settings$acute_f_mean, sd = settings$acute_f_sd) ) #can have only positive values
set.seed(n); n=n+1; simpopulation$acute_f_sd_pp <- simpopulation$acute_f_mean_pp * settings$multiplier_acute_f_sd_pp
set.seed(n); n=n+1; simpopulation$acute_intensity_mean_pp <-  abs( rnorm(settings$N_people, mean = settings$acute_intensity_mean, sd = settings$acute_intensity_sd) ) #can have only positive values
set.seed(n); n=n+1; simpopulation$acute_intensity_sd_pp <- simpopulation$acute_intensity_mean_pp * settings$multiplier_acute_intensity_sd_pp
set.seed(n); n=n+1; simpopulation$acute_duration_mean_pp <-  abs( rnorm(settings$N_people, mean = settings$acute_duration_mean, sd = settings$acute_duration_sd) ) #can have only positive values
set.seed(n); n=n+1; simpopulation$acute_duration_sd_pp <- simpopulation$acute_duration_mean_pp * settings$multiplier_acute_duration_sd_pp
#within subject variation (same as 'intra-individual variation')

## parameters that determine how a cortisol profile translates into allostatic load and sickness
set.seed(n); n=n+1; simpopulation$buffer_capacity_pp <- abs( rnorm(settings$N_people, mean = settings$buffer_capacity_mean, sd = settings$buffer_capacity_sd) ) 
set.seed(n); n=n+1; simpopulation$r_recovery_pp <- abs( rnorm(settings$N_people, mean = settings$r_recovery_mean, sd = settings$r_recovery_sd) ) 
set.seed(n); n=n+1; simpopulation$tolerance_pp <- abs( rnorm(settings$N_people, mean = settings$tolerance_mean, sd = settings$tolerance_sd) ) 


#### 4.3 Prepare for simulation ####

# first extend the data.frame data.day to capture daily values of stressors 
data.day$wake_time <- with(settings, rep(NA, days*N_people))
if (settings$multiplier_wake_delay_sd_pp != 0){data.day$wake_delay <- with(settings, rep(NA, days*N_people))}
if (settings$multiplier_wake_intensity_sd_pp != 0){data.day$wake_intensity <- with(settings, rep(NA, days*N_people))}
if (settings$multiplier_wake_duration_sd_pp != 0){data.day$wake_duration <- with(settings, rep(NA, days*N_people))}
if (settings$multiplier_circ_intensity_sd_pp != 0){data.day$circ_intensity <- with(settings, rep(NA, days*N_people))}
if (settings$multiplier_acute_f_sd_pp != 0){data.day$acute_f <- with(settings, rep(NA, days*N_people))}
# note that duration and intensity of each individual stressor cannot be stored in this data frame, as there are multiple values per day
if ( (settings$multiplier_acute_intensity_sd_pp != 0) | (settings$multiplier_acute_duration_sd_pp != 0) ){
  data.stressor <- data.frame(simperson = 'Dummy', day = 0, duration = 0, intensity = 0)
}

# extend the data.frame data.day to store end day values of cortisol etc 
data.day$end_HPA <- with(settings, rep( c(0, rep(NA,days-1)), each = N_people))
data.day$end_cortisol <- with(settings, rep( c(0, rep(NA,days-1)), each = N_people)) 
data.day$end_tension <- with(settings, rep( c(0, rep(NA,days-1)), each = N_people)) 
data.day$end_Aload <- with(settings, rep( c(0, rep(NA,days-1)), each = N_people))


# we want to plot curves of some 5 random people, first 2 days. Define which
set.seed(n); n=n+1; 
if (settings$N_people >= 5){
  psample <- sample(1:settings$N_people, 5, replace = FALSE, prob = NULL) #5 random people
} else {psample <- 1:settings$N_people}



#### 4.4 Simulate #### 

#per day
for (day in seq(settings$days)){
  
  #per person 
  for (sim in seq(settings$N_people)){
    
    #### 4.4.1 define components of the cortisol producing impulses ####
  
    #### 4.4.1.1 of acute stressors ####
    
    #frequency
    if (settings$multiplier_acute_f_sd_pp != 0){ 
      if (data.day$workday[day]==T){
        set.seed(n);  
        f <- round( abs( rnorm(1, mean = simpopulation$acute_f_mean_pp[sim], 
                               sd = simpopulation$acute_f_sd_pp[sim]) )) #can have only positive values
      } else {f <- 0}
      data.day$acute_f[daycount] <- f
    } else {
      if (data.day$workday[day]==T){f <- round(simpopulation$acute_f_mean_pp[sim])} else {f <- 0}}
    n=n+1;
    
    #intensity
    if (settings$multiplier_acute_intensity_sd_pp != 0){ 
      set.seed(n); 
      acute_intensity <- abs(rnorm(f, mean = simpopulation$acute_intensity_mean_pp[sim], 
                                   sd = simpopulation$acute_intensity_sd_pp[sim]))
    } else {acute_intensity <- rep(simpopulation$acute_intensity_mean_pp, f)}
    n=n+1; 
    
    #duration
    if (settings$multiplier_acute_duration_sd_pp != 0){ 
      set.seed(n);
      acute_duration <- abs(rnorm(f, mean = simpopulation$acute_duration_mean_pp[sim], 
                                   sd = simpopulation$acute_duration_sd_pp[sim]))
    } else {acute_duration <- rep(simpopulation$acute_duration_mean_pp, f)}
    n=n+1; 
    
    #start time
    set.seed(n); 
    acute_start <- runif(f,settings$work_start,settings$work_end)
    n=n+1; 
    
    #### 4.4.1.1 of the awakening impulse ####
    
    #intensity
    if (settings$multiplier_wake_intensity_sd_pp != 0){ 
      set.seed(n); 
      wake_intensity <- abs(rnorm(1, mean = simpopulation$wake_intensity_mean_pp[sim], 
                                   sd = simpopulation$wake_intensity_sd_pp[sim]))
      data.day$wake_intensity[daycount] <- wake_intensity
    } else {wake_intensity <- simpopulation$wake_intensity_mean_pp[sim]}
    n=n+1; 
    
    #duration
    if (settings$multiplier_wake_duration_sd_pp != 0){ 
      set.seed(n); 
      wake_duration <- abs(rnorm(1, mean = simpopulation$wake_duration_mean_pp[sim], 
                                  sd = simpopulation$wake_duration_sd_pp[sim]))
      data.day$wake_duration[daycount] <- wake_duration
    } else {wake_duration <- simpopulation$wake_duration_mean_pp[sim]}
    n=n+1; 
    
    #wake time
    if (settings$wake_time_sd_pp != 0){ 
      set.seed(n); 
      wake_time <- abs(rnorm(1, mean = simpopulation$wake_time_mean_pp[sim], 
                                 sd = simpopulation$wake_time_sd_pp[sim]))
    } else {wake_time <- simpopulation$wake_time_mean_pp[sim]}
    data.day$wake_time[daycount] <- wake_time
    n=n+1; 
    
    #start of the awakening impulse
    if (settings$multiplier_wake_delay_sd_pp != 0){ 
      set.seed(n); 
      wake_start <- abs(wake_time - rnorm(1, mean = simpopulation$wake_delay_mean_pp[sim], 
                             sd = simpopulation$wake_delay_sd_pp[sim]))
      data.day$wake_start[daycount] <- wake_start
    } else {wake_start <- wake_time - simpopulation$wake_delay_mean_pp[sim]}
    n=n+1; 

    #### 4.4.1.3 of the circadian rhythm ####
    
    #intensity
    if (settings$multiplier_circ_intensity_sd_pp != 0){ 
      set.seed(n); 
      circ_intensity <- abs(rnorm(1, mean = simpopulation$circ_intensity_mean_pp[sim], 
                                  sd = simpopulation$circ_intensity_sd_pp[sim]))
      data.day$circ_intensity[daycount] <- circ_intensity
    } else {circ_intensity <- simpopulation$circ_intensity_mean_pp[sim]}
    n=n+1; 

    
    #### 4.4.2 create HPA activity ####
  
    #create acute stress responses
    HPA_s <- c() #initialize

    if (f>0){
      HPA_s <- matrix(NA,N,f) #initialize a column for each stressor, so they can be added later
      #create each individual stressor
      for (i in seq(f)){
        HPA_s[,i] <- blockburst(acute_intensity[i], acute_duration[i], acute_start[i], N, settings$dt)
      }
      HPA_s <- apply(HPA_s,1,max) #if stress impulses overlap, the strongest impulse is used. They are not summed
    } else{HPA_s <- rep(0,N)} #if f <= 0 it produces a row of zeros

    #create awakening inpuls
    HPA_a <- blockburst(wake_intensity, wake_duration, wake_start, N, settings$dt)
  
    #create a circadian sinus fluctuation
    HPA_c <- circ(circ_intensity, simpopulation$circ_peaktime_pp[sim], N)
    

    #### 4.4.3 calculate cortisol, tension and allostatic load from the HPA activity ####
    
    #intialize
    data <- data.frame(
      simperson=rep(simpopulation$simperson[sim],N),
      time=(1:N)*settings$dt + 24*(day-1),
      HPA=HPA_s + HPA_a + HPA_c,
      cortisol=as.numeric(c(rep(NA,N))),
      tension=as.numeric(c(rep(NA,N))),
      Aload=as.numeric(c(rep(NA,N))),
      stringsAsFactors=FALSE)
  

    ## calculate cortisol and allostatic load
    
    #asign input for the function
    dt <- settings$dt 
    halftime <- simpopulation$halftime_pp[sim]
    buffer_capacity <- simpopulation$buffer_capacity_pp[sim]
    r_recovery <- simpopulation$r_recovery_pp[sim]
    integration <- settings$integration
    
    if (day == 1){
      HPA_0 <- 0; 
      cortisol_0 <- 0; 
      tension_0 <- 0; 
      Aload_0 <- 0
      
      data$cortisol <- CTA(HPA_0,cortisol_0,tension_0,Aload_0,data$HPA,dt,halftime,buffer_capacity,r_recovery,integration
      )$cortisol
      data$tension <-  CTA(HPA_0,cortisol_0,tension_0,Aload_0,data$HPA,dt,halftime,buffer_capacity,r_recovery,integration
      )$tension
      data$Aload <-    CTA(HPA_0,cortisol_0,tension_0,Aload_0,data$HPA,dt,halftime,buffer_capacity,r_recovery,integration
      )$Aload
      
    } else {
      HPA_0 <- data.day$end_HPA[daycount-settings$N_people]; 
      cortisol_0 <- data.day$end_cortisol[daycount-settings$N_people];
      tension_0 <- data.day$end_tension[daycount-settings$N_people];
      Aload_0 <- data.day$end_Aload[daycount-settings$N_people]
      
      data$cortisol <- CTA(HPA_0,cortisol_0,tension_0,Aload_0,data$HPA,dt,halftime,buffer_capacity,r_recovery,integration
      )$cortisol
      data$tension <-  CTA(HPA_0,cortisol_0,tension_0,Aload_0,data$HPA,dt,halftime,buffer_capacity,r_recovery,integration
      )$tension
      data$Aload <-    CTA(HPA_0,cortisol_0,tension_0,Aload_0,data$HPA,dt,halftime,buffer_capacity,r_recovery,integration
      )$Aload
    }
    
    #### 4.4.4 Store results ####
    
    #day results
    data.day$end_HPA[daycount] <- data$HPA[N]    
    data.day$end_cortisol[daycount] <- data$cortisol[N]
    data.day$end_tension[daycount] <- data$tension[N]
    data.day$end_Aload[daycount] <- data$Aload[N]
    data.day$sick[daycount] <- data$Aload[N]>simpopulation$tolerance_pp[sim]
    
    ## store data of the first day for creating the quantile plot of cortisol profiles as in Miller et al.
    if (day == 1){
      if (exists("Qdata")){
        Qdata <- rbind(Qdata,data.frame( time = data$time[which(data$time>wake_time)]-wake_time, cortisol = data$cortisol[which(data$time>wake_time)]))
      } else {Qdata <- data.frame(time = data$time[which(data$time>wake_time)]-wake_time, cortisol = data$cortisol[which(data$time>wake_time)])}}
    
    ## store the complete data of a few people for creating example plots of fluctuations of cortisol and allotic load within a day 
    if (sim %in% psample & day <= 2){
      if (exists("plotdata")){
        plotdata <- rbind(plotdata,data)
      } else {plotdata <- data}}
    
    
    daycount <- daycount + 1
  } # completes the day simulation
    
} # completes the person simulation

end.time <- Sys.time()
time.sim <- end.time - start.time
time.sim


################################### 5 MAKE CURVES ###################################

# example graph - plot the parts constituting cortisol production (HPA activity)  
if (sim == 1){plot(cbind(
    ts(HPA_s+HPA_a+HPA_c), ts(HPA_s), ts(HPA_a), ts(HPA_c)
))}

#### 5.1 Example day profiles of e.g. cortisol dynamics ####

#HPA output
ggplot(plotdata, aes(x=time, y=HPA)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Cortisol production (= HPA-axis activity) of 5 randomly selected people (first 2 days)") +
  ylab("Cortisol production [mmol/L]") + xlab("time [h]")

#cortisol profile
ggplot(plotdata, aes(x=time, y=cortisol)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Bloodstream cortisol of 5 randomly selected people (first 2 days)") +
  ylab("Cortisol [mmol/L]") + xlab("time [h]")

#dynamic load profile
ggplot(plotdata, aes(x=time, y=tension)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Dynamic load of 5 randomly selected people (first 2 days)") +
  ylab("Dynamic load") + xlab("time [h]")

#allostatic load profile
ggplot(plotdata, aes(x=time, y=Aload)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Allostatic load of 5 randomly selected people (first 2 days)") +
  ylab("Allostatic load") + xlab("time [h]")

#### 5.2 Average cortisol day profile since awakening #### 
#create plot of quantiles since awakening, as in Miller et al.

df.qt <- Qdata %>% group_by(time) %>% summarise(
  Q05=quantile(cortisol,probs=.05),
  Q25=quantile(cortisol,probs=.25),
  Q50=quantile(cortisol,probs=.50),
  Q75=quantile(cortisol,probs=.75),
  Q95=quantile(cortisol,probs=.95))

with(df.qt, ts.plot(
  ts(Q05, start = 0, end = (length(Q05)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q25, start = 0, end = (length(Q25)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q50, start = 0, end = (length(Q50)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q75, start = 0, end = (length(Q75)-1)*settings$dt, frequency = 1/settings$dt),
  ts(Q95, start = 0, end = (length(Q95)-1)*settings$dt, frequency = 1/settings$dt)
))

#### 5.3 Allostatic load ~ average stressor frequency ####

#choose subset of days to plot
if(settings$days<=10){plotdays <- 1:settings$days
} else {plotdays <- seq( as.integer(settings$days/10), settings$days, as.integer(settings$days/10) )}

#define average number of stressors per person
if("acute_f" %in% colnames(data.day)){
  data.day$acute_f_mean <- rep(NA, nrow(data.day))
  for (i in 1:nrow(data.day)){
    data.day$acute_f_mean[i] <- with(data.day, mean(acute_f[which(simperson==simperson[i] & workday == T)]))
    }
} else {
  data.day$acute_f_mean <- rep(round(simpopulation$acute_f_mean_pp),times = settings$days)
}

ggplot(subset(data.day, day %in% plotdays), aes(x=acute_f_mean, y=end_Aload)) +
  geom_point() + facet_wrap( ~ day, ncol=5) + geom_smooth(method="auto", se=F) +
  #theme(legend.position="none") + #remove legend
  ggtitle(paste("Allostatic load vs. work stress, N =",settings$N_people)) +
  ylab("Allostatic load") + xlab("average frequency of stressors")

#### 5.4 Odd of disease per additional stressor ####

#choose subset of days to plot
if(settings$days<=10){plotdays <- 1:settings$days
} else {plotdays <- seq( as.integer(settings$days/10), settings$days, as.integer(settings$days/10) )}

#make a dataframe to contain the odds ratios and confidence intervals for each day that the will plot
odds <- data.frame(day = as.ordered(plotdays), OR = rep(NA,length(plotdays)), lower = rep(NA,length(plotdays)), upper = rep(NA,length(plotdays)))

for (aday in plotdays){
  odds.data <- data.day[which(data.day$day == aday),]
  
  if (length(which(odds.data$sick==T))==0 | length(which(odds.data$sick==F))==0){
    odds$OR[which(odds$day == aday)] <- 1
    odds$lower[which(odds$day == aday)] <- 1
    odds$upper[which(odds$day == aday)] <- 1
  } else {
    #compute odds ratios
    odds_model <- with(odds.data, glm(sick~acute_f_mean, family = binomial()))
    
    # #assess assumptions
    # dwt(odds_model) #Independence of errors (Durbin-Watson)
    # densityplot(rstandard(odds_model)) #check for outliers
    # hist(dfbeta(odds_model)) #influential cases
    # densityplot(hatvalues(odds_model)) #influential cases
    # plot(odds_model)
    # 
    # #look at model
    # summary(odds_model)
    # 
    # #model fit by Chi2
    # Chi_sq <- odds_model$null.deviance-odds_model$deviance
    # chidf <- odds_model$df.null - odds_model$df.residual
    # chisq.prob <- 1 - pchisq(Chi_sq, chidf)
    # Chi_sq; chidf; chisq.prob
    # 
    # #model fit by R2
    # R2.hl<-Chi_sq/odds_model$null.deviance
    # R.cs<-1-exp((odds_model$deviance-odds_model$null.deviance)/nrow(odds.data))
    # R.n<-R.cs/(1-(exp(-(odds_model$null.deviance/nrow(odds.data)))))
    # R2.hl; R.cs; R.n
    
    #Store odds ratios
    odds$OR[which(odds$day == aday)] <- exp(odds_model$coefficients)[2] #position [1] is the intercept...
    odds$lower[which(odds$day == aday)] <- exp(confint(odds_model))[2,1]
    odds$upper[which(odds$day == aday)] <- exp(confint(odds_model))[2,2]
  }
}

#plot odds ratio's
ggplot(odds, aes(y= OR, x = day)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  geom_hline(yintercept = 1, linetype=2) +
  coord_flip() +
  labs(title = "OR per day", x = "Day", y = "OR") +
  theme_bw()

end.time <- Sys.time()
time.plotting <- end.time - start.time
time.plotting

################################### 6 Some finishing action ###################################

if (save_settings == T){
  
  ## create and set saving directory
  mainDir <- getwd()
  subDir <- paste("simulation", format(Sys.time(),"%Y-%m-%d_%H-%M"))
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  write.csv(settings, "settings.csv")
  write.csv(simpopulation, "simpopulation.csv")
  write.csv(data.day, "data.day.csv")
  write.csv(plotdata, "plotdata.csv")
  write.csv(Qdata, "Qdata.csv")
}

## determine how long the script took
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
time.sim
time.plotting


#20.000 people for 10 days took ~30min

## to improve:
# present odds ratios such that they are comparable to Kivimaki. Use standardized frequency? How is lifetime risk determined; only people that have died? How to define 'experiencing job strain' in my model? Based on an absolute frequency? Which reflects the number of people experiencing job strain? (determine this from ESS?)
# allow calculations with large numbers of people (20'000) for ca. 1000 days
# create a (data driven) means to determine the 'right' load recovery, load tolerance, disease threshold






