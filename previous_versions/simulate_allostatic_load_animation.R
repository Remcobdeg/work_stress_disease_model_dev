
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
  settings <- list(dt = NA) #initialize dataframe
  
  ## simulation settings
  settings$dt <- dt <- 1/30 #sampling time, in hours
  settings$N <- N <- 24/dt #number of samples in one day
  settings$days <- days <- 10
  settings$people <- people <- 1
  settings$integration <- "RK4"
  
  ## parameters that determine the cortisol profile
  
  #cortisol decay
  settings$halftimeM <- 80/60 #half time of cortisol is 80 minutes or 4/3hour
  settings$halftimeSD <- 8/60 #between person variation

  #awakening impulse
  settings$wakeDelay <- .5 #time after onset of the awakening impulse that people (report) waking up 
  settings$wakefM <- 40 
  settings$wakefSD <- 10
  settings$waketime <- 6.5

  #stress stressors
  settings$stressfM <- 10; 
  settings$stressfSD <- 5 #frequency of stress stressors
  settings$worktime <- list(8,17)

  ## parameters that determine how a cortisol profile translates into allostatic load and sickness
  settings$thLM <- 15 #threshold of dynamic load becoming allostatic load 
  settings$thLSD <- 0
  settings$rLM <- .3 #rate at which the body recovers form the dynamic load
  settings$rLSD <- 0
  settings$thAM <- 100 #threshold of allostatic load causing disease
  settings$thASD <- 0
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

#### 3.1 CALCULATE CORTISOL & ALLOSTATIC LOAD ####
CTA <- function(HPA_0, C_0, L_0, A_0, HPA, dt, 
                halftime, thL, rL, integration = c("Euler", "RK4")[1]){
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
  C <- as.numeric(c(C_0, rep(NA,N)))
  L <- as.numeric(c(L_0, rep(NA,N))) #initialize
  A <- as.numeric(c(A_0, rep(NA,N))) #initialize
  
  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           for (t in 1:N){
             
             C[t+1] <- Euler(Y1 = C[t], func = lingrowth, param = c(r, HPA[t]), h = dt) 
             
             L[t+1] <- Euler(Y1 = L[t], func = lingrowth, param = c(-rL, C[t]), h = dt) 
             #L[t+1] <- L[t] - L[t] * rL * dt + C[t] * dt 
             
             if (L[t] > thL){
               A[t+1] <- A[t] + (L[t] - thL)*dt}
             else {A[t+1] <- A[t]}
           },
         
         "RK4" =   # Runge-Kutta 4th Order Integration
         for (t in 1:N){
           
           C[t+1] <- RK4(Y1 = C[t], func = lingrowth, param = c(r, HPA[t]), h = dt) #RK4 numerical integration
           
           L[t+1] <- RK4(Y1 = L[t], func = lingrowth, param = c(-rL, C[t]), h = dt) 
           
           if (L[t] > thL){A[t+1] <- A[t] + (L[t] - thL)*dt}
           else {A[t+1] <- A[t]}
         }
         )

  return(data.frame(C[2:(N+1)], L[2:(N+1)], A[2:(N+1)])) #remove the start value 
  
}

#### 3.2 Euler integration ####
Euler <- function(Y1, func, param, h){
  ## This function performs one iteration using the Euler Integration

  # Iterative process
  Y2 <- Y1 + func(Y1, param)*h
  
  return(Y2)
}

#### 3.3 Runge-Kutta 4th Order Integration ####
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

#### 3.4 Linear growth (used in Euler and RK4 functions) ####
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

n <- 1 #seed counter; everytime a random number is drawn, the next random number is drawn using a new seed, but every execution of this script causes the same random numbers to be used
daycount <- 1 #is increased every time a new day is simulated

## Create a dataframe to store information for each person per day 
daydata <- data.frame( simperson = rep(as.character(paste0('SimPerson_',1:settings$people)), settings$days) )

#add day information
daydata$day <- with(settings, rep(1:days, each = people))

#add stressor frequency information
daydata$day <- with(settings, rep(NA, people*days))

#add holiday information
daydata$workday <- with(settings, rep(workday[1:days], each = people))
rm(workday)


#### 4.2 Create people from population characteristics ####

#create data frame to store all information about individuals
simpopulation <- data.frame(simperson = as.character(paste0('SimPerson_',1:settings$people)) )

## parameters that determine the cortisol profile

##Cortisol decay
set.seed(n); n=n+1; simpopulation$halftime_pp <- abs( rnorm(settings$people, mean = settings$halftimeM, sd = settings$halftimeSD) )

##awakening impulse
#indivual mean & within subject variation (same as 'intra-individual variation')
set.seed(n); n=n+1; simpopulation$wakefM_pp <- round( abs( rnorm(settings$people, mean = settings$wakefM, sd = settings$wakefSD) )) #can have only positive values

##stress stressors
#indivual mean
set.seed(n); n=n+1; simpopulation$stressfM_pp <- round( abs( rnorm(settings$people, mean = settings$stressfM, sd = settings$stressfSD) )) #can have only positive values

## parameters that determine how a cortisol profile translates into allostatic load and sickness
set.seed(n); n=n+1; simpopulation$thL_pp <- abs( rnorm(settings$people, mean = settings$thLM, sd = settings$thLSD) ) 
set.seed(n); n=n+1; simpopulation$rL_pp <- abs( rnorm(settings$people, mean = settings$rLM, sd = settings$rLSD) ) 
set.seed(n); n=n+1; simpopulation$thA_pp <- abs( rnorm(settings$people, mean = settings$thAM, sd = settings$thASD) ) 


#### 4.3 Prepare for simulation ####


# extend the data.frame daydata to store end day values of cortisol etc 
daydata$end_HPA <- with(settings, rep( c(0, rep(NA,days-1)), each = people))
daydata$end_C <- with(settings, rep( c(0, rep(NA,days-1)), each = people)) 
daydata$end_L <- with(settings, rep( c(0, rep(NA,days-1)), each = people)) 
daydata$end_A <- with(settings, rep( c(0, rep(NA,days-1)), each = people))


# # we want to plot curves of some 5 random people, first 2 days. Define which
# set.seed(n); n=n+1; 
# if (settings$people >= 5){
#   psample <- sample(1:settings$people, 5, replace = FALSE, prob = NULL) #5 random people
# } else {psample <- 1:settings$people}



#### 4.4 Simulate #### 

#per day
for (day in seq(days)){
  
  #per person 
  for (sim in seq(people)){
    
    #### 4.4.2 create HPA activity ####
  
    #create stress stress responses
    if (daydata$workday[day]==T){f <- round(simpopulation$stressfM_pp[sim])} else {f <- 0}
    set.seed(n); n=n+1; HPA_s <- runif(f,min = settings$worktime[[1]], max = settings$worktime[[2]])

    #create awakening inpuls
    wakef <- daydata$stressfM <- simpopulation$wakefM_pp
    set.seed(n); n=n+1; HPA_a <- abs(-rexp(wakef) + with(settings, waketime + wakeDelay))

    HPA <- 1:N
    for(i in 1:N){HPA[i]<-sum(round(c(HPA_s,HPA_a)/dt)==i)}
    
    #### 4.4.3 calculate cortisol, L and allostatic load from the HPA activity ####
    
    #intialize
    data <- data.frame(
      simperson=rep(simpopulation$simperson[sim],N),
      time=(1:N)*dt + 24*(day-1),
      HPA=HPA, #HPA_s + HPA_a + HPA_c,
      C=as.numeric(c(rep(NA,N))),
      L=as.numeric(c(rep(NA,N))),
      A=as.numeric(c(rep(NA,N))),
      stringsAsFactors=FALSE)
  
    ## calculate cortisol and allostatic load
    
    #asign input for the function
    dt <- dt 
    halftime <- simpopulation$halftime_pp[sim]
    thL <- simpopulation$thL_pp[sim]
    rL <- simpopulation$rL_pp[sim]
    integration <- settings$integration
    
    if (day == 1){
      HPA_0 <- 0; 
      C_0 <- 0; 
      L_0 <- 0; 
      A_0 <- 0
      
      data$C <- CTA(HPA_0,C_0,L_0,A_0,HPA,dt,halftime,thL,rL,integration
      )$C
      data$L <-  CTA(HPA_0,C_0,L_0,A_0,HPA,dt,halftime,thL,rL,integration
      )$L
      data$A <-    CTA(HPA_0,C_0,L_0,A_0,HPA,dt,halftime,thL,rL,integration
      )$A
      
    } else {
      HPA_0 <- daydata$end_HPA[daycount-settings$people]; 
      C_0 <- daydata$end_C[daycount-settings$people];
      L_0 <- daydata$end_L[daycount-settings$people];
      A_0 <- daydata$end_A[daycount-settings$people]
      
      data$C <- CTA(HPA_0,C_0,L_0,A_0,HPA,dt,halftime,thL,rL,integration
      )$C
      data$L <-  CTA(HPA_0,C_0,L_0,A_0,HPA,dt,halftime,thL,rL,integration
      )$L
      data$A <-    CTA(HPA_0,C_0,L_0,A_0,HPA,dt,halftime,thL,rL,integration
      )$A
    }
    
    #### 4.4.4 Store results ####
    
    #day results
    daydata$end_HPA[daycount] <- data$HPA[N]    
    daydata$end_C[daycount] <- data$C[N]
    daydata$end_L[daycount] <- data$L[N]
    daydata$end_A[daycount] <- data$A[N]
    daydata$sick[daycount] <- data$A[N]>simpopulation$thA_pp[sim]
    
    ## store data of the first day for creating the quantile plot of cortisol profiles as in Miller et al.
    # if (day == 1){
    #   if (exists("Qdata")){
    #     Qdata <- rbind(Qdata,data.frame(time = data$time, C = data$C))
    #   } else {Qdata <- data.frame(time = data$time, C = data$C)}}
    # 
    ## store the complete data of a few people for creating example plots of fluctuations of cortisol and allotic load within a day 
    # if (sim %in% psample & day <= 2){
    #   if (exists("plotdata")){
    #     plotdata <- rbind(plotdata,data)
    #   } else {plotdata <- data}}
    
    #animation data
    if (sim ==1 & day <= 10){
      if (exists("animdata")){
        animdata <- rbind(animdata,data)
      } else {animdata <- data}}
    
    daycount <- daycount + 1
  } # completes the day simulation
    
} # completes the person simulation

end.time <- Sys.time()
time.sim <- end.time - start.time
time.sim


################################### 5 MAKE CURVES ###################################

# if (file.exists(mainDir)){setwd(mainDir)}
# library(rmarkdown)
# render("1-example.Rmd")

# example graph - plot the parts constituting cortisol production (HPA activity)  
HPAs <- 1:N; for(i in 1:N){HPAs[i]<-sum(round(HPA_s/dt)==i)}
HPAa <- 1:N; for(i in 1:N){HPAs[i]<-sum(round(HPA_a/dt)==i)}
if (sim == 1){plot(cbind(
    ts(HPA), ts(HPAs), ts(HPAa)
))}

#### 5.6 Animation data ####

#HPA output
ggplot(animdata, aes(x=time, y=HPA)) +
  geom_line() + 
  ggtitle("Cortisol production (= HPA-axis activity) of 5 randomly selected people (first 2 days)") +
  ylab("Cortisol production [mmol/L]") + xlab("time [h]")



#cortisol profile
ggplot(animdata, aes(x=time, y=C)) +
  geom_line() + 
  ggtitle("Bloodstream cortisol of 5 randomly selected people (first 2 days)") +
  ylab("Cortisol [mmol/L]") + xlab("time [h]")

#dynamic load profile
ggplot(animdata, aes(x=time, y=L)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Dynamic load of 5 randomly selected people (first 2 days)") +
  ylab("Dynamic load") + xlab("time [h]")

#allostatic load profile
ggplot(animdata, aes(x=time, y=A)) +
  geom_line() + facet_wrap( ~ simperson, ncol=1) +
  ggtitle("Allostatic load of 5 randomly selected people (first 2 days)") +
  ylab("Allostatic load") + xlab("time [h]")

#package("devtools")
#devtools::install_github("dgrtwo/gganimate")

animdata1 <- animdata[(1):(30),] 
#for (i in 1:(nrow(animdata)-30)){animdata1 <- rbind(animdata1,animdata[(i+1):(i+30),])}
for (i in 1:(50-30)){animdata1 <- rbind(animdata1,animdata[(i+1):(i+30),])}
animdata1$snap <- rep(1:(i+1), each = 30)

p <- ggplot(animdata1, aes(time, HPA, frame = snap)) + geom_line() 

library(gganimate)

gganimate(p, interval = .2)




################################### 6 Some finishing action ###################################

if (save_settings == T){
  
  ## create and set saving directory
  mainDir <- getwd()
  subDir <- paste("simulation", format(Sys.time(),"%Y-%m-%d_%H-%M"))
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  write.csv(settings, "settings.csv")
  write.csv(simpopulation, "simpopulation.csv")
  write.csv(daydata, "daydata.csv")
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
# create a (data driven) means to determine the 'right' load recovery, load thA, disease threshold






