##############################################
# This script ..............................
# # ----CREATING CORTISOL PROFILE BY SIMULATING SUDDEN CORTISOL PRODUCTION BURST USING SIMPLE FORMULA----
# # The circardian cortisol profile will be constructed as a periodically varying cortisol production of the HPA-axis,
# minus decay of cortisol using:
# ..............................
##############################################

rm(list = ls()) #clear the environment

#### 1 LOAD PACKAGES & SET DIRECTORY #####

#function to load a library and install it when it is not yet installed
package <- function(pack) { #function to load a library and install it when it is not yet installed
  if (!pack %in% installed.packages()) {install.packages(pack) }
  library(pack,character.only=TRUE)}

package("lattice")
package("plyr")
package("dplyr")
package("ggplot2")
package("tcltk")
package("car")

## set saving directories

#checks if a dedicated data output folder is available in the working directory and, if not, ask the user to define the directory
winDir <- "D:/Surfdrive/BS28A - Major project/simulations"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 
if (!file.exists("simoutput")){setwd(tk_choose.dir(getwd(), "Choose folder for storing simulation data"))
  if (!file.exists("simoutput")){dir.create("simoutput")}} 
maindir <- getwd()
setwd(file.path(getwd(),"simoutput"))


######## 2 DEFINE PARAMETERS AND SIMULATION SETTINGS ########

#### 2.1 Load parameters from file ####

load_from_file <- F 

if (load_from_file == T){
  settings <- read.csv(file.choose())}

#### 2.2 Set parameters manually ####

if (load_from_file == F){
  
  settings <- list(dt = NA) #initialize list
  
  ## simulation settings
  settings$fs <- fs <- 30 #samples per hour
  settings$N <- 24*fs #number of samples in one day
  settings$days <- 10
  settings$people <- 1000
  settings$integration <- "RK4"
  
  ## off-work days
  settings$dayoff <- c(0,6,7) #weekday 0 == weekday 7; for computational reasons, always include both
  
  ## parameters that determine the cortisol profile
  
  #cortisol decay
  settings$halftimeM <- 80/60 #half time of cortisol is 80 minutes or 4/3hour
  settings$halftimeSD <- 8/60 #between person variation
  
  #awakening impulse
  settings$CAR <- .5 #time after onset of the awakening impulse that people (report) waking up 
  settings$wakefM <- 10 
  settings$wakefSD <- 3
  settings$waketime <- 6.5
  
  #work stressors
  settings$stressfM <- 10; 
  settings$stressfSD <- 5 #frequency of stress stressors
  settings$worktime <- c(8,17)
  
  ## parameters that determine how a cortisol profile translates into allostatic load and sickness
  settings$thLM <- 15 #threshold of dynamic load becoming allostatic load 
  settings$thLSD <- 0
  settings$rLM <- .3 #rate at which the body recovers form the dynamic load
  settings$rLSD <- 0
  settings$thAM <- 100 #threshold of allostatic load causing disease
  settings$thASD <- 0
  
  ## number of days to save and display
  settings$savedays <- 10 
}


##### 3 DEFINE FUNCTIONS #####

#### 3.1 CALCULATE CORTISOL & ALLOSTATIC LOAD ####
CTA <- function(HPA, fs, halftime, thL, rL, integration = c("Euler", "RK4")[1]){
  #HPA is the vector that describes cortisol production by the HPA-axis
  #dt (=1/fs) is the time interval between samples
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
  C <- L <- A <- as.numeric(c(0, rep(NA,N-1))) #initialize
  
  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- Euler(Y1 = C[t], func = lingrowth, param = c(r, HPA[t]), h = 1/fs)),
         "RK4" =   # Runge-Kutta 4th Order Integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- RK4(Y1 = C[t], func = lingrowth, param = c(r, HPA[t]), h = 1/fs)))
  switch(integration,
         "Euler" =    # Euler numerical integration
           sapply(seq_along(L), function(t) L[[t+1]] <<- Euler(Y1 = L[t], func = lingrowth, param = c(-rL, C[t]), h = 1/fs)),
         "RK4" =   # Runge-Kutta 4th Order Integration
           sapply(seq_along(L), function(t) L[[t+1]] <<- RK4(Y1 = L[t], func = lingrowth, param = c(-rL, C[t]), h = 1/fs)))
  sapply(seq_along(A), function(t) A[[t+1]] <<- ifelse(L[t] > thL, 
                                                       A[t] + (L[t] - thL)*1/fs, 
                                                       A[t]))
  
  return(data.frame(C,L,A))
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

#### 3.6 probability HPAA activity ####
#this function creates a vector with probability values that HPA activity will occur at any possible sample
makepHPA <- function(settings,fsleep,fwork){
  
  #define which of the N samples per day relate to night time HPA activity and which to work time HPA activity
  night <- with(settings, 1:(N*(waketime+CAR)/24)); work <- with(settings, (N*worktime[1]/24+1):(N*worktime[2]/24))
  
  #for non-working days, there is only a linear increasing probability of HPA activity throughout the night, with the cummulative probability = fsleep
  pHPAo <- rep(0,settings$N); pHPAo[night] <- night/sum(night)*fsleep
  
  #for working days, there is also a constrant probability of HPA activity throughout the work day, with the cummulative probability = fwork
  pHPAw <- rep(0,settings$N); pHPAw[night] <- night/sum(night)*fsleep; pHPAw[work] <- 1/length(work)*fwork
  
  #make a vector that combines days of work and non-work, depending of which days are off-days
  pHPA <- stack(as.data.frame(sapply(1:settings$days, function(day) if((day - floor(day/7)*7) %in% settings$dayoff){pHPAo}else{pHPAw})))[1]
  return(pHPA)
}

## 3.7 executing the simulation per person ####
dosim <- function(settings, simperson, n){   
  
  ## create HPA activity
  
  #create a vector with probability values that HPA activity will occur at any possible sample
  pHPA <- makepHPA(settings,fsleep=simperson$wakefM_pp,fwork=simperson$stressfM_pp)
  
  #create HPA activity by comparing each probability value with a random generated number [0:1]
  set.seed(n); HPA <- (ifelse(runif(settings$N*settings$days)<pHPA,1,0))  
  
  ## calculate cortisol, L and allostatic load from the HPA activity
  
  data <- CTA(HPA, settings$fs, simperson$halftime_pp, simperson$thL_pp, simperson$rL_pp, settings$integration)
  
  ## Store results
  
  savedays <- settings$savedays #for practical use...
  
  #store data for creating an ODD's ratio of becoming sick 
  write.csv(data.frame(simperson=rep(simperson$simperson,length(savedays)),days=savedays,stressfM=rep(simperson$stressfM_pp,length(savedays)),Aload=data$A[savedays],Sick=data$A[savedays]>simperson$thA_pp)
            ,file.path(getwd(), subDir, "daydata",paste0("daydata",simperson$simperson)))
  
  # store data for creation of the average cortisol profiles
  write.csv(
    cbind(time=0:(settings$N-1)/settings$N*24,data[((simperson$sampleday-1)*settings$N+1):(simperson$sampleday*settings$N),]),
    file.path(getwd(), subDir, "Qdata",paste0("Qdata",simperson$simperson)))
  
  # store data for of a sample of persons for creating sample profiles
  if(simperson$simperson %in% settings$samplepeople){write.csv(
    cbind(simperson=rep(simperson$simperson,settings$N*7),
          time=0:(settings$N*7-1)/settings$N*24,
          HPA=HPA[1:(settings$N*7)],data[1:(settings$N*7),]),
    file.path(getwd(), subDir, "plotsample",paste0("plotsample",simperson$simperson))
  )}
}

#### 3.8 combine files into a larger file #### 
bindfiles = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE) #identify the files
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)}) #read the files into a list of data frames
  write.csv(
    Reduce(function(x,y) {rbind(x,y)}, datalist), #combine the data frames
    paste0(mypath,".csv")) #store the new file
  unlink(mypath, recursive = T, force = T) #delete folder with originals
}


################################### 4 EXECUTION PART ###################################

#### 4.1 Initialize ####

## time stamp the start of the script execution, so we can determine how long the whole thing took
start.time <- Sys.time()

n <- 1 #seed counter; everytime a random number is drawn, the next random number is drawn using a new seed, but every execution of this script causes the same random numbers to be used


#### 4.2 Create people from population characteristics ####

#create data frame to store all information about individuals
simpopulation <- data.frame(simperson = as.character(paste0('SimPerson_',1:settings$people)) )

## parameters that determine the cortisol profile

##Cortisol decay
set.seed(n); n=n+1; simpopulation$halftime_pp <- abs( rnorm(settings$people, mean = settings$halftimeM, sd = settings$halftimeSD) )

##number of awakening impulses (indivual mean)
set.seed(n); n=n+1; simpopulation$wakefM_pp <- abs( rnorm(settings$people, mean = settings$wakefM, sd = settings$wakefSD) ) #can have only positive values

##number of stress stressors (indivual mean)
set.seed(n); n=n+1; simpopulation$stressfM_pp <- abs( rnorm(settings$people, mean = settings$stressfM, sd = settings$stressfSD) ) #can have only positive values

## parameters that determine how a cortisol profile translates into allostatic load and sickness
set.seed(n); n=n+1; simpopulation$thL_pp <- abs( rnorm(settings$people, mean = settings$thLM, sd = settings$thLSD) ) 
set.seed(n); n=n+1; simpopulation$rL_pp <- abs( rnorm(settings$people, mean = settings$rLM, sd = settings$rLSD) ) 
set.seed(n); n=n+1; simpopulation$thA_pp <- abs( rnorm(settings$people, mean = settings$thAM, sd = settings$thASD) ) 

#### 4.3 Prepare for file saving ####

if (load_from_file == F){
  #select 5 random people of whom the first 7 consequent days will be plotted
  settings$samplepeople <- sample(simpopulation$simperson, ifelse(settings$people >= 5,5,settings$people), replace = FALSE) 
  
  #choose subset of days to plot
  if(settings$days<=settings$savedays){settings$savedays <- 1:settings$days
  } else {settings$savedays <- seq( as.integer(settings$days/settings$savedays), settings$days, as.integer(settings$days/settings$savedays) )}
}

#randomly choose on which day in the live a person is sampled for his daily cortisol fluctuation
set.seed(n); n=n+1; 
simpopulation$sampleday <- sample(
  which(sapply(1:settings$days, function(day) !((day - floor(day/7)*7) %in% settings$dayoff))==T), #only working days
  settings$people, replace = T) 

#create folders for storing data
subDir <- paste("simulation", format(Sys.time(),"%Y-%m-%d_%H-%M"))
dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
dir.create(file.path(getwd(), subDir, 'Qdata'), showWarnings = FALSE)
dir.create(file.path(getwd(), subDir, 'daydata'), showWarnings = FALSE)
dir.create(file.path(getwd(), subDir, 'plotsample'), showWarnings = FALSE)

settings$preptime <- Sys.time() - start.time

#### 4.4 Simulate #### 

start.time <- Sys.time() #initiate time recording

#per person 
l_ply(seq(settings$people), function(id) dosim(settings = settings, simperson = simpopulation[id,], n = n+id))

settings$simtime <- Sys.time() - start.time

################################### 5 STORE DATA ###################################

#combine individual files and delete originals
bindfiles(file.path(getwd(), subDir, "Qdata")) 
bindfiles(file.path(getwd(), subDir, "daydata"))
bindfiles(file.path(getwd(), subDir, "plotsample"))
#write settings 
save(settings, file = file.path(getwd(), subDir, "settings.RData"))
write.csv(simpopulation, file.path(getwd(), subDir, "simpopulation.csv"))

setwd(maindir)

# ################################### 6 MAKE CURVES ###################################

package("rmarkdown")
package("knitr")
render("makegraphs.Rmd",output_file = paste0("simoutput/graphs",subDir,".html"),params = list(datadir = file.path(getwd(), 'simoutput',subDir)))
