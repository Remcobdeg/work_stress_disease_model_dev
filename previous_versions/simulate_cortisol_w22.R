##############################################
# This script create has the goal of showing the impact of various population distributions
# and stressor frequencies on the population average
##############################################

rm(list = ls()) #clear the environment

######## 0.1 DEFINE PARAMETERS VARYING BETWEEN SIMULATIONS ########

#****USER INSTR.: CHANGE HERE THE PARAMETERS THAT VARY BETWEEN SIMULATIONS****#

paramsets <- list(
  list(#simulation 1
    worktime = c(9,17),
    distribution = c("population","uniform","chisq")[3],
  )#,list(#simulation 2
  #   worktime = c(8,17),
  #   distribution = c("population","uniform")[2],
  #   dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  # ),list(#simulation 3
  #   worktime = c(8,19),
  #   distribution = c("population","uniform")[2],
  #   dayoff = c(0,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  # )
)

######## 0.2 DEFINE DEFAULT SETTINGS ########

#****USER INSTR.: SET HERE THE PARAMETERS AS THEY OCCUR FOR ALL SIMULATIONS****#

defaults <- list( #initialize list
  
  ## simulation settings
  fs = 10, #samples per hour
  N = NULL, #number of samples in one day - dummy, this value is defined later
  days = 1,
  people = 1000,
  integration = "RK4",
  
  #awakening impulse
  wakefM = 20, #per day
  wakefSD = 20,
  waketime = 7,
  CAR = .5, #duration (in hours) of the cortisol awakening response
  nightgrowxp = 1,#exponential growth factor 
  
  #work stressors
  distribution = c("population","uniform","chisq")[3],
  stressfM = 20, #per 8h
  stressfSD = 20, #frequency of stress stressors
  stressunif = c(0,180),
  worktime = c(9,17),
  anticipation = c("YES","NO")[1],
  dayoff = c(0,6), #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  
  ## physiological parameters
  kHPA = NA, #set to NA to have to code determine the scaling factor
  halftimeM = 80/60, #half time of cortisol is 80 minutes or 4/3hour
  halftimeSD = 8/60, #between person variation
  thLM = 20, #threshold of dynamic load becoming allostatic load 
  thLSD = 0,
  rLM = .6, #rate at which the body recovers form the dynamic load
  rLSD = 0,
  thAM = 1000, #threshold of allostatic load causing disease
  thASD = 0,
  
  ## number of days to save and display
  savedays = 5
)

#### 1 LOAD PACKAGES & SET DIRECTORY #####

start.time = Sys.time()

#function to load a library and install it when it is not yet installed
package <- function(pack) { #function to load a library and install it when it is not yet installed
  if (!as.character(pack) %in% installed.packages()) {install.packages(as.character(pack)) }
  library(pack,character.only=TRUE)}

package("lattice")
package("plyr")
package("dplyr")
package("tidyr")
package("ggplot2")
package("tcltk")
package("car")
package("rmarkdown")
package("knitr")



#checks if a dedicated data output folder is available in the working directory and, if not, ask the user to define the directory
winDir <- "D/Surfdrive/BS28A - Major project/simulations"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 
if (!file.exists("simoutput")){setwd(tk_choose.dir(getwd(), "Choose folder for storing simulation data"))
  if (!file.exists("simoutput")){dir.create("simoutput")}} 

#define data storage locations
if(length(paramsets)<2){
  datadir=file.path("simoutput",paste0("profile", format(Sys.time(),"%Y-%m-%d.%Hh%M")))
}else{
  compdir<-file.path("simoutput",paste0("profiles", format(Sys.time(),"%Y-%m-%d.%Hh%M")))
  datadir<-file.path(compdir,paste0("sim",seq(length(paramsets))))
  dir.create(compdir)
}

##### 3 DEFINE EMBEDDED FUNCTIONS #####

#### 3.1 CALCULATE CORTISOL & ALLOSTATIC LOAD ####
HPAA <- function(HPA, fs, halftime, integration = c("Euler", "RK4")[1]){
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
  C <- as.numeric(c(0, rep(NA,N-1))) #initialize
  
  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- Euler(Y1 = as.numeric(C[t]), func = lingrowth, param = c(r, HPA[t]), h = 1/fs)),
         "RK4" =   # Runge-Kutta 4th Order Integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- RK4(Y1 = as.numeric(C[t]), func = lingrowth, param = c(r, HPA[t]), h = 1/fs)))
  return(C)
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
  Y2 <- as.numeric(Y1)*as.numeric(r) + as.numeric(input)
  return(as.numeric(Y2))
}

#### 3.6 create day HPA activity ####
#this function creates a vector with probability values that HPA activity will occur at any possible sample
makeHPA <- function(day, simperson, settings,n){
  wordayfactor <- (settings$worktime[[2]]-settings$worktime[[1]])/8 #corrects in the workstressors for longer or shorter working days
  set.seed(n)
  if((day - floor(day/7)*7) %in% settings$dayoff) #for non-working days there is only night time HPA activity
  {set.seed(n); spikes <- abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR))
  }else{ #for working days, there is also HPA activity during work time
    spikes <- switch(settings$anticipation,
                     NO = c(abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR)),
                            runif(simperson$stressfM_pp*wordayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])),
                     YES = c(abs(-rexp(simperson$wakefM_pp+simperson$stressfM_pp*wordayfactor, rate = settings$nightgrowxp) + with(settings, waketime + CAR)),
                             runif(simperson$stressfM_pp*wordayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])))
  }
  spikes <- ceiling(spikes*settings$fs) #corresponding sample moments of the spikes (and prevent spikes at time point 0)
  
  HPA<-rep(0, settings$N)
  HPA[unique(spikes)]<-sapply(unique(spikes), function(spike) sum(spikes==spike))
  return(HPA)
}

## 3.7 executing the simulation per person ####
dosim <- function(settings, simperson, n, subDir){   
  
  ## create HPA activity
  
  #create a vector with probability values that HPA activity will occur at any possible sample
  HPA <- stack(as.data.frame(sapply(seq(settings$days), function(day) makeHPA(day, simperson, settings,n=n+simperson$simperson*day))))[,1]
  
  ## calculate cortisol
  return(HPAA(HPA, settings$fs, simperson$halftime_pp, settings$integration))
}



######## 2 EXECUTION FUNCTION ########

exesim <- function(paramset,defaults,subDir){

  #### 2.1 Initialization steps ####
  
  ## time stamp the start of the script execution, so we can determine how long the whole thing took
  start.time <- Sys.time()
  
  n <- 1 #seed counter; everytime a random number is drawn, the next random number is drawn using a new seed, but every execution of this script causes the same random numbers to be used
  
  ######## 2.2.2 DEFINE SETTINGS FROM PROVIDED PARAMETERS AND DEFAUTLS ########
  
  #define the settings
  settings <- sapply(names(defaults), function(param) ifelse(param %in% names(paramset),paramset[param],defaults[param]))
  
  #add N for practical purposes
  settings$N <- 24*settings$fs #number of samples in one day
  
  #### 2.3 Create individuals ####
  
  #create data frame to store all information about individuals
  simpopulation <- data.frame(simperson = seq(settings$people) )
  
  ## parameters that determine the cortisol profile
  
  ##Cortisol decay
  set.seed(n); n=n+1; simpopulation$halftime_pp <- abs( rnorm(settings$people, mean = settings$halftimeM, sd = settings$halftimeSD) )
  
  ##number of awakening impulses (indivual mean)
  set.seed(n); n=n+1; 
  simpopulation$wakefM_pp <- switch(settings$distribution, 
                                    population = abs( rnorm(settings$people, mean = settings$wakefM, sd = settings$wakefSD)-1 )+1, #can have only positive values
                                    uniform = settings$wakefM/3*rchisq(settings$people, df=4)+1,
                                    chisq = settings$wakefM/3*rchisq(settings$people, df=4)+1)
  
  ##number of stress stressors (indivual mean)
  set.seed(n); n=n+1; 
  simpopulation$stressfM_pp <- switch(settings$distribution, 
                                      population = abs( rnorm(settings$people, mean = settings$stressfM, sd = settings$stressfSD) ), #can have only positive values
                                      uniform = runif(settings$people, min = settings$stressunif[1], max = settings$stressunif[2]),
                                      chisq = settings$stressfM/3*rchisq(settings$people, df=4))
  
  
  #### 2.4 Prepare for file saving ####
  
  simlog <- list(preptime=NULL, simtime = NULL) #create a data frame for logging information about the simulation
  
  #select 5 random people of whom individual profiles will be plotted
  settings$samplepeople <- sample(simpopulation$simperson, ifelse(settings$people >= 5,5,settings$people), replace = FALSE) 
  
  simlog$preptime <-  difftime(Sys.time(), start.time, unit="mins")

  #create folders for storing data
  dir.create(subDir)


  #### 2.5 Simulate #### 
  
  #perform the simulation per person and return several evaluation values
  Qdata <- as.data.frame(sapply(seq(settings$people), function(id) dosim(settings = settings, simperson = simpopulation[id,], n = n+id, subDir = subDir)))
  Qdata <- cbind(time=0:(settings$N)/settings$N*24,Qdata)
  Qdata <- gather(Qdata, key = person, value = C, -time) #long data
  
  simlog$simtime <- difftime(Sys.time(), start.time, unit="mins")
  
  ################################### 2.6 STORE DATA ###################################
  
  save(settings, file = file.path(getwd(), subDir, "settings.RData"))
  save(simlog, file = file.path(getwd(), subDir, "simlog.RData"))
  save(simpopulation, file = file.path(getwd(), subDir, "simpopulation.RData"))
  save(Qdata, file = file.path(getwd(), subDir, "Qdata.RData"))
}

#### 4 RUN THE SIMULATION FOR THE GIVEN PARAMETER SETS ####

start.time<-Sys.time()

sapply(seq(length(paramsets)), function(i) exesim(paramsets[[i]],defaults,datadir[i]))

difftime(Sys.time(), start.time, unit="mins")



################################### 2.7 MAKE CURVES  ###################################

if(length(datadir)>1){  
  save(paramsets, file = file.path(compdir, "paramsets.RData"))
  save(datadir, file = file.path(compdir, "datadir.RData")) 
  render("makeprofiles.Rmd",output_file = paste0(compdir,".html"),params = list(datadir = datadir, paramsets=paramsets))
} else {render("makeprofiles.Rmd",output_file = paste0(datadir,".html"),params = list(datadir = datadir, paramsets=paramsets))}


