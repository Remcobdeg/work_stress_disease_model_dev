##############################################
# This script create has the goal of showing the impact of various population distributions
# and stressor frequencies on the collective cortisol profile that is simulated
##############################################

rm(list = ls()) #clear the environment

######## 0.1 DEFINE PARAMETERS VARYING BETWEEN SIMULATIONS ########

#****USER INSTR.: CHANGE HERE THE PARAMETERS THAT VARY BETWEEN SIMULATIONS****#

paramsets <- list(
  list(#simulation 1
  #   wakefM = 8, #per day
  #   dshape.wake = 1.2,
  #   anticipation = c("YES","NO")[1],
  #   stressfM = 6, #per 8h
  #   dshape.stress = 3
  # ),list(#simulation 2
    wakefM = 14, #per day
    dshape.wake = 3,
    anticipation = c("YES","NO")[2],
    stressfM = 6, #per 8h
    dshape.stress = 3
  ),list(#simulation 2
    wakefM = 8, #per day
    dshape.wake = 1,
    anticipation = c("YES","NO")[1],
    stressfM = 6, #per 8h
    dshape.stress = 3
  # ),list(#simulation 2
  #   wakefM = 8, #per day
  #   dshape.wake = 1.1,
  #   anticipation = c("YES","NO")[1],
  #   stressfM = 6, #per 8h
  #   dshape.stress = 2
    )
)

######## 0.2 DEFINE DEFAULT SETTINGS ########

#****USER INSTR.: SET HERE THE PARAMETERS AS THEY OCCUR FOR ALL SIMULATIONS****#

defaults <- list( #initialize list
  
  ## simulation settings
  fs = 2, #samples per hour
  days = 1,
  people = 10000,
  integration = "RK4",
  
  #awakening impulse
  wakefM = 3, #per day
  dshape.wake = 2,
  waketime = 7,
  nightgrowxp = 1,#exponential growth factor 
  
  #work stressors
  stressfM = 3, #per 8h
  dshape.stress = 2,
  worktime = c(9,17),
  anticipation = c("YES","NO")[1],
  dayoff = c(0,6), #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  
  ## physiological parameters
  delay = .5,
  kHPA = NA, #set to NA to have to code determine the scaling factor
  halftimeM = 80/60, #half time of cortisol is 80 minutes or 4/3hour
  halftimeSD = 8/60 #between person variation
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
macDir <- "/Users/remcobenthemdegrave/OneDrive - Newcastle University/Courses and learning/BS28A - Major project/simulations"
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
HPAA <- function(HPA, fs, delay, halftime, integration = c("Euler", "RK4")[1]){
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
  
  #first introduce the delay between the impulses to the HPAA and its cortisol output (for practical reasons, we do this outside the formula below)
  HPAout <- c(rep(0,round(delay*fs)),HPA[1:(length(HPA)-round(delay*fs))]) #

  #initiate array and add values start values (last values of the previous day)
  N <- length(HPAout)
  C <- as.numeric(c(0, rep(NA,N-1))) #initialize

  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- Euler(Y1 = as.numeric(C[t]), func = lingrowth, param = c(r, HPAout[t]), h = 1/fs)),
         "RK4" =   # Runge-Kutta 4th Order Integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- RK4(Y1 = as.numeric(C[t]), func = lingrowth, param = c(r, HPAout[t]), h = 1/fs)))
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
makeHPA <- function(day, fwake, fstress, settings,n){
  wordayfactor <- (settings$worktime[[2]]-settings$worktime[[1]])/8 #corrects in the workstressors for longer or shorter working days
  set.seed(n)
  if((day - floor(day/7)*7) %in% settings$dayoff) #for non-working days there is only night time HPA activity
  {set.seed(n); spikes <- abs(-rexp(fwake, rate = settings$nightgrowxp) + settings$waketime)
  }else{ #for working days, there is also HPA activity during work time
    spikes <- switch(settings$anticipation,
                     NO = c(abs(-rexp(fwake, rate = settings$nightgrowxp) + settings$waketime),
                            runif(fstress*wordayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])),
                     YES = c(abs(-rexp(fwake+fstress*wordayfactor, rate = settings$nightgrowxp) + settings$waketime),
                             runif(fstress*wordayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])))
  }
  
  spikes <- ceiling(spikes*settings$fs) #corresponding sample moments of the spikes (and prevent spikes at time point 0)
  
  HPA<-rep(0, settings$fs*24)
  HPA[unique(spikes)]<-sapply(unique(spikes), function(spike) sum(spikes==spike))
  return(HPA)
}

## 3.7 executing the simulation per person ####
dosim <- function(settings, simperson, n, subDir){   
  
  ## creating HPA activity
  
  #determine the number of impulses that the simulated fictive individual experiences on a particular day
  set.seed(n)
  fwake <- rpois(settings$days, simperson$wakefM_pp)
  set.seed(n+1e6)
  fstress <- rpois(settings$days, simperson$stressfM_pp)
  
  #create a vector with probability values that HPA activity will occur at any possible sample
  HPA <- stack(as.data.frame(sapply(seq(settings$days), function(day) makeHPA(day, fwake[day], fstress[day], settings,n=n+simperson$simperson*day))))[,1]
  
  ## calculate cortisol
  return(HPAA(HPA, settings$fs, settings$delay, simperson$halftime_pp, settings$integration))
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
  
  #### 2.3 Create individuals ####
  
  #create data frame to store all information about individuals
  simpopulation <- data.frame(simperson = seq(settings$people) )
  
  ## parameters that determine the cortisol profile
  
  ##Cortisol decay
  set.seed(n); n=n+1; simpopulation$halftime_pp <- abs( rnorm(settings$people, mean = settings$halftimeM, sd = settings$halftimeSD) )
  
  ##number of awakening impulses (indivual mean)
  set.seed(n); n=n+1; 
  simpopulation$wakefM_pp <- rgamma(settings$people, shape = settings$dshape.wake, rate = settings$dshape.wake/settings$wakefM)

  ##number of stress stressors (indivual mean)
  set.seed(n); n=n+1; 
  simpopulation$stressfM_pp <- rgamma(settings$people, shape = settings$dshape.stress, rate = settings$dshape.stress/settings$stressfM)
  
  #### 2.4 Prepare for file saving ####
  
  simlog <- list(preptime=NULL, simtime = NULL) #create a data frame for logging information about the simulation
  
  simlog$preptime <-  difftime(Sys.time(), start.time, unit="mins")

  #create folders for storing data
  dir.create(subDir)

  #### 2.5 Simulate #### 
  
  #perform the simulation per person and return several evaluation values
  Qdata <- as.data.frame(sapply(seq(settings$people), function(id) dosim(settings = settings, simperson = simpopulation[id,], n = n+id, subDir = subDir)))
  Qdata <- cbind(time=seq(0,24*settings$days,1/settings$fs),Qdata)
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

l_ply(seq(length(paramsets)), function(i) exesim(paramsets[[i]],defaults,datadir[i]))

difftime(Sys.time(), start.time, unit="mins")

################################### 2.7 MAKE CURVES  ###################################

if(length(datadir)>1){  
  save(paramsets, file = file.path(compdir, "paramsets.RData"))
  save(datadir, file = file.path(compdir, "datadir.RData")) 
  render("makeprofiles.Rmd",output_file = paste0(compdir,".html"),params = list(datadir = datadir, paramsets=paramsets))
} else {render("makeprofiles.Rmd",output_file = paste0(datadir,".html"),params = list(datadir = datadir, paramsets=paramsets))}


