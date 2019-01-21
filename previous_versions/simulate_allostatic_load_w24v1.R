##############################################
# This script ..............................
# # ----CREATING CORTISOL PROFILE BY SIMULATING SUDDEN CORTISOL PRODUCTION BURST USING SIMPLE FORMULA----
# # The circardian cortisol profile will be constructed as a periodically varying cortisol production of the HPA-axis,
# minus decay of cortisol using:
# ..............................
##############################################

rm(list = ls()) #clear the environment

######## 0.1 DEFINE PARAMETERS VARYING BETWEEN SIMULATIONS ########

# General explaination of part 0.1
# this allows performing simulations with different parameters and thus compare the effect of changing the parameters
# each list within 'paramsets' refers to a complete simulation of a number of people (e.g. 10'000) and days (e.g. 1000)

#****USER INSTR.: CHANGE HERE THE PARAMETERS THAT SHOULD VARY BETWEEN SIMULATIONS****#


paramsets <- list(
  list(#simulation 1
    worktime = c(9,17),
    anticipation = c("YES","NO")[1],
    dayoff = c(0,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 2
    worktime = c(9,19),
    anticipation = c("YES","NO")[1],
    dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
   ))


######## 0.2 DEFINE DEFAULT SETTINGS ########

#****USER INSTR.: SET HERE THE PARAMETERS AS THEY OCCUR FOR ALL SIMULATIONS****#

# General explaination of part 0.2:
# this defines the default settings for each simulation. Variables that have been defined above in 'paramsets' 
# overwrite specific default values. Those variables not defined in 'paramsets' will receive their default value

defaults <- list( #initialize list

    ## simulation settings
    fs = 2, #samples per hour
    days = 10, #simulated number of days per person
    people = 500, #number of persons simulated
    integration = "RK4", #numerical integration method used

    #night impulses
    night.f.mean = 8, #average number of night impulses 
    night.dshape = 1, #shape factor of the gamma distribution of impulse frequencies
    waketime = 7,

    #work stressors
    distribution = c("gamma","uniform")[1], #distribution type describing the variation in the number of work impulses 
    #between individuals (gamma being representative of the population; a uniform distribution representing an equal amount of people for each number of work impulses)
    work.f.mean = 6, #average number of work impulses per 8h, only relevant when the gamma distribution is chosen
    work.dshape = 3, #shape factor of the gamma distribution of impulse frequencies - only relevant for gamma
    work.unif = c(1,50), #range of work impulses - only relevant for uniform distribution
    worktime = c(9,17),
    anticipation = c("YES","NO")[1], #whether or not anticipation of work day stress is assumed
    dayoff = c(0,6), #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'

    ## physiological parameters
    kHPA = 2.193871, #multiplication constant describing the cortisol release in response to a neural impulse 
    delay = .5, #delay between neural impulse and the cortisol release
    halflifeM = 80/60, #average half life of cortisol is 80 minutes or 4/3hour
    halflifeSD = 8/60, #between person variation
    eps.A = 20, #threshold of elastic load becoming allostatic load 
    thLSD = 0, #TO BE REMOVED
    rho.E = .6, #the elasticity constant - rate of recovery from elastic load
    rLSD = 0, #TO BE REMOVED
    eps.D = 400, #threshold of allostatic load causing disease
    thASD = 0, #TO BE REMOVED

    ## number of days to save and display
    savedays = 5
)


#### 1 LOAD PACKAGES & SET DIRECTORY #####

start0<-Sys.time() #to track the duration of the computations

#function to load a library and install it when it is not yet installed
package <- function(pack) { #function to load a library and install it when it is not yet installed
  if (!as.character(pack) %in% installed.packages()) {install.packages(as.character(pack)) }
  library(pack,character.only=TRUE)}

package("lattice")
package("plyr")
package("dplyr")
package("ggplot2")
package("tcltk")
package("car")
package("rmarkdown")
package("knitr")

## set saving directories

#checks if a dedicated data output folder is available in the working directory and, if not, ask the user to define the directory
winDir <- "D/Surfdrive/BS28A - Major project/simulations"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 
if (!file.exists("simoutput")){setwd(tk_choose.dir(getwd(), "Choose folder for storing simulation data"))
  if (!file.exists("simoutput")){dir.create("simoutput")}} 

#define data storage locations
if(length(paramsets)<2){
  datadir=file.path("simoutput",paste0("sim", format(Sys.time(),"%Y-%m-%d.%Hh%M")))
}else{
  compdir<-file.path("simoutput",paste0("comp", format(Sys.time(),"%Y-%m-%d.%Hh%M")))
  datadir<-file.path(compdir,paste0("sim",seq(length(paramsets))))
  dir.create(compdir)
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
  set.seed(n); n=n+1; simpopulation$halflife_pp <- abs( rnorm(settings$people, mean = settings$halflifeM, sd = settings$halflifeSD) )
  
  ##number of awakening impulses (indivual mean)
  set.seed(n); n=n+1; 
  simpopulation$night.f.mean_pp <- rgamma(settings$people, shape = settings$night.dshape, rate = settings$night.dshape/settings$night.f.mean)
  
  ##number of stress stressors (indivual mean)
  set.seed(n); n=n+1; 
  simpopulation$work.f.mean_pp <- switch(settings$distribution, 
                                      gamma = rgamma(settings$people, shape = settings$work.dshape, rate = settings$work.dshape/settings$work.f.mean), 
                                      uniform = runif(settings$people, min = settings$stressunif[1], max = settings$stressunif[2]))
  
  ## parameters that determine how a cortisol profile translates into allostatic load and sickness
  set.seed(n); n=n+1; simpopulation$thL_pp <- abs( rnorm(settings$people, mean = settings$eps.A, sd = settings$thLSD) ) 
  set.seed(n); n=n+1; simpopulation$rL_pp <- abs( rnorm(settings$people, mean = settings$rho.E, sd = settings$rLSD) ) 
  set.seed(n); n=n+1; simpopulation$thA_pp <- abs( rnorm(settings$people, mean = settings$eps.D, sd = settings$thASD) ) 
  
  #### 2.4 Prepare for file saving ####
  
  simlog <- list(preptime=NULL, simtime = NULL, filemergetime = NULL) #create a data frame for logging information about the simulation
  
  
  #select 5 random people of whom the first 7 consequent days will be plotted
  settings$samplepeople <- sample(simpopulation$simperson, ifelse(settings$people >= 5,5,settings$people), replace = FALSE) 
  
  #choose subset of days to plot
  if(settings$days<=settings$savedays){settings$savedays <- 1:settings$days
  } else {settings$savedays <- seq( as.integer(settings$days/settings$savedays), settings$days, as.integer(settings$days/settings$savedays) )}
  
  #randomly choose on which day in the live a person is sampled for his daily cortisol fluctuation
  set.seed(n); n=n+1; 
  simpopulation$sampleday <- sample(
    which(sapply(1:settings$days, function(day) !((day - floor(day/7)*7) %in% settings$dayoff))==T), #only working days
    settings$people, replace = T) 
  
  #create folders for storing data
  dir.create(subDir)
  dir.create(file.path(subDir, 'Qdata'))
  dir.create(file.path(subDir, 'daydata'))
  dir.create(file.path(subDir, 'plotsample'))
  
  simlog$preptime <-  difftime(Sys.time(), start.time, unit="mins")
  
  #### 2.5 Simulate #### 
  
  start.time <- Sys.time() #initiate time recording
  
  #determine the scaling factor for translating (reference to the average peak at 9mmol/L in Miller)
  if(is.na(settings$kHPA)){settings$kHPA<-8/mean(ldply(seq(settings$people), function(id) getscale(settings = settings, simperson = simpopulation[id,], n = n+id))[,1])}
  
  #perform the simulation per person and return several evaluation values
  evalv <- sapply(seq(settings$people), function(id) dosim(settings = settings, simperson = simpopulation[id,], n = n+id, subDir = subDir))
  
  simlog$kHPA<-settings$kHPA
  simlog$recovery<-mean(unlist(evalv[1,])) # ratio of days on which full recovery occurs
  simlog$aloadratio<-mean(unlist(evalv[2,])) # ratio people developing allostatic load
  simlog$sickratio<-mean(unlist(evalv[3,])) # ratio people developing allostatic load
  
  simlog$simtime <- difftime(Sys.time(), start.time, unit="mins")
  
  ################################### 2.6 STORE DATA ###################################
  
  #combine individual files and delete originals
  start.time <- Sys.time()
  bindfiles(file.path(getwd(), subDir, "Qdata"))
  bindfiles(file.path(getwd(), subDir, "daydata"))
  bindfiles(file.path(getwd(), subDir, "plotsample"))
  write.csv(simpopulation, file.path(getwd(), subDir, "simpopulation.csv"))
  simlog$filemergetime <- difftime(Sys.time(), start.time, unit="mins")
  
  #write settings
  save(settings, file = file.path(getwd(), subDir, "settings.RData"))
  #write log
  save(simlog, file = file.path(getwd(), subDir, "simlog.RData"))
  
  ################################### 2.7 MAKE CURVES (SINGLE SIMULATION) ###################################
  
  #only create a markdown file for the single simulation if this execution concerns only a single simulation
  if(!grepl("comp",subDir, fixed = T)){  render("makegraphs.Rmd",output_file = paste0(subDir,".html"),params = list(datadir = subDir, paramsets = NULL, totalsimtime=difftime(Sys.time(), start0, unit="mins")))}
}

##### 3 DEFINE EMBEDDED FUNCTIONS #####

#### 3.1 CALCULATE CORTISOL & ALLOSTATIC LOAD ####
CTA <- function(HPA, fs, delay, halflife, thL, rL, integration = c("Euler", "RK4")[1], kHPA=1){
  #HPA is the vector that describes cortisol production by the HPA-axis
  #dt (=1/fs) is the time interval between samples
  #halflife is the cortisol half-time
  #!! note: make sure that dt and halflife are given in the same unit, e.g. h
  
  #calculate decay rate
  #standard decay function: dY/dt = Y*r (decay is a function of the current concentration)
  #one can rewrite this as: dY/Y = r*dt. Integrating over time gives ln(Y(t))-ln(Y(0)) = r*t or Y(t)/Y(0) = exp(r*t) or Y(t) = Y(0)*exp(r*t)
  #filling in: Y(t=1/2) = 1/2*Y(0) = Y(0)*exp(r*t_halflife) or -ln(2) = r*t_halflife or r = -ln(2)/t_halflife
  r <- -log(2)/halflife
  
  #calculate Y[t+1] from dY/dt = Y*r and additional input (cortisol produced by the HPA axis)
  #dY/dt = Y*r + input or dY = (Y*r + input)*dt.
  #this can be integrated numerically using Euler: Y[t+1] = (Y[t]*r + input[t])*dt + Y[t], but also using RK4 (see below)
  
  #first introduce the delay between the impulses to the HPAA and its cortisol output (for practical reasons, we do this outside the formula below)
  HPAout <- c(rep(0,round(delay*fs)),HPA[1:(length(HPA)-round(delay*fs))])*kHPA #
  
  #initiate and add values start values (last values of the previous day)
  N <- length(HPAout)
  C <- L <- A <- as.numeric(c(0, rep(NA,N-1))) #initialize
  
  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- Euler(Y1 = C[t], func = lingrowth, param = c(r, HPAout[t]), h = 1/fs)),
         "RK4" =   # Runge-Kutta 4th Order Integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- RK4(Y1 = C[t], func = lingrowth, param = c(r, HPAout[t]), h = 1/fs)))
  
  if(!is.na(kHPA)){#only if the scaling factor is determined will the calculation of L and A continue, otherwise values of C are first use to determine this factor
    switch(integration,
           "Euler" =    # Euler numerical integration
             sapply(seq_along(L), function(t) L[[t+1]] <<- Euler(Y1 = L[t], func = lingrowth, param = c(-rL, C[t]), h = 1/fs)),
           "RK4" =   # Runge-Kutta 4th Order Integration
             sapply(seq_along(L), function(t) L[[t+1]] <<- RK4(Y1 = L[t], func = lingrowth, param = c(-rL, C[t]), h = 1/fs)))
    sapply(seq_along(A), function(t) A[[t+1]] <<- ifelse(L[t] > thL, 
                                                         A[t] + (L[t] - thL)*1/fs, 
                                                         A[t]))
    
    return(data.frame(C,L,A))
  }else{return(C)}
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
  Y2 <- as.numeric(Y1)*r + as.numeric(input)
  return(Y2)
}

#### 3.6 create day HPA activity ####
#this function creates a vector with probability values that HPA activity will occur at any possible sample
makeHPA <- function(day, night.f, fstress, settings,n){
  workdayfactor <- (settings$worktime[[2]]-settings$worktime[[1]])/8 #corrects in the workstressors for longer or shorter working days
  set.seed(n)
  if((day - floor(day/7)*7) %in% settings$dayoff) #for non-working days there is only night time HPA activity
  {set.seed(n); spikes <- abs(-rexp(night.f) + with(settings, waketime))
  }else{ #for working days, there is also HPA activity during work time
    spikes <- switch(settings$anticipation,
                     NO = c(abs(-rexp(night.f) + with(settings, waketime)),
                            runif(fstress*workdayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])),
                     YES = c(abs(-rexp(night.f+fstress*workdayfactor) + with(settings, waketime)),
                             runif(fstress*workdayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])))
  }
  
  spikes <- ceiling(spikes*settings$fs) #corresponding sample moments of the spikes (and prevent spikes at time point 0)
  
  HPA<-rep(0, settings$fs*24)
  HPA[unique(spikes)]<-sapply(unique(spikes), function(spike) sum(spikes==spike))
  return(HPA)
}

## 3.7 executing the simulation per person ####
dosim <- function(settings, simperson, n, subDir){   
  
  ## create HPA activity
  
  #determine the number of impulses that the simulated fictive individual experiences on a particular day
  set.seed(n)
  night.f <- rpois(settings$days, simperson$night.f.mean_pp)
  set.seed(n+1e6)
  fstress <- rpois(settings$days, simperson$work.f.mean_pp)
  
  #create a vector with probability values that HPA activity will occur at any possible sample
  HPA <- stack(as.data.frame(sapply(seq(settings$days), function(day) makeHPA(day, night.f[day], fstress[day], settings,n=n+simperson$simperson*day))))[,1]
  
  ## calculate cortisol, L and allostatic load from the HPA activity
  
  data <- CTA(HPA, settings$fs, settings$delay, simperson$halflife_pp, simperson$thL_pp, simperson$rL_pp, settings$integration, settings$kHPA)

  ## Store results
  
  savedays <- settings$savedays #for practical use...
  
  #store data for creating an ODD's ratio of becoming sick 
  write.csv(data.frame(simperson=rep(simperson$simperson,length(savedays)),days=savedays,work.f.mean=rep(simperson$work.f.mean_pp,length(savedays)),Aload=data$A[savedays*settings$fs*24],Sick=data$A[savedays*settings$fs*24]>simperson$thA_pp)
            ,file.path(getwd(), subDir, "daydata",paste0("daydata",simperson$simperson)))
  
  # store data for creation of the average cortisol profiles
  write.csv(
    cbind(time=seq(0,24,1/settings$fs),data[1:(settings$fs*24+1),]),
    file.path(getwd(), subDir, "Qdata",paste0("Qdata",simperson$simperson)))
  
  # store data for of a sample of persons for creating sample profiles
  if(simperson$simperson %in% settings$samplepeople){write.csv(
    cbind(simperson=rep(simperson$simperson,settings$fs*24*7+1),
          time=seq(0,24*7,1/settings$fs),
          HPA=HPA[1:(settings$fs*24*7+1)],data[1:(settings$fs*24*7+1),]),
    file.path(getwd(), subDir, "plotsample",paste0("plotsample",simperson$simperson))
  )}
  
  #count the ratio of days where there was dynamic load remaining at the end of the day
  recovery <- sum(data$L[seq(settings$days)*settings$fs*24] < rep(max(data$L)/20,settings$days))/settings$days
  
  #whether or not the person developed allostatic load during the simulation
  aloadYN <- data$A[nrow(data)]>0
  
  sickYN <- data$A[length(data$A)]>simperson$thA_pp
  
  return(data.frame(recovery,aloadYN,sickYN))
}


## 3.8 calibrating kHPA ####
getscale <- function(settings, simperson, n){   
  
  ## create HPA activity
  
  #create a vector with probability values that HPA activity will occur at any possible sample
  HPA <- makeHPA(day=1, night.f[day], fstress[day], settings,n=n+simperson$simperson)
  #HPA <- stack(as.data.frame(sapply(seq(1), function(day) makeHPA(day, simperson, settings,n=n+simperson$simperson*day))))[1]
  
  ## calculate cortisol from the HPA activity
  
  data <- CTA(HPA[1:(with(settings, (waketime+delay)*fs)+1)], settings$fs, simperson$halflife_pp, simperson$thL_pp, simperson$rL_pp, settings$integration, settings$kHPA)
  #data <- CTA(HPA, settings$fs, simperson$halflife_pp, simperson$thL_pp, simperson$rL_pp, settings$integration, settings$kHPA)
  #print(which(data==max(data)))
  
  #return(data[with(settings, (waketime+CAR)*fs)])
  return((data)[length(data)])
}

#### 3.9 combine files into a larger file #### 
bindfiles = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE) #identify the files
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)}) #read the files into a list of data frames
  write.csv(
    Reduce(function(x,y) {rbind(x,y)}, datalist), #combine the data frames
    paste0(mypath,".csv")) #store the new file
  unlink(mypath, recursive = T, force = T) #delete folder with originals
}

#### 4 RUN THE SIMULATION FOR THE GIVEN PARAMETER SETS ####

start.time<-Sys.time()

sapply(seq(length(paramsets)), function(i) exesim(paramsets[[i]],defaults,datadir[i]))

difftime(Sys.time(), start.time, unit="mins")

totalsimtime=round(as.numeric(difftime(Sys.time(), start0, unit="mins"))) #in minutes

#if it concerns a comparison, store comparison information and render comparison graph
if(length(datadir)>1){  
  save(paramsets, file = file.path(compdir, "paramsets.RData"))
  save(datadir, file = file.path(compdir, "datadir.RData"))
  render("makegraphs.Rmd",output_file = paste0(compdir,".html"),params = list(datadir = datadir, paramsets=paramsets, totalsimtime=totalsimtime))
  }
