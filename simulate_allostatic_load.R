##############################################
# This script ..............................
# # ----CREATING CORTISOL PROFILE BY SIMULATING SUDDEN CORTISOL PRODUCTION BURST USING SIMPLE FORMULA----
# # The circardian cortisol profile will be constructed as a periodically varying cortisol production of the HPA-axis,
# minus decay of cortisol using:
# ..............................
# updates from previous version (simulate_allostatic_load_w24v2.R):
# - updated macDir 
##############################################

rm(list = ls()) #clear the environment

######## 0.1 DEFINE PARAMETERS VARYING BETWEEN SIMULATIONS ########

# General explaination of part 0.1
# this allows performing simulations with different parameters and thus compare the effect of changing the parameters
# each list within 'paramsets' refers to a complete simulation of a number of people (e.g. 10'000) and days (e.g. 1000)

#****USER INSTR.: CHANGE HERE THE PARAMETERS THAT SHOULD VARY BETWEEN SIMULATIONS****#


paramsets <- list(
  list(#simulation 1
    eps.A = 10,  
    rho.E = .3, 
    worktime = c(9,9+40/5),
    dayoff = c(0,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 2
    eps.A = 10,  
    rho.E = .3, 
    worktime = c(9,9+40/5),
    dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 3
    eps.A = 10,  
    rho.E = .3, 
    worktime = c(9,9+40/4),
    dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 4
    eps.A = 10,  
    rho.E = .3, 
    worktime = c(9,9+40/3),
    dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 5
    eps.A = 10,  
    rho.E = .3, 
    worktime = c(9,9+40/6),
    dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 6
    eps.A = 10,  
    rho.E = .3, 
    worktime = c(9,9+40/7),
    dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'

  ),list(#simulation 7
    eps.A = 10,  
    rho.E = .6, 
    worktime = c(9,9+40/5),
    dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 8
    eps.A = 10,  
    rho.E = .6, 
    worktime = c(9,9+40/5),
    dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 9
    eps.A = 10,  
    rho.E = .6, 
    worktime = c(9,9+40/4),
    dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 10
    eps.A = 10,  
    rho.E = .6, 
    worktime = c(9,9+40/3),
    dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 11
    eps.A = 10,  
    rho.E = .6, 
    worktime = c(9,9+40/6),
    dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 12
    eps.A = 10,  
    rho.E = .6, 
    worktime = c(9,9+40/7),
    dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    
  ),list(#simulation 13
    eps.A = 10,  
    rho.E = .9, 
    worktime = c(9,9+40/5),
    dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 8
    eps.A = 10,  
    rho.E = .9, 
    worktime = c(9,9+40/5),
    dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 9
    eps.A = 10,  
    rho.E = .9, 
    worktime = c(9,9+40/4),
    dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 10
    eps.A = 10,  
    rho.E = .9, 
    worktime = c(9,9+40/3),
    dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 11
    eps.A = 10,  
    rho.E = .9, 
    worktime = c(9,9+40/6),
    dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  ),list(#simulation 18
    eps.A = 10,  
    rho.E = .9, 
    worktime = c(9,9+40/7),
    dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'

  ),list(#simulation 19
      eps.A = 20,  
      rho.E = .3, 
      worktime = c(9,9+40/5),
      dayoff = c(0,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 2
      eps.A = 20,  
      rho.E = .3, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 3
      eps.A = 20,  
      rho.E = .3, 
      worktime = c(9,9+40/4),
      dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 4
      eps.A = 20,  
      rho.E = .3, 
      worktime = c(9,9+40/3),
      dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 5
      eps.A = 20,  
      rho.E = .3, 
      worktime = c(9,9+40/6),
      dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 24
      eps.A = 20,  
      rho.E = .3, 
      worktime = c(9,9+40/7),
      dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
      
    ),list(#simulation 25
      eps.A = 20,  
      rho.E = .6, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 8
      eps.A = 20,  
      rho.E = .6, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 9
      eps.A = 20,  
      rho.E = .6, 
      worktime = c(9,9+40/4),
      dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 10
      eps.A = 20,  
      rho.E = .6, 
      worktime = c(9,9+40/3),
      dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 11
      eps.A = 20,  
      rho.E = .6, 
      worktime = c(9,9+40/6),
      dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 30
      eps.A = 20,  
      rho.E = .6, 
      worktime = c(9,9+40/7),
      dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
      
    ),list(#simulation 31
      eps.A = 20,  
      rho.E = .9, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 8
      eps.A = 20,  
      rho.E = .9, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 9
      eps.A = 20,  
      rho.E = .9, 
      worktime = c(9,9+40/4),
      dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 10
      eps.A = 20,  
      rho.E = .9, 
      worktime = c(9,9+40/3),
      dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 11
      eps.A = 20,  
      rho.E = .9, 
      worktime = c(9,9+40/6),
      dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 36
      eps.A = 20,  
      rho.E = .9, 
      worktime = c(9,9+40/7),
      dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
      
    ),list(#simulation 37
      eps.A = 30,  
      rho.E = .3, 
      worktime = c(9,9+40/5),
      dayoff = c(0,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 2
      eps.A = 30,  
      rho.E = .3, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 3
      eps.A = 30,  
      rho.E = .3, 
      worktime = c(9,9+40/4),
      dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 4
      eps.A = 30,  
      rho.E = .3, 
      worktime = c(9,9+40/3),
      dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 5
      eps.A = 30,  
      rho.E = .3, 
      worktime = c(9,9+40/6),
      dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 42
      eps.A = 30,  
      rho.E = .3, 
      worktime = c(9,9+40/7),
      dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
      
    ),list(#simulation 43
      eps.A = 30,  
      rho.E = .6, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 8
      eps.A = 30,  
      rho.E = .6, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 9
      eps.A = 30,  
      rho.E = .6, 
      worktime = c(9,9+40/4),
      dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 10
      eps.A = 30,  
      rho.E = .6, 
      worktime = c(9,9+40/3),
      dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 11
      eps.A = 30,  
      rho.E = .6, 
      worktime = c(9,9+40/6),
      dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 48
      eps.A = 30,  
      rho.E = .6, 
      worktime = c(9,9+40/7),
      dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
      
    ),list(#simulation 49
      eps.A = 30,  
      rho.E = .9, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 8
      eps.A = 30,  
      rho.E = .9, 
      worktime = c(9,9+40/5),
      dayoff = c(0,3) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 9
      eps.A = 30,  
      rho.E = .9, 
      worktime = c(9,9+40/4),
      dayoff = c(0,3,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 10
      eps.A = 30,  
      rho.E = .9, 
      worktime = c(9,9+40/3),
      dayoff = c(0,3,5,6) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 11
      eps.A = 30,  
      rho.E = .9, 
      worktime = c(9,9+40/6),
      dayoff = c(0) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    ),list(#simulation 54
      eps.A = 30,  
      rho.E = .9, 
      worktime = c(9,9+40/7),
      dayoff = c(NULL) #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
      
  ))


######## 0.2 DEFINE DEFAULT SETTINGS ########

#****USER INSTR.: SET HERE THE PARAMETERS AS THEY OCCUR FOR ALL SIMULATIONS****#

# General explaination of part 0.2:
# this defines the default settings for each simulation. Variables that have been defined above in 'paramsets' 
# overwrite specific default values. Those variables not defined in 'paramsets' will receive their default value

defaults <- list( #initialize list

    ## simulation settings
    fs = 2, #samples per hour
    days = 200, #simulated number of days per person
    people = 5000, #number of persons simulated
    integration = "RK4", #numerical integration method used

    #night impulses
    night.f.mean = 8, #average number of night impulses 
    night.dshape = 1, #shape factor of the gamma distribution of impulse frequencies
    waketime = 7,

    #work stressors
    distribution = c("gamma","uniform")[2], #distribution type describing the variation in the number of work impulses 
    #between individuals (gamma being representative of the population; a uniform distribution representing an equal amount of people for each number of work impulses)
    work.f.mean = 6, #average number of work impulses per 8h, only relevant when the gamma distribution is chosen
    work.dshape = 3, #shape factor of the gamma distribution of impulse frequencies - only relevant for gamma
    work.unif = c(1,50), #range of work impulses - only relevant for uniform distribution
    worktime = c(9,17),
    anticipation = c("YES","NO")[1], #whether or not anticipation of work day stress is assumed
    dayoff = c(0,6), #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'

    ## physiological parameters
    kHPA = 2.194, #multiplication constant describing the cortisol release in response to a neural impulse 
    delay = .5, #delay between neural impulse and the cortisol release
    halflifeM = 80/60, #average half life of cortisol is 80 minutes or 4/3hour
    halflifeSD = 8/60, #between person variation
    eps.A = 20, #threshold of elastic load becoming allostatic load 
    rho.E = .6, #the elasticity constant - rate of recovery from elastic load
    eps.D = 400, #threshold of allostatic load causing disease

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
winDir <- "C:/Users/b8058356/OneDrive - Newcastle University/Courses and learning/BS28A - Major project/simulations/simoutput/samples/"
macDir <- "/Users/remcobenthemdegrave/OneDrive - Newcastle University/Courses and learning/BS28A - Major project/simulations"
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



##### 3 FUNCTIONS #####

## 3.1 TOP LEVEL FUNCTION / EXECUTION FUNCTION ###
exesim <- function(paramset,defaults,subDir){
  #this function executes the individual simulations (the number of simulations that is performed is defined by the number of lists in 'paramsets')
  #'paramset' is one of the lists within 'paramsets' 
  #'defaults' are the default settings as defined in 0.2
  #'subDir' is the lower level folder where the data from a simulation is stored
  
  #### 3.1.1 Initialization steps ####
  
  ## time stamp the start of the script execution, so we can determine how long the whole thing took
  start.time <- Sys.time()
  
  n <- 1 #seed counter; everytime a random number is drawn, the next random number is drawn using a new seed, but every execution of this script causes the same random numbers to be used
  
  ######## 3.1.2 DEFINE SETTINGS FROM PROVIDED PARAMETERS AND DEFAUTLS ########
  
  #define the settings - this line of code takes the default settings and overwrites those that were defined in 'paramset'
  settings <- sapply(names(defaults), function(param) ifelse(param %in% names(paramset),paramset[param],defaults[param]))
  
  #### 3.1.3 Create individuals ####
  
  #create data frame to store all information about individuals and create a row for each individual
  simpopulation <- data.frame(simperson = seq(settings$people) )
  
  ## define an individual cortisol decay value and add it to the data frame
  set.seed(n); n=n+1; simpopulation$halflife_pp <- abs( rnorm(settings$people, mean = settings$halflifeM, sd = settings$halflifeSD) )
  
  ## define the individual number of night impulses (indivual mean) from a gamma distribution and add it to the data frame
  set.seed(n); n=n+1; 
  simpopulation$night.f.mean_pp <- rgamma(settings$people, shape = settings$night.dshape, rate = settings$night.dshape/settings$night.f.mean)
  
  ## define the individual number of work impulses (indivual mean) from a gamma (or uniform) distribution and add it to the data frame
  set.seed(n); n=n+1; 
  simpopulation$work.f.mean_pp <- switch(settings$distribution, 
                                      gamma = rgamma(settings$people, shape = settings$work.dshape, rate = settings$work.dshape/settings$work.f.mean), 
                                      uniform = runif(settings$people, min = settings$work.unif[1], max = settings$work.unif[2]))

  n=n+3; 

  #### 3.1.4 Prepare for file saving ####
  
  #create a data frame for logging information about the simulation
  simlog <- list(preptime=NULL, simtime = NULL, filemergetime = NULL) 
  
  #randomly select 5 people of whom the first 7 consequent days will be plotted
  settings$samplepeople <- sample(simpopulation$simperson, ifelse(settings$people >= 5,5,settings$people), replace = FALSE) 
  
  #determine for which days in the simulation data is stored for graphical representation
  if(settings$days<=settings$savedays){settings$savedays <- 1:settings$days
  } else {settings$savedays <- seq( as.integer(settings$days/settings$savedays), settings$days, as.integer(settings$days/settings$savedays) )}
  
  #create folders for storing data
  dir.create(subDir)
  dir.create(file.path(subDir, 'Qdata'))
  dir.create(file.path(subDir, 'daydata'))
  dir.create(file.path(subDir, 'plotsample'))
  
  #this end the preparation, the duration of the preparation is saved
  simlog$preptime <-  difftime(Sys.time(), start.time, unit="mins")
  
  #### 3.1.5 Simulate #### 
  
  start.time <- Sys.time() #initiate new time recording to determine the duration of the simulation
  
  #####-----DELETE THIS -------#######
  #determine the scaling factor for translating (reference to the average peak at 9mmol/L in Miller)
  #if(is.na(settings$kHPA)){settings$kHPA<-8/mean(ldply(seq(settings$people), function(id) getscale(settings = settings, simperson = simpopulation[id,], n = n+id))[,1])}
  #####-----DELETE THIS -------#######
  
  #perform the simulation per person, while returning several evaluation values
  evalv <- sapply(seq(settings$people), function(id) ppsim(settings = settings, simperson = simpopulation[id,], n = n+id, subDir = subDir))
  
  simlog$recovery<-mean(unlist(evalv[1,])) # ratio of days on which 'full' recovery occurs
  simlog$aloadratio<-mean(unlist(evalv[2,])) # ratio people developing allostatic load
  simlog$sickratio<-mean(unlist(evalv[3,])) # ratio people developing disease
  
  simlog$simtime <- difftime(Sys.time(), start.time, unit="mins") #this end the simulation, register the time  
  
  #### 3.1.6 MERGE AND STORE DATA ####
  
  start.time <- Sys.time() #start a new time recording to determine how long file merging takes
  
  #merge selected information from simulated individuals into a new file and delete the originals
  bindfiles(file.path(getwd(), subDir, "Qdata"))
  bindfiles(file.path(getwd(), subDir, "daydata"))
  bindfiles(file.path(getwd(), subDir, "plotsample"))
  
  simlog$filemergetime <- difftime(Sys.time(), start.time, unit="mins") #this end the merge, register the time
  
  write.csv(simpopulation, file.path(getwd(), subDir, "simpopulation.csv")) #store characteristics of the individuals
  save(settings, file = file.path(getwd(), subDir, "settings.RData")) #store settings
  save(simlog, file = file.path(getwd(), subDir, "simlog.RData"))  #store the simulation log
  
  #### 3.1.7 MAKE CURVES (SINGLE SIMULATION) ####
  
  #if only one simulation is done, this line of code call upon the markdown script to create graphics 
  if(!grepl("comp",subDir, fixed = T)){  render("makegraphs.Rmd",output_file = paste0(subDir,".html"),params = list(datadir = subDir, paramsets = NULL, totalsimtime=difftime(Sys.time(), start0, unit="mins")))}
} 



## 3.2 SIMULATION INDIVIDUALS ####
ppsim <- function(settings, simperson, n, subDir){   
  # performs the simulation per individual (as defined in 'simpopulation') - called by 'exesim()'
  # simperson is the row of data in 'simpopulation' referring to a specific individual
  
  ## 3.2.1 create individual HPA activity (define per individual the impulses to the HPA axis)
  
  #determine the number of impulses that the simulated individual experiences on a particular day using draws from a poisson distribution
  set.seed(n); night.f <- rpois(settings$days, simperson$night.f.mean_pp)
  set.seed(n+1); work.f <- rpois(settings$days, simperson$work.f.mean_pp)
  
  #based on night.f and work.f determine the HPA activity of a person over the whole time course of the simulation
  HPA <- stack(as.data.frame(sapply(seq(settings$days), function(day) makeHPA(day, night.f[day], work.f[day], settings,n=n+simperson$simperson*day))))[,1]
  
  ## 3.2.2 calculate cortisol, elastic load and allostatic load through numeric integration
  data <- CEA(HPA, settings$fs, settings$delay, simperson$halflife_pp, settings$eps.A, settings$rho.E, settings$integration, settings$kHPA)
  
  ## 3.2.3 Store results
  
  savedays <- settings$savedays #for practical use
  
  #store data for making a graph of ODD's ratio of becoming sick 
  write.csv(data.frame(simperson=rep(simperson$simperson,length(savedays)),days=savedays,work.f.mean=rep(simperson$work.f.mean_pp,length(savedays)),Aload=data$A[savedays*settings$fs*24],Sick=data$A[savedays*settings$fs*24]>settings$eps.D)
            ,file.path(getwd(), subDir, "daydata",paste0("daydata",simperson$simperson)))
  
  # store data for making a graph of the average cortisol profiles
  write.csv(
    cbind(time=seq(0,24,1/settings$fs),data[1:(settings$fs*24+1),]),
    file.path(getwd(), subDir, "Qdata",paste0("Qdata",simperson$simperson)))
  
  # store data to create graphs of sample profiles of HPA, C, E, A and D
  if(simperson$simperson %in% settings$samplepeople){write.csv(
    cbind(simperson=rep(simperson$simperson,settings$fs*24*7+1),
          time=seq(0,24*7,1/settings$fs),
          HPA=HPA[1:(settings$fs*24*7+1)],data[1:(settings$fs*24*7+1),]),
    file.path(getwd(), subDir, "plotsample",paste0("plotsample",simperson$simperson))
  )}
  
  ## 3.2.4 Calculate evaluation values
  
  #count the ratio of days where there was dynamic load remaining at the end of the day
  recovery <- sum(data$L[seq(settings$days)*settings$fs*24] < rep(max(data$L)/20,settings$days))/settings$days
  
  #whether or not the person developed allostatic load during the simulation
  aloadYN <- data$A[nrow(data)]>0
  
  sickYN <- data$A[length(data$A)]>settings$eps.D
  
  return(data.frame(recovery,aloadYN,sickYN))
}



#### 3.3 CALCULATE CORTISOL & ALLOSTATIC LOAD ####
CEA <- function(HPA, fs, delay, halflife, eps.A, rho.E, integration = c("Euler", "RK4")[1], kHPA=1){
  #this function calculated cortisol, elastic load and allostatic load by integrating over time
  #called by ppsim()
  #HPA is the vector describing the neural impulses to the HPA axis
  #fs is the time interval between samples (sample frequency)
  #halflife is the cortisol half-life
  
  ## 3.3.1 preparation

  #before we start we convert the halflife value into a decay rate: 
  #standard decay function: dY/dt = Y*r (decay is a function of the current concentration)
  #one can rewrite this as: dY/Y = r*dt. Integrating over time gives ln(Y(t))-ln(Y(0)) = r*t or Y(t)/Y(0) = exp(r*t) or Y(t) = Y(0)*exp(r*t)
  #filling in: Y(t=1/2) = 1/2*Y(0) = Y(0)*exp(r*halflife) or -ln(2) = r*halflife or r = -ln(2)/halflife
  r <- -log(2)/halflife
  
  #introduce the delay between the impulses to the HPA axis and its cortisol output (for practical reasons, we do this outside the formula below)
  HPAout <- c(rep(0,round(delay*fs)),HPA[1:(length(HPA)-round(delay*fs))])*kHPA #
  
  #initiate and add values start values for the integration
  N <- length(HPAout)
  C <- L <- A <- as.numeric(c(0, rep(NA,N-1))) #initialize
  
  ## 3.3.2 performing the numerical integration
  
  #calculate cortisol values
  switch(integration,
         "Euler" =    # Euler numerical integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- Euler(Y1 = C[t], func = lingrowth, param = c(r, HPAout[t]), h = 1/fs)),
         "RK4" =   # Runge-Kutta 4th Order Integration
           sapply(seq_along(C), function(t) C[[t+1]] <<- RK4(Y1 = C[t], func = lingrowth, param = c(r, HPAout[t]), h = 1/fs)))
  switch(integration,
           "Euler" =    # Euler numerical integration
             sapply(seq_along(L), function(t) L[[t+1]] <<- Euler(Y1 = L[t], func = lingrowth, param = c(-rho.E, C[t]), h = 1/fs)),
           "RK4" =   # Runge-Kutta 4th Order Integration
             sapply(seq_along(L), function(t) L[[t+1]] <<- RK4(Y1 = L[t], func = lingrowth, param = c(-rho.E, C[t]), h = 1/fs)))
  sapply(seq_along(A), function(t) A[[t+1]] <<- ifelse(L[t] > eps.A, A[t] + (L[t] - eps.A)*1/fs, A[t]))
    
  return(data.frame(C,L,A))
}



#### 3.4 Euler integration ####
Euler <- function(Y1, func, param, h){
  #This function performs one iteration using the Euler Integration
  #called by CEA()
  
  # Iterative process
  Y2 <- Y1 + func(Y1, param)*h
  
  return(Y2)
}



#### 3.5 Runge-Kutta 4th Order Integration ####
RK4 <- function(Y1, func, param, h){
  #This function performs one iteration using the Runge-Kutta 4th Order Integration
  #called by CEA()
  
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
#Called by Euler() and RK4()
lingrowth <- function(Y1,param){
  r <- param[1]; input <- param[2]
  Y2 <- as.numeric(Y1)*r + as.numeric(input)
  return(Y2)
}



#### 3.7 create day HPA activity ####
makeHPA <- function(day, night.f, work.f, settings,n){
  #executed for each individual person, this function determines the HPA activity for a given day
  
  #as the value of f.work is based on an 8hour workday, adapt this value for the true number of worked hours 
  workdayfactor <- (settings$worktime[[2]]-settings$worktime[[1]])/8 
  
  #determine the time points of HPA activity
  set.seed(n)
  if((day - floor(day/7)*7) %in% settings$dayoff) #for non-working days there is only night time HPA activity
  {spikes <- abs(-rexp(night.f) + with(settings, waketime))
  }else{ #for working days, there is also HPA activity during work time
    spikes <- switch(settings$anticipation,
                     NO = c(abs(-rexp(night.f) + settings$waketime),
                            runif(work.f*workdayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])),
                     YES = c(abs(-rexp(night.f+work.f*workdayfactor) + settings$waketime),
                             runif(work.f*workdayfactor,min = settings$worktime[[1]], max = settings$worktime[[2]])))
  }
  
  spikes <- ceiling(spikes*settings$fs) #corresponding sample moments of the spikes (and prevent spikes at time point 0)
  
  #sum spikes that happen at the same sample moment
  HPA<-rep(0, settings$fs*24)
  HPA[unique(spikes)]<-sapply(unique(spikes), function(spike) sum(spikes==spike))
  
  return(HPA)
}



#### 3.8 combine files into a larger file #### 
bindfiles = function(mypath){
  #called by exesim()
  
  filenames=list.files(path=mypath, full.names=TRUE) #identify the files
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)}) #read the files into a list of data frames
  write.csv(
    Reduce(function(x,y) {rbind(x,y)}, datalist), #combine the data frames
    paste0(mypath,".csv")) #store the new file
  unlink(mypath, recursive = T, force = T) #delete folder with originals
}


#### 4 RUN THE SIMULATION FOR THE GIVEN PARAMETER SETS ####

start.time<-Sys.time() #start a new time recording to calculate how long the simulations take

sapply(seq(length(paramsets)), function(i) exesim(paramsets[[i]],defaults,datadir[i]))

difftime(Sys.time(), start.time, unit="mins")

totalsimtime=round(as.numeric(difftime(Sys.time(), start0, unit="mins"))) #in minutes

#### MAKE CURVES (COMPARISON) ####

#if it concerns a comparison, store comparison information and render comparison graph
if(length(datadir)>1){  
  save(paramsets, file = file.path(compdir, "paramsets.RData"))
  save(datadir, file = file.path(compdir, "datadir.RData"))
  render("makegraphs.Rmd",output_file = paste0(compdir,".html"),params = list(datadir = datadir, paramsets=paramsets, totalsimtime=totalsimtime))
  }
