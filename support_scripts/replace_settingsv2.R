### this script changes the layout and names of an old settings version to the new style

winDir <- "C:/Users/b8058356/OneDrive - Newcastle University/Courses and learning/BS28A - Major project/simulations/simoutput/samples/"
macDir <- "/Users/remcobenthemdegrave/OneDrive - Newcastle University/Courses and learning/BS28A - Major project/simulations"
if (file.exists(winDir)){setwd(winDir)}; 
if (file.exists(macDir)){setwd(macDir)}; 

#get the paths to the settings files to change
choicedir<-tk_choose.dir(getwd(), "Show folder with settings file(s) to change")

#list the file paths to the different simulations in the folder (using a bit of a complex way to get there)
subdirlist <- unlist(lapply(list.dirs(choicedir), function(i) 
  unlist(strsplit(i,"/"))[length(unlist(strsplit(i,"/")))])) #for each folder path...
#...abstract the name of the last folder in the folderpath
subdirlist <- subdirlist[subdirlist %in% dir(choicedir)] #keep only those that are a direct subpath of...
#...the choicedir
subdirlist <- subdirlist[grepl('sim', subdirlist)] #test for each subdir if it contains 'sim' ...
#...and keep only those
datadir <- file.path(choicedir,subdirlist)
if(length(datadir)==0){datadir<-choicedir} #if the choicefolder already contains the simulation data, ...
#...then keep this folder







changesettings <- function(filepath){
  
  load(paste0(filepath,"/settings.Rdata"))
  
  #new settings as on 27-05-2018
  newsettings <- list( #initialize list
    
    ## simulation settings
    fs = 10, #samples per hour
    days = 10,
    people = 500,
    integration = "RK4",
    
    #awakening impulse
    night.f.mean = 20, 
    night.dshape = 20,
    waketime = 6,
    
    #work stressors
    distribution = "uniform",
    work.f.mean = 40, 
    work.dshape = 3, 
    work.unif = c(1, 50),
    worktime = list(8,17),
    anticipation = c("YES","NO")[1],
    dayoff = list(0,6), #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
    
    ## physiological parameters
    kHPA = NA, #set to NA to have to code determine the scaling factor
    delay = .5,
    halflifeM = 80/60, #half time of cortisol is 80 minutes or 4/3hour
    halflifeSD = 8/60, #between person variation
    eps.A = 20, #threshold of dynamic load becoming allostatic load 
    rho.E = .6, #rate at which the body recovers form the dynamic load
    eps.D = 1000, #threshold of allostatic load causing disease
  
    ## number of days to save and display
    savedays = 5,
    samplepeople = list(1,2,3,4,5)
  )
  
  #**** USER MESSAGE: these variables have a new name ****#
  names(settings)[!(names(settings) %in% names(newsettings))]
  
  #**** USER MESSAGE: these are the new names (not ordered!) ****#
  names(newsettings)[!(names(newsettings) %in% names(settings))]
  differences = sum(!(names(newsettings) %in% names(settings)))
  
  
  #copy the settings values with the same name in newsettings to newsettings
  for (name in names(newsettings)[(names(newsettings) %in% names(settings))]){
    print(name)
    print(newsettings[name])
    print(settings[name])
    newsettings[name] = settings[name]
  }
  
  #**** USER INSTRUCTION: assure that the new variable is added to the list, as the example below ****#
  changed = 0
  if(with(settings, exists("scaleC"))){newsettings$kHPA = settings$scaleC; changed = changed +1}
  
  if(with(settings, exists("wakefM"))){newsettings$night.f.mean = settings$wakefM; changed = changed +1}
  if(with(settings, exists("dshape.wake"))){newsettings$night.dshape = settings$dshape.wake; changed = changed +1}
  if(with(settings, exists("stressfM"))){newsettings$work.f.mean = settings$stressfM; changed = changed +1}
  if(with(settings, exists("dshape.stress"))){newsettings$work.dshape = settings$dshape.stress; changed = changed +1}
  if(with(settings, exists("stressunif"))){newsettings$work.unif = settings$stressunif; changed = changed +1}
  if(with(settings, exists("CAR"))){newsettings$delay = settings$CAR; changed = changed +1}
  if(with(settings, exists("halftimeM"))){newsettings$halflifeM = settings$halftimeM; changed = changed +1}
  if(with(settings, exists("halftimeSD"))){newsettings$halflifeSD = settings$halftimeSD; changed = changed +1}
  if(with(settings, exists("thLM"))){newsettings$eps.A = settings$thLM; changed = changed +1}
  if(with(settings, exists("rLM"))){newsettings$rho.E = settings$rLM; changed = changed +1}
  if(with(settings, exists("thAM"))){newsettings$eps.D = settings$thAM; changed = changed +1}
  
  if(differences == changed){
    settings = newsettings
    save(settings, file = paste0(filepath,"/settings.Rdata"))
    print(paste(filepath, "\n",changed, "parameters changed."))
  } else {print("WARNING! Settings is not properly converted and therefore not saved!")}
}

#get settings for each simulation
lapply(datadir, function(filepath) changesettings(filepath))

