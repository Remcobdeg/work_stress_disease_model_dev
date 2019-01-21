### this script changes the layout and names of an old settings version to the new style

winDir <- "\\\\Cnas.ru.nl/u276198/Surfdrive/BS28A - Major project/simulations/simoutput/"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 

#**** USER INSTRUCTION: select the settings file to be changed ****#
filepath = file.choose()

load(filepath)

#new settings as on 27-05-2018
newsettings <- list( #initialize list
  
  ## simulation settings
  fs = 10, #samples per hour
  N = NULL, #number of samples in one day - dummy, this value is defined later
  days = 10,
  people = 500,
  integration = "RK4",
  
  #awakening impulse
  wakefM = 20, 
  wakefSD = 20,
  waketime = 6,
  CAR = .5, #duration (in hours) of the cortisol awakening response
  nightgrowxp = 1,#exponential growth factor 
  
  #work stressors
  stressfM = 40, 
  stressfSD = 40, #frequency of stress stressors
  worktime = list(8,17),
  anticipation = c("YES","NO")[1],
  dayoff = list(0,6), #from 0-6: 'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'
  
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
  savedays = 5,
  samplepeople = list(1,2,3,4,5)
)

#**** USER MESSAGE: these variables have a new name ****#
names(settings)[!(names(settings) %in% names(newsettings))]

#**** USER MESSAGE: these are the new names (not ordered!) ****#
names(newsettings)[!(names(newsettings) %in% names(settings))]

#**** USER INSTRUCTION: assure that the new variable is added to the list, as the example below ****#
if(with(settings, exists("scaleC"))){settings$kHPA = settings$scaleC}
  
#rename the settings
settings = sapply(names(newsettings), function(name) settings[[name]]) #set both in the same order
names(settings)==names(newsettings) #rename

#change structure of 'offday' and 'worktime' if necessary
settings$worktime=unlist(settings$worktime)
settings$dayoff=unlist(settings$dayoff)

save(settings, file = filepath)


# newstyle=c(
#   "fs",
#   "N",            
#   "days",         
#   "people",       
#   "integration",  
#   "CAR"         ,
#   "wakefM"       ,
#   "wakefSD"      ,
#   "waketime"     ,
#   "nightgrowxp"  ,
#   "stressfM"     ,
#   "stressfSD"   ,
#   "worktime"     ,
#   "anticipation" ,
#   "dayoff"       ,
#   "kHPA"         ,
#   "halftimeM"    ,
#   "halftimeSD"  ,
#   "thLM"         ,
#   "thLSD"        ,
#   "rLM"          ,
#   "rLSD"         ,
#   "thAM"         ,
#   "thASD"       ,
#   "savedays"  
# )