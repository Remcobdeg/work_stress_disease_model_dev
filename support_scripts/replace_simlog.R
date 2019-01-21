### this script changes the layout and names of an old settings version to the new style

winDir <- "\\\\Cnas.ru.nl/u276198/Surfdrive/BS28A - Major project/simulations/simoutput/"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 

#**** USER INSTRUCTION: select the simlog file to be changed ****#
filepath = file.choose()

load(filepath)

#new settings as on 27-05-2018

#do the same for simlog

newsimlog = list(
  preptime = NULL,
  simtime = NULL,
  filemergetime = NULL,
  kHPA = NULL,
  recovery = NULL,
  aloadratio = NULL,
  sickratio = NULL
)


#**** USER MESSAGE: these variables have a new name ****#
names(simlog)[!(names(simlog) %in% names(newsimlog))]

#**** USER MESSAGE: these are the new names (not ordered!) ****#
names(newsimlog)[!(names(newsimlog) %in% names(simlog))]

#load settings.Rdata that belongs to the simlog file
filepatha=unlist(strsplit(filepath, "/"))
load(paste(c(getwd(),filepatha[(which(filepatha=="simoutput")+1):(length(filepatha)-1)],"settings.Rdata"),collapse = '/'))

#**** USER INSTRUCTION: assure that the new variable is added to the list, as the example below ****#
if(with(settings, exists("scaleC"))){simlog$kHPA = settings$scaleC}
if(with(settings, exists("kHPA"))){simlog$kHPA = settings$kHPA}

#rename the settings
simlog = sapply(names(newsimlog), function(name) simlog[[name]]) #set both in the same order

simlog = lapply(names(newsimlog)[1:7], function(name) simlog[[name]]) #set both in the same order
names(simlog)<-names(newsimlog) #rename

save(simlog, file = filepath)

