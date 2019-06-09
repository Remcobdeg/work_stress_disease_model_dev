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


changedaydata <- function(ddir){
  
  daydata <- read.csv(file.path(ddir,"daydata.csv"))[,-c(1)]
  names(daydata)[4] <- 'work.f.mean'
  write.csv(daydata, file.path(ddir,"daydata.csv"))

}

#get settings for each simulation
lapply(datadir, function(filepath) changedaydata(filepath))

