#merge simulation results


winDir <- "\\\\Cnas.ru.nl/u276198/Surfdrive/BS28A - Major project/simulations/output/"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/output/"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 

# library(tcltk)
# dir1<-list.dirs(path = tk_choose.dir(getwd(), "Show folder contents of 1st simulation (set) to merge"), full.names = TRUE, recursive = TRUE)[-1]
# dir2<-list.dirs(path = tk_choose.dir(getwd(), "Show folder contents of 2nd simulation (set) to merge"), full.names = TRUE, recursive = TRUE)[-1]


getdatadir <- function(filepath){load(filepath);return(datadir)}

#**** USER INSTRUCTION: select the 'datadir.Rdata' file that contains the 1st (set of) path(s)  ****#
dir1<-getdatadir(filepath = file.choose())
#**** USER INSTRUCTION: select the 'datadir.Rdata' file that contains the 2nd (set of) path(s)  ****#
dir2<-getdatadir(filepath = file.choose())

datadir<-c(dir1,dir2)

#**** USER INSTRUCTION: or provide manual input  ****#
datadir<-c(paste0("/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/comp2018-05-30.14h02/sim",1:9))


#get settings for each simulation
set <- function(ddir){load(file.path(ddir,"settings.Rdata"));return(settings)}
settingsset <- lapply(datadir, function(ddir) set(ddir))

#test which settings are different between 2 simulations
b<-unique(unlist(lapply(seq_along(settingsset), function(i)
  which(lapply(names(settingsset[[1]]), function(name) unique(sort(settingsset[[1]][[name]]==settingsset[[i]][[name]]),decreasing = F)[1])==F)
)))

#create the merged paramsets
paramsets=lapply(settingsset, function(sett) sett[b])
#paramsets=as.data.frame(settingsset[[1]][b])
#paramsets=as.data.frame(t(matrix(unlist(lapply(settingsset, function(sett) sett[b])),nrow = length(b)))); names(paramsets)=names(settingsset[[1]][b])

#define data storage locations
filename<-paste0("merge", format(Sys.time(),"%Y-%m-%d.%Hh%M"))
mergedir<-file.path("simoutput",filename)
dir.create(mergedir)

save(paramsets, file = file.path(mergedir, "paramsets.RData"))
save(datadir, file = file.path(mergedir, "datadir.RData"))
render("makegraphs.Rmd",output_file = paste0(getwd(),"/simoutput/",filename,".html"),params = list(datadir = datadir, paramsets=paramsets))

setwd("/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/")
# 
# 
# choicedir<-tk_choose.dir(getwd(), "Show folder with data to publish")
# if(file.exists(file.path(choicedir,"datadir.RData"))){load(file.path(choicedir,"datadir.RData"));load(file.path(choicedir,"paramsets.Rdata"))}else{datadir<-choicedir; paramsets<-NULL}