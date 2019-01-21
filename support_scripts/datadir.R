#simlog<-list()

winDir <- "D:/Surfdrive/BS28A - Major project/simulations"
macDir <- "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations"
if (file.exists(winDir)){setwd(winDir)}; if (file.exists(macDir)){setwd(macDir)}; 
if (!file.exists("simoutput")){setwd(tk_choose.dir(getwd(), "Choose folder for storing simulation data"))
  if (!file.exists("simoutput")){dir.create("simoutput")}} 

simlog$maindir <- getwd()
simlog$subdir1 <- "simoutput"
simlog$subdir2 <- "comparison 1"
simlog$subdir3 <- c("simulation 2018-05-22_18-22",
"simulation 2018-05-22_18-23",
"simulation 2018-05-22_18-24",
"simulation 2018-05-22_18-25",
"simulation 2018-05-22_18-26",
"simulation 2018-05-22_18-28",
"simulation 2018-05-22_18-29",
"simulation 2018-05-22_18-30",
"simulation 2018-05-22_18-31")

i=1
setwd(with(simlog, file.path(maindir,subdir1,subdir2,subdir3[i])))
setwd(with(simlog, file.path(maindir)))

datadir<-with(simlog, sapply(subdir3, function(i) file.path(maindir,subdir1,subdir2,i)))

#sapply(datadir, function(i) print(i))
           
load("/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/simulation 2018-05-23_14-56/simlog.RData")

a<-as.data.frame(unlist(simlog$recovery))
a[,2]<-round(runif(500))

colMeans(as.data.frame(simlog$recovery))
colMeans(a)

load("/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/simulation 2018-05-23_15-20/simlog.RData")
mean(unlist(simlog$recovery[1,]))
mean(unlist(simlog$recovery[2,]))

load("/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/simulation 2018-05-23_15-27/simlog.RData")


row.names(paramsets)<-paste("sim",seq(nrow(paramsets)))

