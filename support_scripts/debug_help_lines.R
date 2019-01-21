paramset=paramsets[[1]]
subDir=datadir[1]
simperson = simpopulation[1,]
day =1
ddir=datadir[1]



head(Qdata)

head(daydata)



unique(plotsample$type)
p

12/.45

load('/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/comp2018-06-08.23h56/datadir.Rdata')
load('/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/profiles2018-06-06.23h58/datadir.Rdata')

load('/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/simoutput/keep these/comp2018-05-31.22h21/paramsets.RData')


maindir = "/Users/remcobenthemdegrave/surfdrive/BS28A - Major project/simulations/"
ddir1 = "simoutput/profiles2018-06-08.23h09/sim2" #not correlated
ddir2 = "simoutput/profiles2018-06-06.23h58/sim2" #correlated

datadir = paste0(maindir, c(ddir1,ddir2))


ddir=datadir[1]

load(file.path(ddir,"simpopulation.Rdata"))

library(ggplot2)
ggplot(simpopulation, aes(y=stressfM_pp, x=wakefM_pp)) + geom_point() + expand_limits(x=c(0,90),y=c(0,35)) + scale_x_continuous(breaks=seq(0, 90, 10)) +
labs(x='night impulses', y='work stressors') + theme_bw()


ddir=datadir[2]


load(file.path(ddir,"simpopulation.Rdata"))

ggplot(simpopulation, aes(y=stressfM_pp, x=wakefM_pp+stressfM_pp)) + geom_point() + expand_limits(x=c(0,90),y=c(0,35)) + scale_x_continuous(breaks=seq(0, 90, 10)) +
  labs(x='night impulses', y='work stressors') + theme_bw()




ddir1 = 'simoutput/comp2018-06-07.19h43/sim1'
ddir = paste0(maindir, c(ddir1))


#read data
plotsample <- read.csv(file.path(ddir,"plotsample.csv"))[,-c(1,2)]

#reshape into long format
plotsample=gather(plotsample, key = "type", value = "value", -c(time,simperson))

#make sure that the order is kept and asign labels
plotsample$type <- factor(plotsample$type, levels=c("HPA", "C", "L","A"), labels=c("impulses", "cortisol", "dynamic load","allostatic load"), ordered=TRUE)

#plot

a = subset(plotsample, simperson = unique(simperson)[1])
a = subset(a, time <= 24)

ggplot(a, aes(x=time, y=HPA)) +
  geom_line() +
  ylab("number of impulses") + xlab("time [h]") + theme_bw() 

#plot
ggplot(plotsample, aes(x=time, y=value)) +
  geom_line() + facet_grid(type ~ simperson, scales="free") +
  ylab("value [arbitrary units]") + xlab("time [h]") + theme_bw() 


i=list.dirs(choicedir)[2]

##define datadir

#list the file paths to the different simulations in the folder (using a bit of a complex way to get there)
subdirlist <- unlist(lapply(list.dirs(choicedir), function(i) unlist(strsplit(i,"/"))[length(unlist(strsplit(i,"/")))])) #for each folder path abstract the name of the last folder in the folderpath
subdirlist <- subdirlist[subdirlist %in% dir(choicedir)] #keep only those that are a direct subpath of the choicedir
subdirlist <- subdirlist[grepl('sim', subdirlist)] #test for each subdir if it contains 'sim' and keep only those
datadir <- file.path(choicedir,subdirlist)
if(length(datadir)==0){datadir<-choicedir} #if the choicefolder already contains the simulation data, then keep this folder

##define paramsets

#get settings for each simulation
set <- function(ddir){load(file.path(ddir,"settings.Rdata"));return(settings)}
settingsset <- lapply(datadir, function(ddir) set(ddir))

#test which settings are different between 2 simulations
uparam<-unique(unlist(lapply(seq_along(settingsset), function(i)
  which(lapply(names(settingsset[[1]]), function(name) unique(sort(settingsset[[1]][[name]]==settingsset[[i]][[name]]),decreasing = F)[1])==F)
)))

#create the merged paramsets
if(length(uparam)==0){paramsets<-NULL}else{paramsets=lapply(settingsset, function(sett) sett[uparam])}




#scale stressors
daydataall$stressfM <- (scale(daydataall$stressfM))[,1]

ggplot(subset(daydataall, days == max(days) & A_D == 400), aes(x=stressfM, y=Aload)) +
  geom_point() + geom_smooth(method="lm", se=F, colour="red") +
  facet_grid(rL ~ L_A, scales="free") + 
  labs(x = "individual average work stressors (standardized)", y = "Allostatic load") +
  theme_minimal()

ggsave('allos_no_ant.png', width = 12.5, height = 10, units = "cm")




#plot odds ratio's
ggplot(odds_all, aes(y= OR, x = day, shape = A_D)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  facet_grid(L_A ~ rL, labeller = label_bquote(rows = rho[E] == .(.6), cols = epsilon['A'] == .(20))) +
  geom_hline(yintercept = 1, linetype=2) +
  #coord_flip() +
  coord_trans(y="log10") +
  #coord_cartesian(ylim=c(.01, 100)) + #limit the visible range
  scale_y_continuous(breaks=c(1, 10, 100),minor_breaks=c(seq(1,10,1),seq(20,100,10))) + # Ticks 
  scale_x_continuous(breaks=odds_all$day) +
  scale_shape_discrete(name=expression(epsilon['D'])) +
  #labs(title = expression('Relative disease risk - split by value of '*epsilon['A']*' (horizontal) and '*rho['E']*' (vertical)'), x = "Day", y = "Relative disease risk") +
  labs(x = "Day", y = "Relative disease risk") +
  theme_minimal() + theme(legend.position="bottom")

ggsave('odds_ant.png', width = 12.5, height = 10, units = "cm")



write.csv(simlogs, 'table.csv')



sampleplots <- function(ddir){
  
  #load settings
  load(file.path(ddir,"settings.Rdata"))
  
  #read data
  plotsample <- read.csv(file.path(ddir,"plotsample.csv"))[,-c(1,2)]
  plotsample$D <- plotsample$A>settings$thAM
  
  #reshape into long format
  plotsample=gather(plotsample, key = "type", value = "value", -c(time,simperson))

  #make sure that the order is kept and asign labels
  plotsample$type <- factor(plotsample$type, levels=c("HPA", "C", "L","A","D"), labels=c("Impulses", "Cortisol", "Elastic load","Allostatic load","Diseased"), ordered=TRUE)
  
  #plot
  ggplot(plotsample, aes(x=time, y=value)) +
    geom_line() + facet_grid(type ~ simperson, scales="free") +
    ylab("Value [arbitrary units]") + xlab("time [hours]") + theme_minimal() 
}  

sampleplots(datadir[5])

ggsave('example_no_ant.png', width = 25, height = 20, units = "cm")







sampleplots <- function(ddir){
  
  #load settings
  load(file.path(ddir,"settings.Rdata"))
  
  #read data
  plotsample <- read.csv(file.path(ddir,"plotsample.csv"))[,-c(1,2)]
  plotsample$D <- plotsample$A>settings$thAM
  
  #reshape into long format
  plotsample=gather(plotsample, key = "type", value = "value", -c(time,simperson))
  
  #make sure that the order is kept and asign labels
  plotsample$type <- factor(plotsample$type, levels=c("HPA", "C", "L","A","D"), labels=c("Impulses", "Cortisol", "Elastic load","Allostatic load","Diseased"), ordered=TRUE)
  
  #plot
  ggplot(plotsample, aes(x=time, y=value)) +
    geom_line() + facet_grid(type ~ simperson, scales="free") +
    ylab("Value [arbitrary units]") + xlab("time [hours]") + theme_minimal() 
}  

sampleplots(datadir[5])



#read data that contains impulse frequencies per person
simpopulation <- read.csv(file.path(ddir,"simpopulation.csv"))[,-c(1,2,3,6,7,8,9)]

#plot correlation between work impulses and night impulses
ggplot(simpopulation, aes(x=wakefM_pp,y=stressfM_pp)) + geom_point() + 
  labs(x="night impulses",y="work impulses") + theme_minimal()


simpopulation <- gather(simpopulation,
                        key = impulsetype, value = freq)
simpopulation$impulsetype <- factor(simpopulation$impulsetype, levels = c("wakefM_pp","stressfM_pp"),
                                    labels = c("night impulses","work stressors"))

Mw = round(mean(simpopulation$freq[which(simpopulation$impulsetype=="night impulses")]), digits = 1)
Ms = round(mean(simpopulation$freq[which(simpopulation$impulsetype=="work stressors")]), digits = 1)
SDw = round(sd(simpopulation$freq[which(simpopulation$impulsetype=="night impulses")]), digits = 1)
SDs = round(sd(simpopulation$freq[which(simpopulation$impulsetype=="work stressors")]), digits = 1)
Medw = round(median(simpopulation$freq[which(simpopulation$impulsetype=="night impulses")]), digits = 1)
Meds = round(median(simpopulation$freq[which(simpopulation$impulsetype=="work stressors")]), digits = 1)

#plot impulse distributions
ggplot(simpopulation, aes(freq,color=impulsetype,fill=impulsetype)) + geom_density(alpha=0.55) + 
  labs(title=paste0("Median(SD), wake: ",Medw,"(",SDw,")"," stress: ",Meds,"(",SDs,")"),x="impulse frequency") + 
  labs(x="impulse frequency") + 
  theme_minimal()  + theme(legend.position=c(.8, .8))

ggsave('corr_ant_part2.png', width = 12.5, height = 10, units = "cm")



#plot the distributions of allostatic load
ggplot(subset(daydataall, days == max(days) & A_D == 400), aes(x=scale(stressfM)[,1], y=Aload)) +
  geom_point() + geom_smooth(method="lm", se=F, colour="red") +
  facet_grid(rL ~ L_A, scales="free") + 
  #labs(title = expression('Allostatic load per individual at day 200 - split by value of L'['A']*' (horizontal) and r'['L']*' (vertical)'), x = "average frequency of stressors", y = "Allostatic load") +
  labs(x = "individual average work stressors (standardized)", y = "Allostatic load") +
  theme_minimal()

ggsave('allostatic_ant_part2.png', width = 12.5, height = 10, units = "cm")




ggplot(subset(daydataall, A_D == 400), aes(x=days, y=Aload, color = simperson)) +
  geom_line() +   facet_grid(rL ~ L_A, scales="free") + 
  labs(x = 'time [days]', y = "Allostatic load") + theme_minimal() + theme(legend.position="none") 
  #labs(title = 'allostatic load growth curves of 100 people becoming diseases', x = 'time [days]', y = "Allostatic load")

ggsave('time profiles.png', width = 12.5, height = 10, units = "cm")


tempdir <- datadir
datadir <- tempdir

datadir <- datadir[c(1,12,23,34,45,51)]
paramsets <- paramsets[c(1,12,23,34,45,51)]

a<-datadir[c(1,12,23,34,45,51)]



a<-getdatadir(file.path(maindir,"set2"))

rm(datadir)



#function that creates a density plot of the wake and stressor frequencies
load(file.path(ddir,"simpopulation.Rdata"))
load(file.path(ddir,"settings.Rdata"))

#plot correlation between work impulses and night impulses
ggplot(simpopulation, aes(y=stressfM_pp,x=stressfM_pp+wakefM_pp)) + geom_point() + 
  labs(x="night impulses",y="work impulses") + theme_minimal() + coord_cartesian(xlim=c(0,61) ,ylim=c(0, 31))

ggsave('correlation_impulses.png', width = 12.5, height = 10, units = "cm")

## prepare for making density plot

if(settings$anticipation=='YES'){
  simpopulation$wakefM_pp <- simpopulation$wakefM_pp + simpopulation$stressfM_pp
  simpopulation <- gather(simpopulation[,-c(1,2)],key = impulsetype, value = freq)
  simpopulation$impulsetype <- factor(simpopulation$impulsetype, 
                                      levels = c("wakefM_pp","stressfM_pp"),
                                      labels = c("night impulses","work impulses"))
}else{simpopulation <- gather(simpopulation[,c(3,4)],key = impulsetype, value = freq)
simpopulation$impulsetype <- factor(simpopulation$impulsetype, levels = c("wakefM_pp","stressfM_pp"),
                                    labels = c("night impulses","work impulses"))}

ggplot(simpopulation, aes(freq,color=impulsetype,fill=impulsetype)) + geom_density(alpha=0.55) + 
  labs(x="Impulse frequency") + 
  theme_minimal() + theme(plot.title = element_text(size=10), legend.position=c(.8, .85)) + coord_cartesian(xlim=c(0,60) ,ylim=c(0, .13))

ggsave('density_impulses.png', width = 12.5, height = 10, units = "cm")

pa
ggsave('cortisol_since_awake_20h_no_corr.png', width = 12.5, height = 9, units = "cm")

pn
ggsave('cortisol_night_corr.png', width = 12.5, height = 9, units = "cm")


multiplot(plotlist = graphs, cols = 2)


Qdata.pp <- subset(Qdata, time <= 24 & person == 5221)[c('time','HPA')]
a <- data.frame(time = seq(0,24,1/60),impulse=c(rep(Qdata.pp$HPA[1:nrow(Qdata.pp)-1],each=30),Qdata.pp$HPA[nrow(Qdata.pp)]))

resample <- function(data){return(data.frame(time = seq(0,24,1/60),person=rep(data$person,24*60+1),impulse=c(rep(data$HPA[1:nrow(data)-1],each=30),data$HPA[nrow(data)])))}
Qdata.rs <- ldply(sample(unique(Qdata$person),4), function(pp) resample(subset(Qdata, time <= 24 & person == pp)[c('time','person','HPA')]))

ggplot(Qdata.rs, aes(x=time,y=impulse)) + geom_line() + facet_wrap(~person, ncol=4, labeller = label_both) +
  scale_x_continuous(breaks=seq(0, 24, 4)) + scale_y_continuous(breaks=seq(0, 8, 2)) +
  ylab(expression("Impulses [#]")) + xlab("time [hours]") +
  theme_minimal() + theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())

ggsave('Impulses_example2.png', width = 25, height = 5, units = "cm")




distr <- data.frame(time=abs(-rexp(10000)+7),type=rep("exp",10000))
ggplot(distr, aes(distribution,color=type,fill=type)) + geom_density(alpha=0.55) + 
  labs(x="time [hours]") + 
  theme_minimal()  + theme(legend.position="none")

ggsave('Exp_example.png', width = 12.5, height = 9, units = "cm")



#make examples of distributions from poisson
getpois <- function(simperson){return(data.frame(person=rep(simperson$simperson,200),nnight=rpois(200,simperson$wakefM_pp),nwork=rpois(200,simperson$stressfM_pp)))}
distr <- ldply(sample(unique(simpopulation$simperson),4), function(pp) getpois(subset(simpopulation, simperson == pp)))
distr <- gather(distr, key = "type", value = "freq", -person)
distr$person<-factor(distr$person)
distr$type<-factor(distr$type, labels = c("night impulses","work impulses"))

ggplot(distr, aes(x=freq,fill=type)) + geom_histogram(alpha=0.65) + 
  facet_wrap(~person, ncol=4, labeller = label_both) + labs(x="Frequency of impulses") + 
  theme_minimal()  + theme(legend.position="top")

ggsave('pois_example.png', width = 25, height = 6, units = "cm")

