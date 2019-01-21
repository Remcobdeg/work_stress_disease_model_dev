choicedir<-tk_choose.dir(getwd(), "Show folder with data to publish")
if(file.exists(file.path(choicedir,"datadir.RData"))){load(file.path(choicedir,"datadir.RData"));load(file.path(choicedir,"paramsets.Rdata"))}else{datadir<-choicedir; paramsets<-NULL}

ddir=datadir[1]

getsick <- unique(daydata$simperson[which(daydata$Sick==T & daydata$days == max(daydata$days))])
healthy <- unique(daydata$simperson[which(daydata$Sick==F & daydata$days == max(daydata$days))])

length(getsick)/length(healthy)


sim1=unique(daydataall$sim.nr)[1]
sim3=unique(daydataall$sim.nr)[3]
a=daydataall[which(daydataall$sim.nr %in% c(sim1,sim3)),]

sum(sapply(unique(a$sim.nr), function(nr) length(unique(a$sick[which(a$sim.nr==nr)]))<2))>0
sum(sapply(unique(a$sim.nr), function(nr) length(unique(a$sick[which(a$sim.nr==sim3)]))<2))<2

t(matrix(
  unlist(lapply(paramsets, function(pset) repl(pset))),nrow=length(names(paramsets[[1]]))
  ,dimnames = list(names(paramsets[[1]]),paste("sim",seq_along(paramsets)))))


as.data.frame(lapply(paramsets, function(pset) repl(pset))[[1]])
show=data.frame(matrix(nrow = length(paramsets),ncol = length(names(paramsets[[1]])),
                       dimnames = list(paste("sim",seq_along(paramsets)),names(paramsets[[1]]))))
str(show[,i])<-str(paramsets[[1]][i])
for(i in length(paramsets)){show[i,]<-repl(paramsets[[i]])}

i=1
i=2
i=3

paramssetsshow=repl(paramsets[[1]])
for(i in 2:length(paramsets)){paramssetsshow<-rbind(paramssetsshow,repl(paramsets[[i]]))}


a=as.data.frame(lapply(paramsets, function(pset) repl(pset))[[1]])
rbind(a,as.data.frame(lapply(paramsets, function(pset) repl(pset))[[2]]))


as.character(lapply(paramsets, function(pset) repl(pset)))

as.data.frame(t(lapply(paramsets, function(pset) repl(pset))),row.names = paste("sim",seq_along(paramsets)))

kable(lapply(paramsets, function(pset) repl(pset)))

