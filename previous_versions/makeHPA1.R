spiket <- function(spike,N){
  moments <- rep(0,settings$N)
  moments[spike] <- 1
  return(moments)}

makespikes <- function(day, simperson, settings,n){
  set.seed(n)
  if((day - floor(day/7)*7) %in% settings$dayoff) #for non-working days there is only night time HPA activity
  {set.seed(n); spikes <- abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR))
  }else{ #for working days, there is also HPA activity during work time
    spikes <- c(abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR)),
                runif(simperson$stressfM_pp,min = settings$worktime[[1]], max = settings$worktime[[2]]))}
  spikes <- round(spikes*settings$fs) #corresponding sample moments of the spikes
  return(spikes)
}

a<-rowSums(as.data.frame(sapply(makespikes(day, simperson, settings,n+id*day), function(spike) spiket(spike, settings$N))))

rowSums(as.data.frame(sapply(makespikes(day, simperson, settings,n), function(spike) spiket(spike, settings$N))))

# makeHPA <- function(day, simperson, settings,n){
#   set.seed(n)
#   if((day - floor(day/7)*7) %in% settings$dayoff) #for non-working days there is only night time HPA activity
#   {set.seed(n); spikes <- abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR))
#   }else{ #for working days, there is also HPA activity during work time
#     spikes <- c(abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR)),
#                 runif(simperson$stressfM_pp,min = settings$worktime[[1]], max = settings$worktime[[2]]))}
#   spikes <- round(spikes*settings$fs) #corresponding sample moments of the spikes
# 
#   HPA <- rowSums(as.data.frame(sapply(spikes, function(spike) spiket(spike,settings$N))))
#   return(HPA)
# }



id=1
simperson = simpopulation[id,]
start.time<-Sys.time()
sapply(settings$days, function(day) rowSums(as.data.frame(sapply(makespikes(day, simperson, settings,n+id*day), function(spike) spiket(spike, settings$N)))))
duration<-c(duration,Sys.time()-start.time)

rm(moments)
rm(day)
day=1
