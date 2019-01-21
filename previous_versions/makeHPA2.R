makeHPA <- function(day, simperson, settings,n){
  set.seed(n)
  if((day - floor(day/7)*7) %in% settings$dayoff) #for non-working days there is only night time HPA activity
  {set.seed(n); spikes <- abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR))
  }else{ #for working days, there is also HPA activity during work time
    spikes <- c(abs(-rexp(simperson$wakefM_pp, rate = settings$nightgrowxp) + with(settings, waketime + CAR)),
                runif(simperson$stressfM_pp,min = settings$worktime[[1]], max = settings$worktime[[2]]))}
  spikes <- round(spikes*settings$fs) #corresponding sample moments of the spikes
  
  HPA<-rep(0, settings$N)
  HPA[unique(spikes)] <- sapply(unique(spikes), function(spike) sum(spikes==spike))
  return(HPA)
}

HPA <- makeHPA(day, simperson, settings,n)
