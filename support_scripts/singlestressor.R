#observe time range of 2h with fs = fs
#10000 people with one stressor at 1=0
#duration defined by rgamma of .5
#halftime given by script

library(lattice)

  
  fs = 60
  people = 1000
  halftime <- abs(rnorm(people, mean = 60/60, sd = 8/60))
  dur = abs(rnorm(people, mean = .5, sd = .2))
  #dur = rgamma(people, shape = 120, scale = .5/120) #rep(.5, people)
  #densityplot(rgamma(people, shape = 40, scale = .5/40))
  #densityplot(rnorm(people, mean = .5, sd = .2))
  
  
  spop=data.frame(person=seq(people),halftime,dur = dur)
  
  all<-ldply(spop$person, function(pp) sim(spop[pp,],fs))
  ggplot(all, aes(time,C,group=person)) + geom_line() 
    
  
  #calculate quantiles in data
  qq <- all %>% group_by(time) %>% summarise(
    Q05=quantile(C,probs=.10),
    Q25=quantile(C,probs=.25),
    Q50=quantile(C,probs=.50),
    Q75=quantile(C,probs=.75),
    Q95=quantile(C,probs=.90))
  
  #reshape to long data frame for ggplot
  qq <- gather(qq, key=quantile, value=C, -time)
  
  
  ggplot(qq, aes(time,C,group=quantile)) + geom_line(aes(linetype=quantile)) + 
    scale_linetype_manual(values=c("dotted","twodash","solid","twodash","dotted"),
                          name="Percentile",
                          breaks=c("Q90", "Q75", "Q50", "Q25", "Q10"),
                          labels=c("90th", "75th", "50th", "25th", "10th")) +
    labs(x="time [h]", y="simulated cortisol [arbitrary units]") +
    theme_minimal() + theme(legend.justification=c(1,1), legend.position=c(1,1), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
  
  





#function
sim <- function(sp,fs){
  I <- rep(0, 2*fs)
  I[1:(sp$dur*fs)] <- 1
  
  return(data.frame(person = rep(sp$person,fs*2), time=seq(2*fs), I=I, C=HPAA(I, fs, sp$halftime)))
}





plot(ts(HPAA(I, fs, halftime)))

#### 3.1 CALCULATE CORTISOL & ALLOSTATIC LOAD ####
HPAA <- function(I, fs, halftime){
  r <- -log(2)/halftime

  N <- length(I)
  C <- as.numeric(c(0, rep(NA,N-2))) #initialize
  
  #calculate cortisol values
           sapply(seq_along(C), function(t) C[[t+1]] <<- RK4(Y1 = as.numeric(C[t]), func = lingrowth, param = c(r, I[t]), h = 1/fs))
  return(C/max(C))
}

#### 3.3 Runge-Kutta 4th Order Integration ####
RK4 <- function(Y1, func, param, h){
  ## This function performs one iteration using the Runge-Kutta 4th Order Integration
  
  k1 <- func(Y1, param)
  k2 <- func(Y1 + 1/2*h*k1, param)
  k3 <- func(Y1 + 1/2*h*k2, param)
  k4 <- func(Y1 + h*k3, param)
  
  # Iterative process
  Y2 <- Y1 + (1/6)*h*(k1+2*k2+2*k3+k4)
  
  return(Y2)
}

#### 3.4 Linear growth (used in Euler and RK4 functions) ####
#This function performs one iteration of linear growth
lingrowth <- function(Y1,param){
  r <- param[1]; input <- param[2]
  Y2 <- as.numeric(Y1)*as.numeric(r) + as.numeric(input)
  return(as.numeric(Y2))
}