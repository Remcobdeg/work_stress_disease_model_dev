package('Hmisc')

x=3
das <- data.frame(anim=1:15,
                  wt=c(181,179,180.5,201,201.5,245,246.4,
                       189.3,301,354,369,205,199,394,231.3))
a <- split(das, cut(das$wt, x))

#arrange people by their number of stressors
#split the data in x equal sized groups. x being min(ceil(settings$people/50), 20)

#for each group, divide the number of sick people by the total number of people in the group

b <- as.data.frame(t(sapply(a, function(pa) sum(pa$anim)/sum(pa$wt))))
names(b) <- paste('level',seq(x)) 
b

b <- as.data.frame(sapply(a, function(pa) sum(pa$anim)/sum(pa$wt)))
names(b)[1] <- 'odds'
b$level <- paste('level',seq(x)) 
names(b) <- paste('level',seq(x)) 
b

ggplot(b, aes(x=level,y=odds)) + geom_col()


ldply(unique(daydata$days, function(day)) split(day, cut))  

ldply(split(daydata, daydata$days), function(day) )

c <- split(daydata, daydata$days)[[10]]
split(c, cut(c$stressfM,5))
c <- c[order(c$stressfM),]
c$num <- 1:10
c
