library(lattice)

densityplot(10*rchisq(10000, df=4,ncp = 0))

?rlnorm()

densityplot(rlnorm(1000, meanlog = 1, sdlog = 1))
?rpois()

densityplot(rpois(10000, lambda = 2.1))
?rgamma()

densityplot(rgamma(1000, shape = 2, rate = .1))
mean(rgamma(1000, shape = 2, rate = .1))

rpois(1, 2.1)
