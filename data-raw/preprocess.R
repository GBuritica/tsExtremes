#### Load data
precip2   <- "data-raw/tototiti.txt"
pre       <- read.table(file = precip2, header = TRUE, sep = " ")
########## Season analysis
dd<-as.integer(pre$y.date/10000)*10000
dm<-pre$y.date-dd
dSpring<-(300<dm)&(dm<600)
dSummer<-(600<dm)&(dm<900)
dFall<-(900<dm)&(dm<1200)
dWinter<-((100<dm)&(dm<300))|((1200<dm)&(dm<1300))


pre$SEASON <- pre[,2]
pre$SEASON[dSpring] <- 'SPRING'
pre$SEASON[dSummer] <- 'SUMMER'
pre$SEASON[dFall]   <- 'FALL'
pre$SEASON[dWinter] <- 'WINTER'

rainfall <- pre
## save the data in the package
usethis::use_data(rainfall, overwrite = T)
