#### Create piecewise regression to determine membership categories for dominance info from percolation & conductance analyses.
library(stringr)
library(ggplot2)
library(data.table)
library(segmented)

getwd()
setwd('../Dropbox/SNH/FeedingCompetition/')
IndDomCertNC8 <-read.csv('IndivDomCertainty_NC8.csv')
DyadDomCertNC8 <- read.csv('DyadicDomCertainty_NC8.csv')
DomProbNC8 <- read.csv('DomProbs_NC8.csv')
head(IndDomCertNC8)
head(DomProbNC8)



# Create a variable for only study period (if it's not already in the dataset)
IndDomCertNC8$StudyPeriod <- str_sub(IndDomCertNC8$Initiator.SP,-1,-1)

# Add rank from DomProb to IndDomCert
### First, change the individual dom certainty and dom prob data.frames to 
### data.tables
dtIndivNC8 = data.table(IndDomCertNC8)
dtDomNC8 = data.table(DomProbNC8)

### Then set the ID.studyperiod variable to the key for both data.tables
setkey(dtIndivNC8, Initiator.SP)
setkey(dtDomNC8, ID.SP)


### Now do an inner join on dtIndiv where the key matches the key from dtDom
dtNC8 <- dtIndivNC8[dtDomNC8, nomatch = 0]

### So let's make a data table with only the information we need to do the
### piecewise regression.
dtRegNC8 <- data.frame(cert = dtNC8$Mean, 
                       sp = dtNC8$StudyPeriod
                       )
dyadNC8 <- data.frame(cert = DyadDomCertNC8$RankingCertainty, 
                      SP = DyadDomCertNC8$SP
                      )


##### We need to split this by study period before we can work with it
### For NC8 SP1:
preNC8SP1<- dtRegNC8[dtRegNC8$sp==1,]
NC8SP1 <- preNC8SP1[ do.call(order, preNC8SP1), ]
head(NC8SP1)
NC8SP1$certRank <- c(1:length(NC8SP1$sp))

### For NC8 SP2:
preNC8SP2<- dtRegNC8[dtRegNC8$sp==2,]
NC8SP2 <- preNC8SP2[ do.call(order, preNC8SP2), ]
head(NC8SP2)
NC8SP2$certRank <- c(1:length(NC8SP2$sp))

### And for NC8 SP3:
preNC8SP3<- dtRegNC8[dtRegNC8$sp==3,]
NC8SP3 <- preNC8SP3[ do.call(order, preNC8SP3), ]
head(NC8SP3)
NC8SP3$certRank <- c(1:length(NC8SP3$sp))

### For NC8 Dyads SP1:
preDyadsNC8SP1<- dyadNC8[dyadNC8$SP==1,]
dyadsNC8SP1 <- preDyadsNC8SP1[ do.call(order, preDyadsNC8SP1), ]
head(dyadsNC8SP1)
dyadsNC8SP1$certRank <- c(1:length(dyadsNC8SP1$SP))

### For NC8 Dyads SP2:
preDyadsNC8SP2<- dyadNC8[dyadNC8$SP==2,]
dyadsNC8SP2 <- preDyadsNC8SP2[ do.call(order, preDyadsNC8SP2), ]
head(dyadsNC8SP2)
dyadsNC8SP2$certRank <- c(1:length(dyadsNC8SP2$SP))

### And for NC8 Dyads SP3:
preDyadsNC8SP3<- dyadNC8[dyadNC8$SP==3,]
dyadsNC8SP3 <- preDyadsNC8SP3[ do.call(order, preDyadsNC8SP3), ]
head(dyadsNC8SP3)
dyadsNC8SP3$certRank <- c(1:length(dyadsNC8SP3$SP))

# In order to use the segmented package, we need to have estimates about the
# break points for the breakpoints in our piecewise regression. We can do
# plot the data using ggplot. We want our x-value to be rank order of 
# certainty and our y-value to be dominance mean certainty (sorted ascending, which we did
# in the last step)
ggplot(NC8SP1, aes(x = certRank, y = cert)) +
  geom_point()

ggplot(NC8SP2, aes(x = certRank, y = cert)) +
  geom_point()

ggplot(NC8SP3, aes(x = certRank, y = cert)) +
  geom_point()

### For NC8 dyads:
ggplot(dyadsNC8SP1, aes(x = certRank, y = cert)) +
  geom_point()

ggplot(dyadsNC8SP2, aes(x = certRank, y = cert)) +
  geom_point()

ggplot(dyadsNC8SP3, aes(x = certRank, y = cert)) +
  geom_point()


# Now we'll use the segmented package. First create a linear model using rank
# as a function of certainty from the appropriate dataset, then feed that information
# into segmented.default with seg.Z as the x variable (rank) and psi as the estimated
# breakpoints in a list with the same name as seg.Z
glmNC8SP1 <-glm(cert ~ certRank, data = NC8SP1)
glmNC8SP2 <-glm(cert ~ certRank, data = NC8SP2)
glmNC8SP3 <-glm(cert ~ certRank, data = NC8SP3)

glmdyadsNC8SP1 <-glm(cert ~ certRank, data = dyadsNC8SP1)
glmdyadsNC8SP2 <-glm(cert ~ certRank, data = dyadsNC8SP2)
glmdyadsNC8SP3 <-glm(cert ~ certRank, data = dyadsNC8SP3)


#Run the breakpoint analyses for each SP:
segNC8SP1 <- segmented(glmNC8SP1,
                       seg.Z = ~certRank,
                       psi=list(certRank=c(15, 75)),
                       control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100000)
                        )
segNC8SP2 <- segmented(glmNC8SP2,
                       seg.Z = ~certRank,
                       psi=list(certRank=c(20, 60)),
                       control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100000)
                        )

segNC8SP3 <- segmented(glmNC8SP3,
                       seg.Z = ~certRank,
                       psi=list(certRank=c(20, 60)),
                       control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100000)
                        )

### Now run the segmented analyses for dyads:
segNCdyads8SP1 <- segmented(glmdyadsNC8SP1,
                            seg.Z = ~certRank,
                            psi=list(certRank=c(1500, 3400)),
                            control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100000)
                            )
segdyadsNC8SP2 <- segmented(glmdyadsNC8SP2,
                            seg.Z = ~certRank,
                            psi=list(certRank=c(900, 2600)),
                            control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100000)
                            )

segdyadsNC8SP3 <- segmented(glmdyadsNC8SP3,
                            seg.Z = ~certRank,
                            psi=list(certRank=c(800, 2200)),
                            control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100000)
                            )
# And retrieve the calculated breakpoints using the following:
segNC8SP1$psi
segNC8SP2$psi
segNC8SP3$psi

save(segNC8SP1, file = "NC8sp1_indiv.rda")
save(segNC8SP2, file = "NC8sp2_indiv.rda")
save(segNC8SP3, file = "NC8sp3_indiv.rda")

#For dyadic info too:
segNCdyads8SP1$psi
segdyadsNC8SP2$psi
segdyadsNC8SP3$psi

save(segNCdyads8SP1, file = "NC8sp1_dyad.rda")
save(segdyadsNC8SP2, file = "NC8sp2_dyad.rda")
save(segdyadsNC8SP3, file = "NC8sp3_dyad.rda")
