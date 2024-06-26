################################################################################
#                                  Fabian SL Yii                               #
#                               fabian.yii@ed.ac.uk                            #
################################################################################
library(car)
library(lme4)
library(ggplot2)
library(sjPlot)

rm(list=ls())
setwd("/Users/fabianyii/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Projects/vesselCalibre/")


################################################################################
################################ Preprocessing #################################
################################################################################
# Read dataframe containing vessel metrics derived using VAMPIRE; dataset comprised patients seen at GCU optometry clinic and vision science students
vesselData  <- read.csv("data/GCUvampireVesselResults.csv") 
# Read dataframes containing clinical information for each eye
d            <- read.csv("data/GCUdataLongCombined.csv")
names(d)[18] <- "name"
# Merge datasets
d            <- merge(vesselData, d, c("name", "eye") )
d$ID         <- factor(d$ID)
# Include data from students (not from patient records) 
# 96 right eyes of 96 participants
d <- subset(d, !fromPxRecords & cyclo==TRUE) 
# Remove eyes with pathologic myopia (3 eyes had diffuse chorioretinal atrophy) 
# 93 right eyes of 93 participants left
GCU         <- d[d$remark == "",]
# Remove eyes with keratoconus (n=1), amblyopia (n=4) and strabismus (n=6)
# 82 right eyes of 82 students left
GCU <- subset(GCU, OH=="good") 

# Compare characteristics of included and excluded participants
include   <- subset(d, remark=="" & OH=="good")
exclude   <- subset(d, remark!="" | OH!="good")
d$include <- TRUE
d$include[which(d$remark!= "" | d$OH!="good")] <- FALSE
chisq.test(table(d$sex, d$include), simulate.p.value=TRUE)             # sex distribution
t.test(subset(d, include==TRUE)$age, subset(d, include==FALSE)$age)       # mean age
t.test(subset(d, include==TRUE)$SER, subset(d, include==FALSE)$SER)       # mean SER
t.test(subset(d, include==TRUE)$meanAL, subset(d, include==FALSE)$meanAL) # mean AL


################################################################################
############################## Helper functions ################################
################################################################################
# "BRE2": AL only (https://link.springer.com/content/pdf/10.1007/BF00175988.pdf)
# "BK2" : SER only (https://link.springer.com/article/10.1007/BF00166758); used by all previous studies
# "L1"  : SER and CR (https://www.thieme-connect.com/products/ejournals/abstract/10.1055/s-2008-1055068 cited by https://link.springer.com/content/pdf/10.1007/BF00175988.pdf)

# Bennett, Rudnicka and Edgar abbreviated AL method
BRE2 <- function(calibre, AL, SER, telecentric){
  if(telecentric==TRUE){p <- 1.37}
  if(telecentric==FALSE){p <- 0.015*SER + 1.521}
  trueCalibre <- p * 0.01306*(AL - 1.82) * calibre
  return(trueCalibre) }
# Bengtsson and Krakau formula (used by all previous studies)
BK2 <- function(calibre, SER, telecentric){ 
  if(telecentric==TRUE){p <- 1.37}
  if(telecentric==FALSE){p <- 0.015*SER + 1.521}
  trueCalibre <- calibre * ((1-0.017*SER)/(p*60))
  return(trueCalibre) }
# Littmann formula
L1 <- function(calibre, SER, CR, telecentric){
  a <- 0.01 + 0.00236 * (CR - 8) 
  b <- 0.6126 + 0.0968 * (CR - 8)
  c <- 30.52 + 2.57 * (CR - 8)
  q <- (a * SER^2 - b * SER + c ) / 100
  if(telecentric==TRUE){p <- 1.37}
  if(telecentric==FALSE){p <- 0.015*SER + 1.521}
  trueCalibre <- p * q * calibre
  return(trueCalibre) }


################################################################################
################################# Regression ###################################
################################################################################
GCU$ethnicGroup <- factor(ifelse(GCU$ethnic=="White", "White", "non-White"))

## Not corrected for magnification
tab_model(lm(CRAE_Knudtson ~ meanAL + ethnicGroup + ageFundus + factor(sex), GCU))
tab_model(lm(CRVE_Knudtson ~ meanAL + ethnicGroup + ageFundus + factor(sex), GCU))

## Magnification correction assuming non-telecentricity ##
# Bennett, Rudnicka and Edgar
tab_model(lm(BRE2(CRAE_Knudtson, meanAL, SER, FALSE) ~ meanAL + ethnicGroup + ageFundus + factor(sex), GCU))
tab_model(lm(BRE2(CRVE_Knudtson, meanAL, SER, FALSE) ~ meanAL + ethnicGroup + ageFundus + factor(sex), GCU))
# Bengtsson and Krakau
tab_model(lm(BK2(CRAE_Knudtson, SER, FALSE) ~ meanAL + ethnicGroup + ageFundus + factor(sex), GCU))
tab_model(lm(BK2(CRVE_Knudtson, SER, FALSE) ~ meanAL + ethnicGroup  + ageFundus + factor(sex), GCU))
# Littmann
tab_model(lm(L1(CRAE_Knudtson, SER, meanCR, FALSE) ~ meanAL + ethnicGroup  + ageFundus + factor(sex), GCU))
tab_model(lm(L1(CRVE_Knudtson, SER, meanCR, FALSE) ~ meanAL + ethnicGroup  + ageFundus + factor(sex), GCU))


################################################################################
################################### Plots ######################################
################################################################################
## Figure 2
plotd <- read.csv("data/extractedDataFromGarwayHeathFig1.csv")
plotd <- subset(plotd, Formula != "Garway-Heath")
ggplot(data=plotd, aes(x=x, y=y, group=Formula, col=Formula)) +
  geom_point() + geom_line(size=1.2) +
  xlab("Axial length (mm)") + ylab("Error (%)") +
  geom_hline(yintercept=0, linetype="dashed", size=0.2) +
  theme_minimal() + theme(legend.position="top",
                          legend.title=element_blank(),
                          panel.grid.minor=element_blank(),
                          panel.grid.major=element_blank(),
                          axis.ticks.x=element_line(size=0.7),
                          axis.ticks.y=element_line(size=0.7)) +
  geom_segment(x=27.3, y=-9, xend=27.3, yend=-13, col="gray", size=0.3,
                 arrow = arrow(length = unit(0.5, "cm"))) +
  geom_text(x=25.7, y=-10.8, size=3, col="gray20", label="Underestimation of 'q'") +
  geom_segment(x=22.5, y=7, xend=22.5, yend=11, col="gray", size=0.3,
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_text(x=24, y=9, size=3, col="gray20", label="Overestimation of 'q'")
ggsave("fig/fig2.png", width=5.5, height=4.5)


# cm <- function(SER, a, b){return(a*SER + b) }
# plot(GCU$SER, cm(GCU$SER, 0.01, 1.14), type="l", ylim=c(0.6,2.4))
# lines(GCU$SER, cm(GCU$SER, 0.02, 1.85) )
# lines(GCU$SER, cm(GCU$SER, 0.02, 2.17) )
# lines(GCU$SER, cm(GCU$SER, 0.01, 0.94) )
# lines(GCU$SER, cm(GCU$SER, 0.03, 1.91) )
# lines(GCU$SER, cm(GCU$SER, 0.02, 1.85) )
# lines(GCU$SER, cm(GCU$SER, 0.01, 1.27) )
# lines(GCU$SER, cm(GCU$SER, 0.01, 0.76) )
# lines(GCU$SER, cm(GCU$SER, 0.01, 1.3) )
# lines(GCU$SER, cm(GCU$SER, 0.01, 2.02) )
# a <- (0.01+0.02+0.02+0.01+0.03+0.02+0.01+0.01+0.01+0.01)/10
# b <- (1.14 + 1.85 + 2.17 + 0.94 + 1.91 + 1.85 +1.27 + 0.76 + 1.3 + 2.02)/10
# lines(GCU$SER, cm(GCU$SER, a, b), lwd=3, col="red" )










