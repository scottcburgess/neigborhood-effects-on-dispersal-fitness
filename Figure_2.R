### This code produces Figure 2 in the manuscript: 
### Barnes DK, Burgess SC. Fitness consequences of marine larval dispersal: 
# the role of neighborhood density, spatial arrangement, 
# and genetic relatedness on survival, growth, reproduction, and paternity
# Code finalized Feb 2024
# Code written by Danie Barnes with input from Scott Burgess
# Any comments or error reporting, please contact Scott Burgess: sburgess@bio.fsu.edu

# sessionInfo()
# R version 4.3.2

# Load required libraries
library('glmmTMB') # version 1.1.7
library('dplyr')
library('tidyr')

# Import data
dat <- read.csv("Experiment_1.csv",header=T)
dat$Age.days <- factor(dat$Age.days)
# Check data structure
summary(dat)
head(dat)
tail(dat)

# create new column for colony identity
dat$colonyidentity <- interaction(dat$Unique.ID, dat$Grid.position) 

# Create vector of the Ages that were recorded.
Age.days.vec <- unique(dat$Age.days)

# Run analyses on the probability of SURVIVAL after 8 days in the field
## Generalized linear mixed effects model with a binomial distribution and random effect of plate number (Unique.ID). 
d <- dat %>% filter(Age.days == Age.days.vec[2])
survival1 <- glmmTMB(Survival ~ Density * Relatedness + (1|Unique.ID), 
                     family= "binomial", data=d) #interaction between density and relatedness
survival2 <- glmmTMB(Survival ~ Density + Relatedness + (1|Unique.ID), 
                     family="binomial", data=d) #additive model
survival3 <- glmmTMB(Survival ~ Density + (1|Unique.ID), 
                     family="binomial", data=d)
survival4 <- glmmTMB(Survival ~ Relatedness + (1|Unique.ID), 
                     family="binomial", data=d)

anova(survival1, survival2,test="Chisq") # interaction at 8 days
anova(survival2, survival3,test="Chisq") # relatedness at 8 days 
anova(survival2, survival4,test="Chisq") # density at 8 days 

# Run analyses on the probability of SURVIVAL after 15 days in the field
## Generalized linear mixed effects model with a binomial distribution and random effect of plate number
# At the second time point of collection, there is a negative effect of density on survival (p = 0.003866) and there is no effect of relatedness on survival (P= 0.1242). 
d <- dat %>% filter(Age.days==Age.days.vec[3])
survival5 <- glmmTMB(Survival ~ Density * Relatedness + (1|Unique.ID), 
                     family="binomial", data=d)
survival6 <- glmmTMB(Survival ~ Density + Relatedness + (1|Unique.ID), 
                     family="binomial", data=d)
survival7 <- glmmTMB(Survival ~ Density + (1|Unique.ID), 
                     family="binomial", data=d)
survival8 <- glmmTMB(Survival ~ Relatedness + (1|Unique.ID), 
                     family="binomial", data=d)

anova(survival5, survival6,test="Chisq") # interaction at 15 days
anova(survival6, survival7,test="Chisq") # relatedness at 15 days 
anova(survival6, survival8,test="Chisq") # density at 15 days 


# Calculate the odds of colony survival (using the additive model)
# 8 days
A <- summary(survival2)$coefficients$cond[2,1] 
round(1-exp(A),4) * 100
round(1-exp(A + c(-2,2) * summary(survival2)$coefficients$cond[2,2]),4) *100

# 15 days
B <- summary(survival6)$coefficients$cond[2,1]
round(1-exp(B),4) * 100
round(1-exp(B + c(-2,2) * summary(survival2)$coefficients$cond[2,2]),4) *100



# RELATIVE GROWTH RATE
# Change data to the wide format
dat_wide <- dat %>% pivot_wider(id_cols=Unique.ID:Inside.outside,
                                names_from=Age.days, 
                                values_from=c(Survival,Bifurcations,Zooids))

# Calculate relative growth rate 
## At collection 1
dat_wide$relgrowth.collection1 <- (log(dat_wide$Zooids_25)- log(dat_wide$Zooids_17)) / (25-17)
## At collection 2
dat_wide$relgrowth.collection2 <- (log(dat_wide$Zooids_32)- log(dat_wide$Zooids_17)) / (32-17)

# Run analyses on the relative growth rate (number of zooids produced per zooid per day) after 8 days in the field
## Generalized linear mixed effects model with a Gaussian distribution and random effect of plate number
rgr1 <- glmmTMB(relgrowth.collection1 ~ Density * Relatedness + (1|Unique.ID), data=dat_wide, family= gaussian())
rgr2 <- glmmTMB(relgrowth.collection1 ~ Density + Relatedness + (1|Unique.ID), data=dat_wide, family= gaussian())
rgr3 <- glmmTMB(relgrowth.collection1 ~ Density + (1|Unique.ID), data=dat_wide, family= gaussian())
rgr4 <- glmmTMB(relgrowth.collection1 ~ Relatedness + (1|Unique.ID), data= dat_wide, family= gaussian())
rgr_overall_1 <- glmmTMB(relgrowth.collection1 ~ 1 + (1|Unique.ID), data= dat_wide, family= gaussian())

anova(rgr1,rgr2) # interaction, 8 days
anova(rgr2,rgr3) # relatedness, 8 days
anova(rgr2,rgr4) # density, 8 days

# Run analyses on the relative growth rate (number of zooids produced per zooid per day) after 15 days in the field
## Generalized linear mixed effects model with a Gaussian distribution and random effect of plate number. 
rgr5 <- glmmTMB(relgrowth.collection2 ~ Density * Relatedness + (1|Unique.ID), data=dat_wide, family= gaussian())
rgr6 <- glmmTMB(relgrowth.collection2 ~ Density + Relatedness + (1|Unique.ID), data=dat_wide, family= gaussian())
rgr7 <- glmmTMB(relgrowth.collection2 ~ Density + (1|Unique.ID), data=dat_wide, family= gaussian())
rgr8 <- glmmTMB(relgrowth.collection2 ~ Relatedness + (1|Unique.ID), data= dat_wide, family= gaussian())
rgr_overall_2 <- glmmTMB(relgrowth.collection2 ~ 1 + (1|Unique.ID), data= dat_wide, family= gaussian())

anova(rgr5,rgr6) # interaction, 15 days
anova(rgr6,rgr7) # relatedness, 15 days
anova(rgr6,rgr8) # density, 15 days


# Calculate zooids produced per zooid per day, overall
B <- summary(rgr_overall_1)$coefficients$cond[1,1]; round(B,2)
round(B + c(-2,2) * summary(rgr_overall_1)$coefficients$cond[1,2],2)

B <- summary(rgr_overall_2)$coefficients$cond[1,1]; round(B,2)
round(B + c(-2,2) * summary(rgr_overall_2)$coefficients$cond[1,2],2)

# Make Figure 2
quartz(width=6,height=5) 
par(mfrow=c(2,2),mar=c(4,5,2,0), oma=c(0,0,0,1))

# a)
Density.vec <- seq(0,20,0.1) #creating a density vector, to make continuous 
pred.frame <- with(dat, expand.grid(Density=Density.vec,Relatedness=unique(Relatedness))) #matrix by everything
X <- model.matrix(delete.response(terms(survival2)),pred.frame) #turning it into a model matrix
pred.frame$response <- drop(X %*% fixef(survival2)[['cond']]) 
V <- diag(X %*% vcov(survival2)[['cond']] %*% t(X)) #variance to use for the confidence intervals 
pred.frame$upr <- pred.frame$response + 1.96*sqrt(V) #create confidence intervals 
pred.frame$lwr <- pred.frame$response - 1.96*sqrt(V)

plot(c(0,21),c(0,1),type="n",ylab="",xlab="",yaxt="n",xaxt="n",bty="l") #add axes lines
with(pred.frame[pred.frame$Relatedness=="Non",], polygon(c(Density,rev(Density)),
                                                         c(plogis(lwr),rev(plogis(upr))),col=adjustcolor("dodgerblue",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Non",], lines(Density, plogis(response),
                                                       col="dodgerblue",lwd=2))
with(pred.frame[pred.frame$Relatedness=="Sib",], polygon(c(Density,rev(Density)),
                                                         c(plogis(lwr),rev(plogis(upr))),col=adjustcolor("tomato",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Sib",], lines(Density, plogis(response),
                                                       col="tomato",lwd=2))
axis(side = 1, at= seq(0,20,5), cex = 0.25) #add numbers to the x axis
axis(side = 2, at= seq(0,1,0.1), cex = 0.25, las=1) #add numbers to the y axis
mtext(side= 1, line = 3, cex= 1, expression(paste("Density (colonies per 12", cm^2,")"))) #add title to the x axis
mtext(side =2, line = 2.5, cex = 1, expression("Probability of survival")) #add title to the y axis
legend(1,0.4,legend= c("Siblings", "Unrelated"), lwd=2, col=c("tomato","dodgerblue"), bty="n")
mtext("a) 8 days in the field", side=3, line=0.5, adj=0, las=1, cex=1)


# b)
Density.vec <- seq(0,20,0.1) #creating a density vector, to make continuous 
pred.frame <- with(dat, expand.grid(Density=Density.vec,Relatedness=unique(Relatedness))) #matrix by everything you have
X <- model.matrix(delete.response(terms(survival6)),pred.frame) #turning it into a model matrix
pred.frame$response <- drop(X %*% fixef(survival6)[['cond']]) #
V <- diag(X %*% vcov(survival6)[['cond']] %*% t(X)) #how you get the variance to use for the confidence intervals 
pred.frame$upr <- pred.frame$response + 1.96*sqrt(V) #create confidence intervals: plus 1.96 * standard error 
pred.frame$lwr <- pred.frame$response - 1.96*sqrt(V)
plot(c(0,21),c(0,1),type="n",ylab="",xlab="",yaxt="n",xaxt="n",bty="l")
with(pred.frame[pred.frame$Relatedness=="Non",], polygon(c(Density,rev(Density)),
                                                         c(plogis(lwr),rev(plogis(upr))),col=adjustcolor("dodgerblue",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Non",], lines(Density, plogis(response),
                                                       col="dodgerblue",lwd=2))
with(pred.frame[pred.frame$Relatedness=="Sib",], polygon(c(Density,rev(Density)),
                                                         c(plogis(lwr),rev(plogis(upr))),col=adjustcolor("tomato",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Sib",], lines(Density, plogis(response),
                                                       col="tomato",lwd=2))
axis(side = 1, at= seq(0,20,5), cex = 0.25)
axis(side = 2, at= seq(0,1,0.1), cex = 0.25, las=1)
mtext(side= 1, line = 3, cex= 1, expression(paste("Density (colonies per 12", cm^2,")")))
mtext(side =2, line = 2.5, cex = 1, expression("Probability of survival"))
legend(1,0.4,legend= c("Siblings", "Unrelated"), lwd=2, col=c("tomato","dodgerblue"), bty="n")
mtext("b) 15 days in the field", side=3, line=0.5, adj=0, las=1, cex= 1)


# c) 
Density.vec <- seq(0,20,0.1) #creating a density vector, to make continuous 

pred.frame <- with(dat_wide, expand.grid(Density=Density.vec,Relatedness=unique(Relatedness))) #matrix by everything
X <- model.matrix(delete.response(terms(rgr2)),pred.frame) #turning it into a model matrix
pred.frame$response <- drop(X %*% fixef(rgr2)[['cond']])
V <- diag(X %*% vcov(rgr2)[['cond']] %*% t(X)) #variance to use for the confidence intervals 
pred.frame$upr <- pred.frame$response + 1.96*sqrt(V) #create confidence intervals: plus 1.96 * standard error 
pred.frame$lwr <- pred.frame$response - 1.96*sqrt(V)

plot(c(0,21),c(0,0.25),type="n",ylab="",xlab="",yaxt="n",xaxt="n",bty="l") #create a blank plot
with(pred.frame[pred.frame$Relatedness=="Non",], polygon(c(Density,rev(Density)),
                                                         c((lwr),rev((upr))),col=adjustcolor("dodgerblue",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Non",], lines(Density, (response),
                                                       col="dodgerblue",lwd=2))
with(pred.frame[pred.frame$Relatedness=="Sib",], polygon(c(Density,rev(Density)),
                                                         c((lwr),rev((upr))),col=adjustcolor("tomato",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Sib",], lines(Density, (response),
                                                       col="tomato",lwd=2))
axis(side = 1, at= seq(0,20,5), cex = 0.25)
axis(side = 2, at= seq(0,0.25,0.05), cex = 0.25, las=1)
mtext(side= 1, line = 3, cex= 1, expression(paste("Density (colonies per 12", cm^2,")")))
mtext(side = 2, line = 3.9, cex = 1,adj=0.5,
      expression(paste("Relative Growth Rate")))
mtext(side = 2, line = 2.5, cex = 1,adj=0.5,
      expression(paste("(zooids " ,zooid^-1, day^-1,")")))
legend(1,0.1,legend= c("Siblings", "Unrelated"), lwd=2, col=c("tomato","dodgerblue"), bty="n")
mtext("c)", side=3, line=0.5, adj=0, las=1, cex=1.2)

# d)
pred.frame <- with(dat_wide, expand.grid(Density=Density.vec,Relatedness=unique(Relatedness)))
X <- model.matrix(delete.response(terms(rgr6)),pred.frame) #turning it into a model matrix
pred.frame$response <- drop(X %*% fixef(rgr6)[['cond']])
V <- diag(X %*% vcov(rgr6)[['cond']] %*% t(X)) #variance to use for the confidence intervals 
pred.frame$upr <- pred.frame$response + 1.96*sqrt(V) #create confidence intervals: plus 1.96 * standard error 
pred.frame$lwr <- pred.frame$response - 1.96*sqrt(V)

plot(c(0,21),c(0,0.25),type="n",ylab="",xlab="",yaxt="n",xaxt="n",bty="l")
with(pred.frame[pred.frame$Relatedness=="Non",], polygon(c(Density,rev(Density)),
                                                         c((lwr),rev((upr))),col=adjustcolor("dodgerblue",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Non",], lines(Density, (response),
                                                       col="dodgerblue",lwd=2))
with(pred.frame[pred.frame$Relatedness=="Sib",], polygon(c(Density,rev(Density)),
                                                         c((lwr),rev((upr))),col=adjustcolor("tomato",alpha.f=0.4),border=NA))
with(pred.frame[pred.frame$Relatedness=="Sib",], lines(Density, (response),
                                                       col="tomato",lwd=2))
axis(side = 1, at= seq(0,20,5), cex = 0.25)
axis(side = 2, at= seq(0,0.25,0.05), cex = 0.25, las=1)
mtext(side= 1, line = 3, cex= 1, expression(paste("Density (colonies per 12", cm^2,")")))
mtext(side = 2, line = 3.9, cex = 1,adj=0.5,
      expression(paste("Relative Growth Rate")))
mtext(side = 2, line = 2.5, cex = 1,adj=0.5,
      expression(paste("(zooids " ,zooid^-1, day^-1,")")))
legend(1,0.1,legend= c("Siblings", "Unrelated"), lwd=2, col=c("tomato","dodgerblue"), bty="n")
mtext("d)", side=3, line=0.5, adj=0, las=1, cex=1)

