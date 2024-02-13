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
dat <- read.csv("Experiment_2.csv",header=T)

#create a new column of "Plate ID"
dat$Plate.ID <- sub("_.*","",dat$Unique.ID) 

# Check structure of data
summary(dat) 
head(dat)
tail(dat)

# Create vector of the Ages that were recorded.
Age.days.vec <- unique(dat$Age.days)

# a) Survival
#filter the data to include only the data from the final age, which contains the overall survival data
d <- dat %>% filter(Age.days == Age.days.vec[2]) 

# Generalized linear mixed effects model with a binomial distribution and random effect of plate number (Plate.ID). 
survival1 <- glmmTMB(Survival ~ Relatedness + (1|Plate.ID), family= "binomial", data=d)
survival2 <- glmmTMB(Survival ~ 1 + (1|Plate.ID), family= "binomial", data=d)

anova(survival1, survival2, test="Chisq")
# There is no evidence that neighbor relatedness affected the probability of survival (X2= 0.036, df= 1, p = 0.8495). 

# b) Relative Growth Rate
# Change data to the wide format
dat_wide <- dat %>% pivot_wider(id_cols=Unique.ID:Position,
                                names_from=Age.days, 
                                values_from=c(Survival,Bifurcations,Zooids))

## Calculate Relative Growth Rate 
dat_wide$relgrowth <- (log(dat_wide$Zooids_46)- log(dat_wide$Zooids_26)) / (46-26)

dat_wide$Plate.ID <- sub("_.*","",dat_wide$Unique.ID) #create a new column of "Plate ID"
# ch2exp2_new <- ch2exp2_wide[complete.cases(ch2exp2_wide$relgrowth), ] #remove any rows with missing relative growth values

## Generalized linear mixed effects model with a Gaussian distribution and random effect of plate number (Plate.ID). 
rgr1exp2 <- glmmTMB(relgrowth ~ Relatedness + (1|Plate.ID), data=dat_wide, family= gaussian())
rgr2exp2 <- glmmTMB(relgrowth ~ 1 + (1|Plate.ID), data=dat_wide, family= gaussian())

anova(rgr1exp2, rgr2exp2, test="Chisq")
## There is no evidence that neighbor relatedness affected... relative growth rate (X2= 0.008, df= 1, p= 0.929).



# Reproductive Output (Total Ovicells) (c)
## Generalized linear mixed effects hurdle model with truncated negative binomial distribution
d <- dat %>% filter(Age.days == Age.days.vec[2]) 
# Add zeros ovicells for colonies that died
d$Ovicells.total <- ifelse(is.na(d$Ovicells.total),0,d$Ovicells.total)
# aggregate(Ovicells.total~Plate.ID,mean,data=d)
totov1 <- glmmTMB(Ovicells.total ~ Relatedness + (1|Plate.ID), 
                  ziformula = ~.,
                  data = d, 
                  family = truncated_nbinom1)
totov2 <- glmmTMB(Ovicells.total ~ 1 + (1|Plate.ID), 
                  ziformula = ~.,
                  data = d, 
                  family = truncated_nbinom1)

anova(totov1, totov2, test= "Chisq")


# Position of colonies
s1 <- glmmTMB(Survival ~ Relatedness*Position + (1|Plate.ID), family= "binomial", data=d)
s2 <- glmmTMB(Survival ~ Relatedness+Position + (1|Plate.ID), family= "binomial", data=d)
s3 <- glmmTMB(Survival ~ Relatedness + (1|Plate.ID), family= "binomial", data=d)
anova(s1, s2, test="Chisq")
anova(s2, s3, test="Chisq")

r1 <- glmmTMB(relgrowth ~ Relatedness*Position + (1|Plate.ID), data=dat_wide, family= gaussian())
r2 <- glmmTMB(relgrowth ~ Relatedness+Position + (1|Plate.ID), data=dat_wide, family= gaussian())
r3 <- glmmTMB(relgrowth ~ Relatedness + (1|Plate.ID), data=dat_wide, family= gaussian())
anova(r1, r2, test="Chisq")
anova(r2, r3, test="Chisq")

t1 <- glmmTMB(Ovicells.total ~ Relatedness*Position + (1|Plate.ID), 
                  ziformula = ~.,
                  data = d, 
                  family = truncated_nbinom1)
t2 <- glmmTMB(Ovicells.total ~ Relatedness+Position + (1|Plate.ID), 
              ziformula = ~.,
              data = d, 
              family = truncated_nbinom1)
t3 <- glmmTMB(Ovicells.total ~ Relatedness + (1|Plate.ID), 
              ziformula = ~.,
              data = d, 
              family = truncated_nbinom1)
anova(t1, t2, test= "Chisq")
anova(t2, t3, test= "Chisq")


# Make Figure 3
related.vec <- c("Non", "Sib") #define levels of Relatedness variable

quartz(width=10,height=3.5) 
par(mfrow=c(1,3),mar=c(4,5,2,1), oma=c(0,0,0,1),cex.axis = 1.5)

# Survival (a)
newdat <- expand.grid(Relatedness = related.vec, Plate.ID=NA) #create data frame with all combinations of Relatedness levels
p <- predict(survival1,
             newdata = newdat,
             type="link",
             se.fit=T) #Generate predictions and standard errors using the model
newdat <- cbind.data.frame(newdat,
                           fit=p$fit,
                           lwr=p$fit-2*p$se.fit,
                           upr=p$fit+2*p$se.fit)

plot(c(0.5,2.5),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
with(newdat, points(1:2,plogis(fit),pch=19,cex=2))
with(newdat,segments(1:2,lwr,1:2,upr))
axis(side=1,at=1:2,labels=c("Unrelated","Siblings"))
axis(side=2,at=seq(0,1,0.2),las=1)
mtext("Probability of Survival",side=2,line=3,cex=1.2)
mtext("Relatedness of Neighbors",side=1,line=2.5,cex=1.2)
mtext("a)", side=3,adj=-0.1,cex=1.2)     
     
# Relative Growth Rate (b)
newdat2 <- expand.grid(Relatedness = related.vec, Plate.ID=NA) #create data frame with all combinations of Relatedness levels
p <- predict(rgr1exp2,
             newdata = newdat2,
             type="link",
             se.fit=T) #Generate predictions and standard errors using the model
newdat2 <- cbind.data.frame(newdat2,
                           fit=p$fit,
                           lwr=p$fit-2*p$se.fit,
                           upr=p$fit+2*p$se.fit)


# B <- as.matrix(coef(rgr1exp2)$cond$Plate.ID)
# x1 <- ifelse(grepl("non",rownames(B)),0,1)
# X <- matrix(c(rep(1,nrow(B)),x1),nrow=nrow(B),ncol=2)
# plate.effects2 <- data.frame(Eff = rowSums(X*B),
#                              Rel = substr(rownames(coef(rgr1exp2)$cond$Plate.ID),1,3)) 

nd <- data.frame(Relatedness = c(rep("Non",5),rep("Sib",10)), 
                 Plate.ID=rownames(coef(rgr1exp2)$cond$Plate.ID)) #Create a data frame with combinations of Relatedness and Plate ID

plate.effects2 <- data.frame(nd,fit=predict(rgr1exp2,
                                            newdata = nd,
                                            type="response",
                                            se=F)) #Create a complete data frame with the combinations of relatedness and plate ID, along with the predicted values using the model


plot(c(0.5,2.5),c(0,0.3),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
# Mean of each Plate estimated from the GLMM
set.seed(99); with(plate.effects2[plate.effects2$Relatedness=="Non",], 
                   points(jitter(rep(1,length(Relatedness)),4),fit,pch=19,
                          col=adjustcolor("grey",alpha.f = 0.6)))
set.seed(99); with(plate.effects2[plate.effects2$Relatedness=="Sib",], 
                   points(jitter(rep(2,length(Relatedness)),2),fit,pch=19,
                          col=adjustcolor("grey",alpha.f = 0.6)))

with(newdat2, points(1:2,fit,pch=19,cex=2))
with(newdat2,segments(1:2,lwr,1:2,upr))
axis(side=1,at=1:2,labels=c("Unrelated","Siblings"))
axis(side=2,at=seq(0,1,0.05),las=1)
mtext(side = 2, line = 5.5, cex = 1.2,adj=0.5,
      expression(paste("Relative Growth Rate")))
mtext(side = 2, line = 3.5, cex = 1.2,adj=0.5,
      expression(paste("(zooids " ,zooid^-1, day^-1,")")))
mtext("Relatedness of Neighbors",side=1,line=2.5,cex=1.2)
mtext("b)", side=3,adj=-0.1,cex=1.2)     


# Reproductive Output (c)
# Model predictions
newdat3 = data.frame(Relatedness = related.vec,Plate.ID=NA)
p <- predict(totov1,
             newdata = newdat3,
             type="response",
             se=T) 
newdat3 <- data.frame(newdat3,
                      fit=p$fit,
                      lwr=p$fit-p$se.fit,
                      upr=p$fit+p$se.fit) 

nd <- data.frame(Relatedness = c(rep("Non",5),rep("Sib",10)), 
                 Plate.ID=rownames(coef(totov1)$cond$Plate.ID)) 

plate.effects3 <- data.frame(nd,fit=predict(totov1,
                                            newdata = nd,
                                            type="response",
                                            se=F)) 

plot(c(0.5,2.5),c(0,900),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")

# Mean of each Plate estimated from the GLMM
set.seed(99); with(plate.effects3[plate.effects3$Relatedness=="Non",], 
                   points(jitter(rep(1,length(Relatedness)),4),fit,pch=19,
                          col=adjustcolor("grey",alpha.f = 0.6)))
set.seed(99); with(plate.effects3[plate.effects3$Relatedness=="Sib",], 
                   points(jitter(rep(2,length(Relatedness)),2),fit,pch=19,
                          col=adjustcolor("grey",alpha.f = 0.6)))

with(newdat3, points(1:2,fit,pch=19,cex=2))
with(newdat3,segments(1:2,lwr,1:2,upr))
axis(side=1,at=1:2,labels=c("Unrelated","Siblings"))
axis(side=2,at=seq(0,1000,200),las=1)
mtext("Reproductive Output\n(ovicells per colony)",side=2,line=4,cex=1.2)
mtext("Relatedness of Neighbors",side=1,line=2.5,cex=1.2)
mtext("c)", side=3,adj=-0.1,cex=1.2)
