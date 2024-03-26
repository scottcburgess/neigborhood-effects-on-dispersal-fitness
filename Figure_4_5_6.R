### This code produces Figure 4 in the manuscript: 
### Barnes DK, Burgess SC. Fitness consequences of marine larval dispersal: 
# the role of neighborhood density, spatial arrangement, 
# and genetic relatedness on survival, growth, reproduction, and paternity
# Code finalized Feb 2024
# Code written by Danie Barnes with input from Scott Burgess
# Any comments or error reporting, please contact Scott Burgess: sburgess@bio.fsu.edu

# sessionInfo()
# R version 4.3.2

# Load required libraries

library('glmmTMB')
library('lme4')
library('DHARMa')
library('vegan') # for rarefy
library('dplyr')
library('tidyverse')
library('tidygraph')
library('ggraph')

## Import fecundity data
dat <- read.csv("Experiment_3.csv")
# Only use the summed data over all time periods
dat <- dat %>% filter(Time_days=="all") %>% select(-Time_days, -Bifurcations, -Zooids)
# Add ID for focal colony
focals <- c("A","B","D","G")
dat$focal <- ifelse(dat$Position %in% focals,1,0)
# Set the factor levels and create a vector 
treatment.vec <- c("alone","far","near","both")
dat$Treatment <- factor(dat$Treatment, levels=treatment.vec)

## Import Paternity data
Paternity_dat_0_1 <- read.table("D0_1_Paternity.txt",header=T,sep=",")
Paternity_dat_0_5 <- read.table("D0_5_Paternity.txt",header=T,sep=",")
Paternity_dat_0_9 <- read.table("D0_9_Paternity.txt",header=T,sep=",")
# Add sample data to each data frame
brks <- seq(0,160,0.1)
add_sample_data <- function(d){
  d$MotherID <- rapply(strsplit(d$OffspringID,"_"),function(x) head(x,1))
  d$FatherID <- rapply(strsplit(d$InferredDad1,"_"),function(x) head(x,1))
  tmp <- dat[match(d$MotherID,dat$Colony),-1]
  names(tmp) <- paste0("Mother.",names(tmp))
  d <- cbind.data.frame(d,tmp)
  tmp <- dat[match(d$FatherID,dat$Colony),-1]
  names(tmp) <- paste0("Father.",names(tmp))
  d <- cbind.data.frame(d,tmp)
  
  # Calculate euclidean distance between observed parents
  d$Distance <- NA
  for(i in 1:nrow(d)){
    d$Distance[i] <- dist(
      matrix(d[i,which(names(d) %in% c("Mother.X","Father.X","Mother.Y","Father.Y"))],2,2,byrow=T),method="euclidean")
  }
  
  return(d)
}
Paternity_dat_0_1 <- add_sample_data(Paternity_dat_0_1)
Paternity_dat_0_5 <- add_sample_data(Paternity_dat_0_5)
Paternity_dat_0_9 <- add_sample_data(Paternity_dat_0_9)


## Import inferred fathers data
BestCluster_0_1 <- read.csv("D0_1_BestCluster.csv",header=T) 
BestCluster_0_5 <- read.csv("D0_5_BestCluster.csv",header=T) 
BestCluster_0_9 <- read.csv("D0_9_BestCluster.csv",header=T) 

## Import Full Sib family data 
FSFamily_0_9 <- read.csv("D0_9_BestFSFamily.csv")

## Import Offspring IDs for those that were genotyped
OffspringGenotypes <- read.table("OffspringGenotype.txt",header=T,sep="")[,1] 





# Calculate the number of offspring sampled per mother
d <- rapply(strsplit(OffspringGenotypes,"_"),function(x) head(x,1))
MotherID <- paste0(d,"_parent")
offspring_per_mother <- as.data.frame(table(MotherID))


# How many settlers were genotyped?
with(offspring_per_mother,table(Freq))
nrow(offspring_per_mother) # 25 mothers
sum(offspring_per_mother$Freq) # 619 offspring
nrow(dat) # 32 potential fathers


##### Figure 5 - Inclusion and exclusion probabilities ##########
# Set breaks for plotting
brks <- seq(0,1,0.01)

# Get frequency of Inclusion probabilities
Inc <- with(FSFamily_0_9, hist(Prob.Inc..,breaks=brks,plot=F))
# Get frequency of Exclusion probabilities
Exc <- with(FSFamily_0_9,hist(Prob.Exc..,breaks=brks,plot=F))


quartz(width=5,height=2.5)
par(mfrow=c(1,2),mar=c(2,2,1,1),oma=c(3,3,0,0))

plot(c(0,1),c(0,500),type="n",bty="l",ylab="",xlab="",xaxt="n",yaxt="n")
d <- cbind.data.frame(mids=Inc$mids, counts=Inc$counts)
y <- d %>% filter(counts>0)
axis(side=1,at=seq(0,1,0.1))
axis(side=2,at=seq(0,500,100),las=1)
with(y, segments(mids,
                 rep(0,length(mids)),
                 mids,
                 counts,
                 lend=2,lwd=3,
                 col=ifelse(mids<0.7,"grey","black")))
mtext(side=3,"a)",adj=0)
mtext(side=1,"Inclusion\nprobability",line=3)
mtext(side=2,"Frequency of\nfull sib families",line=3)

plot(c(0,1),c(0,500),type="n",bty="l",ylab="",xlab="",xaxt="n",yaxt="n")
d <- cbind.data.frame(mids=Exc$mids, counts=Exc$counts)
y <- d %>% filter(counts>0)
axis(side=1,at=seq(0,1,0.1))
axis(side=2,at=seq(0,500,100),las=1)
with(y, segments(mids,
                 rep(0,length(mids)),
                 mids,
                 counts,
                 lend=2,lwd=3,col="black"))
mtext(side=3,"b)",adj=0)
mtext(side=1,"Exclusion\nprobability",line=3)
###############################################################




# Remove offspring from full sib families with Inclusion probabilities < 0.7
tmp <- FSFamily_0_9 %>% filter(Prob.Inc.. < 0.7)
offspring_to_exclude <- c(tmp$Member1,tmp$Member2)
# just a check
# Paternity_dat_0_9 %>% filter((OffspringID %in% offspring_to_exclude))

# How many kept and excluded?
unique(BestCluster_0_9$MotherID) # 25 mothers
nrow(FSFamily_0_9) # 538 full sib families
nrow(FSFamily_0_9 %>% filter(Prob.Inc.. > 0.7)) # 525; (525/538)*100
nrow(FSFamily_0_9 %>% filter(Prob.Inc.. == 1)) # 479; (479/538)*100
nrow(tmp) # 14 full sib families excluded
length(offspring_to_exclude)  # 28 offspring excluded


# Calculate the number of unique sires per focal colony
n_fathers <- function(d){
  d$Mother_ID <- rapply(strsplit(d$MotherID,"_"),function(x) head(x,1))
  
  foo1 <- d %>% 
    group_by(Mother_ID,FatherID) %>% 
    summarise(n=n()) %>% 
    mutate(freq=n/sum(n)) %>% 
    ungroup()
  
  foo2 <- foo1 %>% group_by(Mother_ID) %>% 
    summarise(nGenotyped=sum(n),
              UniqueNumberFather=n_distinct(FatherID))
  
  foo2$standardized.n.father <- with(foo1, tapply(n, Mother_ID, rarefy, sample=20))
  foo2$standardized.n.father <- round(foo2$standardized.n.father,3)
    # add sample info 
  foo2 <- cbind.data.frame(foo2,dat[match(foo2$Mother_ID,dat$Colony),c(1:7,12)])
  
  return(foo2)
}

# n.fathers_0_1 <- n_fathers(BestCluster_0_1)
# n.fathers_0_5 <- n_fathers(BestCluster_0_5)
BestCluster_0_9 <- BestCluster_0_9 %>% filter(!(OffspringID %in% offspring_to_exclude))
n.fathers_0_9 <- n_fathers(d=BestCluster_0_9)
sort(n.fathers_0_9$standardized.n.father)
# n.fathers_0_1$UniqueNumberFather
# n.fathers_0_5$UniqueNumberFather
# n.fathers_0_9$UniqueNumberFather

# Calculate proportion for a beta glmm
n.fathers_0_9$standardized.prop.father <- n.fathers_0_9$standardized.n.father / 20
# make the 1's 0.9999 so can fit beta glmm
n.fathers_0_9$standardized.prop.father <- ifelse(n.fathers_0_9$standardized.prop.father==1,0.99999,n.fathers_0_9$standardized.prop.father)


m1 <- glmmTMB(standardized.n.father  ~ Treatment + (1|Block), 
              data = n.fathers_0_9)
m2 <- glmmTMB(standardized.n.father  ~ 1 + (1|Block), 
              data = n.fathers_0_9)
anova(m1,m2)

# plot(simulateResiduals(fittedModel=m1))


newdat_b = data.frame(Treatment = treatment.vec,Block=NA)
newdat_b$Treatment <- factor(newdat_b$Treatment, levels=treatment.vec)
p <- predict(m1,
             newdata = newdat_b,
             type="link",
             se=T)

# n.fathers_0_9 %>% group_by(Treatment) %>% summarize(mean(standardized.prop.father))
# n.fathers_0_9 %>% group_by(Treatment) %>% summarize(mean(standardized.n.father))

newdat_b <- data.frame(newdat_b,
                     fit=p$fit,
                     lwr=p$fit-2*p$se.fit,
                     upr=p$fit+2*p$se.fit) 
# block.means_b <- expand.grid(Treatment = treatment.vec, 
#                            Block=rownames(coef(m1)$cond$Block)) 
# block.means_b  <- data.frame(block.means_b,fit=predict(m1,
#                                                    newdata = block.means_b,
#                                                    type="response",
#                                                    se=F)) 



# Proportion of offspring from an assigned father
# Get the mother ID's with an assigned father
offspring.proportion <- function(a,b){
  d1 <- a %>% 
    distinct(MotherID,.keep_all=T)
  
  b$Mother_ID <- rapply(strsplit(b$MotherID,"_"),function(x) head(x,1))
  b$Father_ID <- rapply(strsplit(b$FatherID,"_"),function(x) head(x,1))
  
  d2 <- b %>% 
    filter(Mother_ID %in% d1$MotherID) %>% 
    group_by(Mother_ID,Father_ID) %>% 
    summarise(n=n()) %>% 
    mutate(freq=n/sum(n)) %>% 
    ungroup()
  
  d3 <- d2 %>% 
    group_by(Mother_ID) %>% 
    mutate(max.freq=max(freq)) %>% 
    ungroup()
  
  d4 <- d3[grep("#",d3$Father_ID,invert=T),]
  d5 <- cbind.data.frame(d4, d1[match(d4$Mother_ID,d1$MotherID),])
  d6 <- d5 %>% select(Mother_ID, Father_ID,
                      n, freq, max.freq,
                      ProbDad1,
                      Mother.Block, Mother.Treatment, Mother.Position, Mother.Direction,
                      Father.Block, Father.Treatment, Father.Position, Father.Direction)
  return(d6)
}

# offspring.proportion_0_1 <- offspring.proportion(Paternity_dat_0_1, BestCluster_0_1)
# offspring.proportion_0_5 <- offspring.proportion(Paternity_dat_0_5, BestCluster_0_5)
offspring.proportion_0_9 <- offspring.proportion(Paternity_dat_0_9, BestCluster_0_9)

# length(grep("#",BestCluster_0_1$FatherID, invert=T)) 
# length(grep("#",BestCluster_0_5$FatherID, invert=T)) 
length(grep("#",BestCluster_0_9$FatherID, invert=T))
# (12/nrow(BestCluster_0_1))*100
# (14/nrow(BestCluster_0_1))*100
(18/nrow(BestCluster_0_1))*100


# Function to make plot for Figure 4c
offspring.proportion.plot <- function(x){
  tmp <- x %>% filter(Mother.Block==Father.Block)
  # far
  y1 <- tmp %>% filter(Mother.Position=="B" & Father.Position == "C" |
                         Mother.Position=="C" & Father.Position == "B")
  # near 
  y2 <- tmp %>% filter(Mother.Position=="G" & Father.Position %in% c("H"))
  # both
  y3 <- tmp %>% filter(Mother.Position=="D" & Father.Position %in% c("E","F"))
  
  plot(c(0.7,3.3),c(0,0.3),type="n",bty="l",ylab="",xlab="",xaxt="n",yaxt="n")
  mtext(side=3,"c)",adj=0,cex.axis=1.2)
  axis(side=1,at=c(1,2,2.8,3.2),labels=c("far","near","",""),cex.axis=1.2)
  text(c(2.8,3.2),c(-0.03,-0.03),c("far","near"),cex=1.2,xpd=T,srt=45,adj=1)
  text(3,-0.08,"both",xpd=T,cex=1.2)
  axis(side=2,at=seq(0,1,0.05),las=1,cex.axis=1.2)
  mtext(side=2,"Proportion of offspring\nsired by nearest colony",line=3.5)

  # far
  with(y1, points(rep(1,length(freq)),
                                freq,pch=19,col=adjustcolor(alpha.f=0.6,"black")))
  with(y1, segments(rep(1,length(freq)),
                                  freq,
                                  rep(1,length(freq)),
                                  max.freq))
  # near
  with(y2, points(rep(2,length(freq)),
                                freq,pch=19,col=adjustcolor(alpha.f=0.6,"black")))
  with(y2, segments(rep(2,length(freq)),
                                  freq,
                                  rep(2,length(freq)),
                                  max.freq))
  # both
  with(y3[y3$Father.Position=="E"], points(rep(2.8,length(freq)),
                                           freq,pch=19,col=adjustcolor(alpha.f=0.6,"black")))
  with(y3, segments(rep(2.8,length(freq)),
                                  freq,
                                  rep(2.8,length(freq)),
                                  max.freq))
  with(y3[y3$Father.Position=="F"], points(rep(3.2,length(freq)),
                                                         freq,pch=19,col=adjustcolor(alpha.f=0.6,"black")))
  with(y3, segments(rep(3.2,length(freq)),
                                  freq,
                                  rep(3.2,length(freq)),
                                  max.freq))
  
  # Ignore error message, it just means there was no data for a certain position
}


# offspring.proportion.plot(offspring.proportion_0_1)
# offspring.proportion.plot(offspring.proportion_0_5)
# offspring.proportion.plot(offspring.proportion_0_9)



# Relative growth rate
m1 <- lm(rgr40 ~ Treatment, data=dat)
m2 <- lm(rgr40 ~ 1, data=dat)
anova(m1,m2,test="F")


# Reproductive output
m1 <- glmmTMB(Offspring ~ Treatment + (1|Block), 
              ziformula = ~.,
              data = dat, 
              family = truncated_nbinom1)
m2 <- glmmTMB(Offspring ~ 1 + (1|Block), 
              ziformula = ~.,
              data = dat, 
              family = truncated_nbinom1)
anova(m1,m2)

# simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
# plot(simulationOutput)

newdat = data.frame(Treatment = treatment.vec,Block=NA)
newdat$Treatment <- factor(newdat$Treatment, levels=treatment.vec)
p <- predict(m1,
             newdata = newdat,
             type="response",
             se=T) 
newdat <- data.frame(newdat,
                     fit=p$fit,
                     lwr=p$fit-p$se.fit,
                     upr=p$fit+p$se.fit) 
block.means <- expand.grid(Treatment = treatment.vec, 
                           Block=rownames(coef(m1)$cond$Block)) 
block.means  <- data.frame(block.means,fit=predict(m1,
                                                   newdata = block.means,
                                                   type="response",
                                                   se=F)) 


# FIGURE 4
quartz(height=3,width=7)
par(mfrow=c(1,3),mar=c(4,5,3,1),oma=c(0,3,0,0))
# Figure 4a
plot(c(0.8,4.2),c(0,2000),type="n",bty="l",ylab="",xlab="",xaxt="n",yaxt="n")
mtext(side=3,"a)",adj=0,cex.axis=1.2)
axis(side=1,at=1:4,labels=treatment.vec,cex.axis=1.2)
# mtext(side=1,"Treatment",line=3)
axis(side=2,at=seq(0,3000,200),cex.axis=1.2,las=1)
mtext(side=2,"Reproductive output\n(number of settlers)",line=4)
set.seed(99); with(block.means,points(jitter(as.numeric(Treatment),0.2),
                                      fit,
                                      pch=19,
                                      col=adjustcolor(alpha.f=0.6,col="grey")))
with(newdat, points(Treatment,fit,
                    pch=19,
                    cex=1.5,
                    col="black"))
with(newdat,segments(as.numeric(Treatment),lwr,as.numeric(Treatment),upr))

# Figure 4b
plot(c(0.8,4.7),c(5,21),type="n",ylab="",xlab="",bty="l",xaxt="n",yaxt="n")
mtext(side=3,"b)",adj=0,cex.axis=1.2)
mtext(side=1,"Treatment",line=3)
axis(side=1,at=1:4,labels=NA)
text(1:4,rep(3,4),labels=levels(n.fathers_0_9$Treatment),cex=1.2,xpd=T)
axis(side=2,at=seq(0,40,2),las=1,cex.axis=1.2)
mtext("Number of unique sires per colony\n(standardized to a sample of 20)",side=2,line=3)
y <- n.fathers_0_9 %>% filter(focal==1)
set.seed(99); with(y, points(jitter(as.numeric(Treatment),0.2),
                             standardized.n.father,
                             col=adjustcolor(alpha.f=0.4,col="grey"),
                             pch=19))
y <- n.fathers_0_9 %>% filter(focal==0)
set.seed(99); with(y, points(jitter(as.numeric(Treatment)+0.1,0.2),
                             standardized.n.father,
                             col=adjustcolor(alpha.f=0.6,col="grey"),
                             pch=1))
with(newdat_b, points(Treatment,fit,
                    pch=19,
                    cex=1.5,
                    col="black"))
with(newdat_b,segments(as.numeric(Treatment),lwr,as.numeric(Treatment),upr))

# Figure 4c
offspring.proportion.plot(x=offspring.proportion_0_9)
# Ignore error message, it just means there was no data for a certain position





# # Figure 5
# # Calculate the frequency of potential and observed distances between parents
brks <- seq(0,200,0.1)
# # Potential distances
mat <- dist(dat[,which(names(dat)%in%c("X","Y"))],diag=T)
Possible_distances <- mat[lower.tri(mat,diag=T)]
a <- hist(Possible_distances,breaks=brks,plot=F)
possible_density <- data.frame(density = a$density / max(a$density),
                               mids = a$mids)
# Observed
density_freq <- function(d){
  b <- hist(d$Distance,breaks=brks,plot=F)
  data.frame(n = b$counts,
             density = b$density / max(b$density),
             mids = b$mids)
}
# observed_density_0_1 <- density_freq(Paternity_dat_0_1)
# observed_density_0_5 <- density_freq(Paternity_dat_0_5)
observed_density_0_9 <- density_freq(Paternity_dat_0_9)

length(OffspringGenotypes) # 619 offspring genotyped
nrow(Paternity_dat_0_9) # 18 (18/619)*100 = 2.91%) were assigned paternity
# nrow(Paternity_dat_0_1) # 12 (12/619)*100
# nrow(Paternity_dat_0_5) # 14 (14/619)*100

with(Paternity_dat_0_9,table(Distance))
(11/18)*100 # 

Paternity_dat_0_9 %>% filter(Distance>139)

dd1=Paternity_dat_0_1[,1:2]
dd2=Paternity_dat_0_5[,1:2]
dd3=Paternity_dat_0_9[,1:2]

cbind.data.frame(dd1[,-1],dd3 %>% filter((OffspringID %in% dd1$OffspringID))) # same fathers
dd3 %>% filter(!(OffspringID %in% dd1$OffspringID))  # 6 offspring not assigned when father prob was 0.1
# (look at the Paternity.xls spreadsheet for scrutiny of genotypes)

# Function to create plot of observed sperm dispersal distances
make_plot <- function(y0,y1,xlims,xby){
  plot(c(0,xlims),c(0,1),type="n",bty="l",ylab="",xlab="",las=1,xaxt="n",yaxt="n")
  axis(side=1,at=seq(0,xlims,xby),cex.axis=1)
  axis(side=2,at=seq(0,1,0.2),cex.axis=1,las=1)
  y0 <- y0 %>% filter(density>0)
  y1 <- y1 %>% filter(density>0)
  offset = 0
  with(y0, segments(mids,rep(0,length(mids)),mids,density,lend=2,lwd=8,
                    col=adjustcolor(alpha.f=0.6,"lightgrey")))
  with(y1, segments(mids+offset,rep(0,length(mids)),mids+offset,density,lend=2,lwd=8,
                    col=adjustcolor(alpha.f=0.6,"tomato")))
}



quartz(width=4,height=3)
par(mfrow=c(2,1),mar=c(2,1,1.5,2),oma=c(2,3,0,0))
# a)
make_plot(y0=possible_density,y1=observed_density_0_9,xlims=155,xby=10)
mtext(side=3,"a) Whole array",adj=0)
# b) within blocks
make_plot(y0=possible_density,y1=observed_density_0_9,xlims=1.2,xby=0.1)
mtext(side=3,"b) Within treatments",adj=0)
# arrows(0,1.6,0,0.8,length=0.1,xpd=T,col="darkgrey")
# arrows(0.6,1.6,9,0.4,length=0.1,xpd=T,col="darkgrey")

mtext("Distance between mother and father (meters)",
      side=1,line=0.5,cex=1,outer=T,adj=0.5)
mtext("Relative frequency", side=2,line=2,cex=1,outer=T)
par(lend=2);legend(0.8,1.4,legend=c("Observed","Potential"),col=c("tomato","lightgrey"),lwd=2,bty="n",cex=0.8,xpd=T)
