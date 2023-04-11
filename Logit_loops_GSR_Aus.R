library(lme4)
library(ggplot2)
library(labeling)
library(ggeffects)
library(patchwork)

setwd("~/Desktop/Kinematics Paper/R things")
filepath = paste(getwd(),"loop-data.csv", sep="/")
beh.data <- read.csv(filepath)
beh.data <- beh.data[!is.na(beh.data$gscm_rel),]

beh.data$loop = as.factor(beh.data$loop)

cp.data <- beh.data[beh.data$genus == "Chrysopelea",]
dp.data <- beh.data[beh.data$genus == "Dendrelaphis",]

cp.logit <- glmer(loop ~ gscm_rel + (1|ID), data = cp.data, family = "binomial")
summary(cp.logit)

m <- max(cp.data$gscm_rel)
n <- min(dp.data$gscm_rel)

#dot plot of each individual's gap sizes and whether they did a loop or not
p2<- ggplot(dp.data, aes(x=gscm_rel,y=svl,color=loop))+geom_point()+xlim(n,m)+theme_minimal()
p3<- ggplot(cp.data, aes(x=gscm_rel,y=ID,color=loop))+geom_point()+xlim(n,m) + theme_minimal()

#as a histogram, instead - separate graphs for each species
p4 <- ggplot(cp.data, aes(x=gscm_rel, fill=loop)) + 
  geom_histogram(binwidth = 5, position = position_stack(reverse = TRUE))

p5 <- ggplot(dp.data, aes(x=gscm_rel, fill=loop)) + 
  geom_histogram(binwidth = 5, position = position_stack(reverse = TRUE))


p#as a histogram, with border indicating species
p6 <- ggplot(beh.data, aes(x=gscm_rel, color=genus, fill=loop)) + 
  geom_histogram(binwidth = 5, position = position_stack(reverse = TRUE))

p6 + 
  scale_fill_manual(values = c("black","white")) + 
  scale_color_manual(values = c("#207345","#0d576e"))


### IQRs
tmpdata <- cp.data
jvalues <- with(cp.data, seq(from=min(gscm_rel), to=max(gscm_rel), length.out=100))
pp <- lapply(jvalues, function(j) {
  tmpdata$gscm_rel <-j 
  predict(cp.logit, newdata = tmpdata, type = "response")
  })

#get mean PPs and IQ range for specific gap sizes
plotdat <- t(sapply(pp, function(x) {
  c(M = mean(x), quantile(x, c(0.25, 0.75)))
  }))

plotdat <- as.data.frame(cbind(plotdat, jvalues))
colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "GapSize")
head(plotdat)

#plot average marginal predicted probabilities
p1<- ggplot(plotdat, aes(x = GapSize, y = PredictedProbability)) + 
  geom_linerange(aes(ymin = Lower, ymax = Upper)) + 
  geom_line(size = 1) + ylim(c(0,1))+
  xlim(n,m)+theme_minimal()

#put the names of the plots you want to show below
p4
