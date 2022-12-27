library(lme4)
library(ggplot2)

#load the data ####
setwd("~/Desktop/Kinematics Paper/R Things") 
filepath = paste(getwd(),"all_data_for_stats.csv", sep="/")
cdata <- read.csv(filepath)

#fix units: all should be as frac SVL (not %SVL, which gap size currently is)
cdata$gscm_rel <- cdata$gscm_rel/100

#center all gap size data around the lowest value of gscm_rel 
dmin = min(cdata$gscm_rel)
cdata$gscm_relA = cdata$gscm_rel - dmin

#center all body size data around average body size of Chrysopelea so intercepts will be comparable
csvl_avg = mean(85,71.5,66.7,95.9,72.5,72.5)
cdata$svl = cdata$svl - csvl_avg 


#separate Chrysopelea and Dendrelaphis data
df<-cdata[!(cdata$genus=="Chrysopelea"),]
cf <-cdata[!(cdata$genus=="Dendrelaphis"),]

#drop snakes 13 and 15: too little data
df <- df[!(df$ID == "2019DP13"),]
df <- df[!(df$ID == "2019DP15"),]

#drop large gap size data in CP
dmax = max(df$gscm_rel)
cf <- cf[(cf$gscm_rel<dmax),]

#MODELS
#MODEL 1: Relative horizontal variation in Dendrelaphis####
#Decide random effects structure: slopes vs intercepts
hv.lmer <- lmer(hv_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE)
hv.lmer1<- lmer(hv_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(hv.lmer,hv.lmer1)
#hv.lmer1 has lower AIC.
summary(hv.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(hv.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#MODEL 2: Relative vertical variation  ####
#Decide random effects structure: slopes vs intercepts
zv.lmer <- lmer(zv_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE)
zv.lmer1<- lmer(zv_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(zv.lmer,zv.lmer1)
#zv.lmer has lower AIC.

summary(zv.lmer)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(zv.lmer, method = "boot", nsim = 1800) # Increased # iterations to make sure getting large number without warnings/singularity

#MODEL 3: Avg velocity  ####
#Decide random effects structure: slopes vs intercepts
av.lmer <- lmer(av_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE)
av.lmer1<- lmer(av_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(av.lmer,av.lmer1)
#av.lmer has lower AIC.

summary(av.lmer)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(av.lmer, method = "boot", nsim = 1800) #inc. number to ensure > 1000 no warnings

#MODEL 4: Max velocity  ####
#Decide random effects structure: slopes vs intercepts
mv.lmer <- lmer(mv_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) #does not converge
#check for singularity
tt <- getME(mv.lmer,"theta")
ll <- getME(mv.lmer,"lower")
min(tt[ll==0])  #somewhat close to 0, model may become singular.
ss <- getME(mv.lmer,c("theta","fixef"))
mv.lmer1 <- update(mv.lmer,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4))) #does not converge
mv.lmer2 <- update(mv.lmer,start=ss,control=lmerControl(optimizer="bobyqa",
                                                        optCtrl=list(maxfun=2e5))) #singular

mv.lmer3 <- lmer(mv_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)


summary(mv.lmer3)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(mv.lmer3, method = "boot", nsim = 1500) # The number of bootstrap iterations = 1500.
#running additional iterations because some do not converge. Result is > 1000 iterations used.

#MODEL 5: Landing velocity ####

#Decide random effects structure: slopes vs intercepts
lv.lmer <- lmer(lv_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
lv.lmer1 <- lmer(lv_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(lv.lmer,lv.lmer1) #lv.lmer1 has lower AIC


summary(lv.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(lv.lmer1, method = "boot", nsim = 1500) # The number of bootstrap iterations = 1500.
#running additional iterations because some do not converge. Result is > 1000 iterations used.

#MODEL 6: Loop depth  ####
#Decide random effects structure: slopes vs intercepts
ld.lmer <- lmer(LD_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
#check for singularity
tt <- getME(ld.lmer,"theta")
ll <- getME(ld.lmer,"lower")
min(tt[ll==0])  #not close to 0
ss <- getME(ld.lmer,c("theta","fixef"))
ld.lmer1 <- update(ld.lmer,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4))) #singular
ld.lmer2 <- lmer(LD_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)

summary(ld.lmer2)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(ld.lmer2, method="boot",nsim=1500)

#MODEL 7: Arc height  ####
#Decide random effects structure: slopes vs intercepts
ah.lmer <- lmer(AH_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
ah.lmer1 <- lmer(AH_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(ah.lmer,ah.lmer1) #ah.lmer has lower AIC

summary(ah.lmer)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(ah.lmer, method = "boot", nsim = 3500) # The number of bootstrap iterations = 3500.
#running additional iterations because some do not converge. Result is > 1000 iterations used.

# MODEL 8: Z position of head at AF ####
#Decide random effects structure: slopes vs intercepts
zpaf.lmer <- lmer(zpaf_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE)  #singular
zpaf.lmer1 <- lmer(zpaf_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)

summary(zpaf.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(zpaf.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000.

# MODEL 9: Z position Max ####
#Decide random effects structure: slopes vs intercepts
zpmx.lmer <- lmer(zpmx_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
zpmx.lmer1<- lmer(zpmx_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(zpmx.lmer,zpmx.lmer1) #zpmx.lmer has lower AIC

summary(zpmx.lmer)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(zpmx.lmer, method = "boot", nsim = 2000) # The number of bootstrap iterations = 2000.


#DIAGNOSTICS ####
#Deviations from normality not an issue due to only using bootstrap methods for hypothesis testing

#fitted vs residual; tests assumption of linearity of residuals. 
plot(hv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(hv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(zv.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(zv.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(av.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(av.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(mv.lmer3, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(mv.lmer3, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(lv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(lv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(ld.lmer2, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(ld.lmer2, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(ah.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(ah.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(zpaf.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(zpaf.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(zpmx.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(zpmx.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))



#COEFFICIENTS BY INDIVIDUAL ####
#to see conditional variances by snake:
ranef(hv.lmer1)
#to see fixed effect coefficients by snake:
coef(hv.lmer1)

#to see conditional variances by snake:
ranef(zv.lmer)
#to see fixed effect coefficients by snake:
coef(zv.lmer)

#to see conditional variances by snake:
ranef(av.lmer)
#to see fixed effect coefficients by snake:
coef(av.lmer)

#to see conditional variances by snake:
ranef(mv.lmer3)
#to see fixed effect coefficients by snake:
coef(mv.lmer3)

#to see conditional variances by snake:
ranef(lv.lmer1)
#to see fixed effect coefficients by snake:
coef(lv.lmer1)

#to see conditional variances by snake:
ranef(ld.lmer2)
#to see fixed effect coefficients by snake:
coef(ld.lmer2)

#to see conditional variances by snake:
ranef(ah.lmer)
#to see fixed effect coefficients by snake:
coef(ah.lmer)

#to see conditional variances by snake:
ranef(zpaf.lmer1)
#to see fixed effect coefficients by snake:
coef(zpaf.lmer1)

#to see conditional variances by snake:
ranef(zpmx.lmer)
#to see fixed effect coefficients by snake:
coef(zpmx.lmer)

# CHRYSOPELEA DATA ####
# Model 1a: Chrysopelea, horizontal extent ####
# decide RE structure
cp.hv.lmer <- lmer(hv_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #singular
cp.hv.lmer1 <- lmer(hv_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)

summary(cp.hv.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.hv.lmer1, method = "boot", nsim = 1500) 

# Model 2a: Chrysopelea, vertical extent ####
#Decide random effects structure: slopes vs intercepts
cp.zv.lmer <- lmer(zv_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE)  
cp.zv.lmer1 <- lmer(zv_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)
anova(cp.zv.lmer,cp.zv.lmer1)
#cp.zv.lmer1 has lower AIC.

summary(cp.zv.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.zv.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

# Model 3a: Chrysopelea, average head speed ####
# decide RE structure
cp.av.lmer <- lmer(av_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE)
cp.av.lmer1 <- lmer(av_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)
anova(cp.av.lmer,cp.av.lmer1) #lmer1 has lower AIC

summary(cp.av.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.av.lmer1, method = "boot", nsim = 1500) # The number of bootstrap iterations = 1500. 

# Model 4a: Chrysopelea, max head speed ####
# decide RE structure
cp.mv.lmer <- lmer(mv_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #singular
cp.mv.lmer1 <- lmer(mv_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)

summary(cp.mv.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.mv.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1500. 

# Model 5a: Chrysopelea, landing head speed ####
# decide RE structure
cp.lv.lmer <- lmer(lv_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #singular
cp.lv.lmer1 <- lmer(lv_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)

summary(cp.lv.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.lv.lmer1, method = "boot", nsim = 2500) 

# Model 6a: Chrysopelea, loop depth####
# decide RE structure
cp.ld.lmer <- lmer(LD_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE)
cp.ld.lmer1 <- lmer(LD_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)
anova(cp.ld.lmer,cp.ld.lmer1) #lmer has lower AIC; using random slopes

summary(cp.ld.lmer)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.ld.lmer, method = "boot", nsim = 2000)

# Model 7a: Chrysopelea, arc height ####
# decide RE structure
cp.ah.lmer <- lmer(AH_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #degenerate
cp.ah.lmer1 <- lmer(AH_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)

summary(cp.ah.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.ah.lmer1, method = "boot", nsim = 1500) 

#Model 8a: Chrysopelea Z position of head at AF ####
#Decide random effects structure: slopes vs intercepts
cp.zpaf.lmer <- lmer(zpaf_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) 
cp.zpaf.lmer1<- lmer(zpaf_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)
anova(cp.zpaf.lmer1,cp.zpaf.lmer)
#same AIC< using simpler model (cp.zpaf.lmer1)

summary(cp.zpaf.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.zpaf.lmer1, method = "boot", nsim = 1500)

#Model 9a: Chrysopelea, Z position of head (max) ####
cp.zpmx.lmer <- lmer(zpmx_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #singular
cp.zpmx.lmer1<- lmer(zpmx_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)

summary(cp.zpmx.lmer1)
set.seed(500) #needed to ensure bootstraps are replicable 
confint(cp.zpmx.lmer1, method = "boot", nsim = 1500) 

#Coefficients/variances by individual ####

#to see conditional variances by snake:
ranef(cp.hv.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.hv.lmer1)

#to see conditional variances by snake:
ranef(cp.zv.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.zv.lmer1)

#to see conditional variances by snake:
ranef(cp.av.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.av.lmer1)

#to see conditional variances by snake:
ranef(cp.mv.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.mv.lmer1)

#to see conditional variances by snake:
ranef(cp.lv.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.lv.lmer1)

#to see conditional variances by snake:
ranef(cp.ld.lmer)
#to see fixed effect coefficients by snake:
coef(cp.ld.lmer)

#to see conditional variances by snake:
ranef(cp.ah.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.ah.lmer1)

#to see conditional variances by snake:
ranef(cp.zpaf.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.zpaf.lmer1)

#to see conditional variances by snake:
ranef(cp.zpmx.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.zpmx.lmer1)

#DIAGNOSTICS####

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.hv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.hv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.zv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.zv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.av.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.av.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.mv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.mv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.lv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.lv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.ld.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.ld.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.ah.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.ah.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.zpaf.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.zpaf.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))

#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.zpmx.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.zpmx.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
