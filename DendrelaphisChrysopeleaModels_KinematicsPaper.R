library(lme4)
library(ggplot2)

set.seed(500) #needed to ensure bootstraps are replicable

#load the data ####
setwd("~/Desktop/Kinematics Paper/R Things") 
filepath = paste(getwd(),"all_data_for_stats.csv", sep="/")
cdata <- read.csv(filepath)

#fix units: all should be as frac SVL (not %SVL, which gap size currently is)
cdata$gscm_rel <- cdata$gscm_rel/100
#center all gap size data around the lowest value of gscm_rel
dmin = min(cdata$gscm_rel)
cdata$gscm_relA = cdata$gscm_rel - dmin

#separate Chrysopelea and Dendrelaphis data
df<-cdata[!(cdata$genus=="Chrysopelea"),]
cf <-cdata[!(cdata$genus=="Dendrelaphis"),]

#drop snakes 13 and 15: too little data
df <- df[!(df$ID == "2019DP13"),]
df <- df[!(df$ID == "2019DP15"),]

#drop large gap size data in CP
dmax = max(df$gscm_rel)
cf <- cf[(cf$gscm_rel<dmax),]

#MODEL 1: Relative horizontal variation in Dendrelaphis####
#Decide random effects structure: slopes vs intercepts
hv.lmer <- lmer(hv_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE)
hv.lmer1<- lmer(hv_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(hv.lmer,hv.lmer1)
#hv.lmer1 has lower AIC.

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(hv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(hv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(hv.lmer1)

summary(hv.lmer1)
confint(hv.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#to see conditional variances by snake:
ranef(hv.lmer1)
#to see fixed effect coefficients by snake:
coef(hv.lmer1)

#plot
ggplot(df,aes(x = gscm_relA, y = hv_rel, colour = ID)) +
  geom_point() + geom_smooth(method = "lm", fill = NA) +facet_wrap( ~ svl)

#MODEL 2: Relative vertical variation  ####
#Decide random effects structure: slopes vs intercepts
zv.lmer <- lmer(zv_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE)
zv.lmer1<- lmer(zv_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(zv.lmer,zv.lmer1)
#zv.lmer has lower AIC.

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(zv.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(zv.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(zv.lmer)

summary(zv.lmer)
confint(zv.lmer, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#to see conditional variances by snake:
ranef(zv.lmer)
#to see fixed effect coefficients by snake:
coef(zv.lmer)

#plot
ggplot(df,aes(x = gscm_relA, y = zv_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

#MODEL 3: Avg velocity  ####
#Decide random effects structure: slopes vs intercepts
av.lmer <- lmer(av_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE)
av.lmer1<- lmer(av_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(av.lmer,av.lmer1)
#av.lmer has lower AIC.

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(av.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(av.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(av.lmer)

summary(av.lmer)
confint(av.lmer, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#to see conditional variances by snake:
ranef(av.lmer)
#to see fixed effect coefficients by snake:
coef(av.lmer)

#plot
ggplot(df,aes(x = gscm_relA, y = av_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl) 


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

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(mv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(mv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(mv.lmer1)

summary(mv.lmer3)
confint(mv.lmer3, method = "boot", nsim = 1500) # The number of bootstrap iterations = 1500.
#running additional iterations because some do not converge. Result is > 1000 iterations used.

#to see conditional variances by snake:
ranef(mv.lmer3)
#to see fixed effect coefficients by snake:
coef(mv.lmer3)

#plot
ggplot(df,aes(x = gscm_relA, y = mv_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

#MODEL 5: Landing velocity ####

#Decide random effects structure: slopes vs intercepts
lv.lmer <- lmer(lv_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
lv.lmer1 <- lmer(lv_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(lv.lmer,lv.lmer1) #lv.lmer1 has lower AIC

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(lv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(lv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(lv.lmer1)

summary(lv.lmer1)
confint(lv.lmer1, method = "boot", nsim = 1500) # The number of bootstrap iterations = 1500.
#running additional iterations because some do not converge. Result is > 1000 iterations used.

#to see conditional variances by snake:
ranef(lv.lmer1)
#to see fixed effect coefficients by snake:
coef(lv.lmer)

#plot
ggplot(df,aes(x = gscm_relA, y = lv_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

#MODEL 6: Loop depth  ####
#Decide random effects structure: slopes vs intercepts
ld.lmer <- lmer(LD_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
#check for singularity
tt <- getME(ld.lmer,"theta")
ll <- getME(ld.lmer,"lower")
min(tt[ll==0])  #not close to 0
ss <- getME(ld.lmer,c("theta","fixef"))
ld.lmer1 <- update(ld.lmer,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4))) #still does not converge
ld.lmer2 <- update(ld.lmer,start=ss,control=lmerControl(optimizer="bobyqa",
                                                        optCtrl=list(maxfun=2e5))) #singular
ld.lmer3 <- lmer(LD_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(ld.lmer3, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(ld.lmer3, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(ld.lmer3)

#plot
ggplot(df,aes(x = gscm_relA, y = LD_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)
ggplot(df,aes(x = gscm_relA, y = LD_rel, colour = ID)) +
  geom_point() 

#MODEL 7: Arc height  ####
#Decide random effects structure: slopes vs intercepts
ah.lmer <- lmer(AH_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
ah.lmer1 <- lmer(AH_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(ah.lmer,ah.lmer1) #ah.lmer has lower AIC

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(ah.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(ah.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(ah.lmer)

#plot
ggplot(df,aes(x = gscm_relA, y = AH_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

# MODEL 8: Z position of head at AF ####
#Decide random effects structure: slopes vs intercepts
zpaf.lmer <- lmer(zpaf_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
#check for singularity
tt <- getME(zpaf.lmer,"theta")
ll <- getME(zpaf.lmer,"lower")
min(tt[ll==0])  
ss <- getME(zpaf.lmer,c("theta","fixef"))
zpaf.lmer1 <- update(zpaf.lmer,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))
zpaf.lmer2<- lmer(zpaf_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(zpaf.lmer1,zpaf.lmer2) #zpaf.lmer1 has lower AIC

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(zpaf.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(zpaf.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(zpaf.lmer1)

summary(zpaf.lmer1)
confint(zpaf.lmer1, method = "boot", nsim = 1500) # The number of bootstrap iterations = 1500.
#running additional iterations because some do not converge. Result is > 1000 iterations used.

#to see conditional variances by snake:
ranef(zpaf.lmer1)
#to see fixed effect coefficients by snake:
coef(zpaf.lmer1)

#plot
ggplot(df,aes(x = gscm_relA, y = zpaf_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

# MODEL 9: Z position Max ####
#Decide random effects structure: slopes vs intercepts
zpmx.lmer <- lmer(zpmx_rel ~ gscm_relA + svl + (gscm_relA|ID), data = df, REML = TRUE) 
zpmx.lmer2<- lmer(zpmx_rel ~ gscm_relA + svl + (1|ID), data = df, REML = TRUE)
anova(zpmx.lmer,zpmx.lmer2) #zpmx.lmer has lower AIC

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(zpmx.lmer, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(zpmx.lmer, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(zpmx.lmer)

summary(zpmx.lmer)
confint(zpmx.lmer, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000.

#to see conditional variances by snake:
ranef(zpmx.lmer)
#to see fixed effect coefficients by snake:
coef(zpmx.lmer)

#plot
ggplot(df,aes(x = gscm_relA, y = zpmx_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)




# CHRYSOPELEA PLOTS ####
# Model 1a: Chrysopelea, horizontal extent ####
# decide RE structure
cp.hv.lmer <- lmer(hv_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #singular
cp.hv.lmer1 <- lmer(hv_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.hv.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.hv.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(cp.hv.lmer1)

summary(cp.hv.lmer1)
confint(cp.hv.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#to see conditional variances by snake:
ranef(cp.hv.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.hv.lmer1)

ggplot(cf,aes(x = gscm_relA, y = hv_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

# Model 2a: Chrysopelea, vertical extent ####
#Decide random effects structure: slopes vs intercepts
cp.zv.lmer <- lmer(zv_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #does not converge. 
#check for singularity
tt <- getME(cp.zv.lmer,"theta")
ll <- getME(cp.zv.lmer,"lower")
min(tt[ll==0])  #not close to 0
ss <- getME(cp.zv.lmer,c("theta","fixef"))
cp.zv.lmer1 <- update(cp.zv.lmer,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4))) #conveges
cp.zv.lmer2<- lmer(zv_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)
anova(cp.zv.lmer1,cp.zv.lmer2)
#cp.zv.lmer2 has lower AIC.

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.zv.lmer2, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.zv.lmer2, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(cp.zv.lmer2)

summary(cp.zv.lmer1)
confint(cp.zv.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#to see conditional variances by snake:
ranef(cp.zv.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.zv.lmer1)

#plot
ggplot(cf,aes(x = gscm_relA, y = zv_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

#Model 8a: Chrysopelea Z position of head at AF ####
#Decide random effects structure: slopes vs intercepts
cp.zpaf.lmer <- lmer(zpaf_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) 
cp.zpaf.lmer1<- lmer(zpaf_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)
anova(cp.zpaf.lmer1,cp.zpaf.lmer)
#same AIC< using simpler model (cp.zpaf.lmer1)

#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.zpaf.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.zpaf.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(cp.zpaf.lmer1)

summary(cp.zpaf.lmer1)
confint(cp.zpaf.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#to see conditional variances by snake:
ranef(cp.zpaf.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.zpaf.lmer1)

#plot
ggplot(cf,aes(x = gscm_relA, y = zpaf_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)

#Model 9a: Chrysopelea, Z position of head (max) ####
cp.zpmx.lmer <- lmer(zpmx_rel ~ gscm_relA + (gscm_relA|ID), data = cf, REML = TRUE) #singular
cp.zpmx.lmer1<- lmer(zpmx_rel ~ gscm_relA + (1|ID), data = cf, REML = TRUE)


#Diagnostics:
#fitted vs residual; tests assumption of linearity of residuals. 
plot(cp.zpmx.lmer1, type=c("p","smooth"), col.line=1)
#scale-response plot; tests assumption of homoskedasticity
plot(cp.zpmx.lmer1, sqrt(abs(resid(.))) ~ fitted(.),type = c("p", "smooth"))
#qqplot; test assumption of normality
lattice::qqmath(cp.zpmx.lmer1)
#looks reasonable.

summary(cp.zpmx.lmer1)
confint(cp.zpmx.lmer1, method = "boot", nsim = 1000) # The number of bootstrap iterations = 1000. 

#to see conditional variances by snake:
ranef(cp.zpmx.lmer1)
#to see fixed effect coefficients by snake:
coef(cp.zpmx.lmer1)

#plot
ggplot(cf,aes(x = gscm_relA, y = zpmx_rel, colour = ID)) +
  geom_point() +facet_wrap( ~ svl)


