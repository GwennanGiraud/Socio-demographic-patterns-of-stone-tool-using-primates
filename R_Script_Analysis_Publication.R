#########################Script R Analyses###############################

rm(list=ls())

#####Import libraries#####
library(dplyr)
library(lme4)
library(car)
library(DHARMa)
library(ggplot2)
library(sjPlot)
library(ANTs)
library(glmmTMB)
library(performance)
library(reshape2)
library(pbapply)
library(patchwork)

#####Dataset#####
###Use the Environment_Analysis provided (with 42 individuals)
load("~/Environment_Analysis_Published.RData")

#########Dataset preparation#################
###Transform variables
str(Dataset)
Dataset$Age.Class<-as.factor(Dataset$Age.Class)
Dataset$Age.Group<-as.factor(Dataset$Age.Group)
Dataset$Sex<-as.factor(Dataset$Sex)
Dataset$Hair.pattern<-as.factor(Dataset$Hair.pattern)
Dataset$Tool.user<-as.factor(Dataset$Tool.user)
str(Dataset)

###Standardize rank variable
table(Dataset$Sex)
Dataset<-Dataset %>% 
  mutate(RankStand=ifelse(Sex=="F", Rank/26, Rank/23)*10) %>% 
  relocate(RankStand, .after = Rank)

###Change the reference level
Dataset$Age.Class <- relevel(Dataset$Age.Class, ref = "Juvenile")



##############Analysis###################

####Socio-demographic Model####
#Is Rank and Sex associated with tool user status?
table(Dataset$Tool.user, Dataset$Age.Class) #All subadults are tool users, so we merged subadults with adults for the analyses.
table(Dataset$Tool.user, Dataset$Sex)

###Merging subadults with adults
Dataset.mature<-Dataset%>%
  mutate(Age.Class = ifelse(Age.Class == "Subadult", "Adult", as.character(Age.Class)),
         Age.Class = as.factor(Age.Class))
table(Dataset.mature$Age.Class)

###Change the reference level
Dataset.mature$Age.Class <- relevel(Dataset.mature$Age.Class, ref = "Juvenile")

###Tool users as response variable----------------
#Full model with interaction Rank*Age
##Check for multi-collinearity
ToolUserfull0<-glm(Tool.user ~ RankStand*Age.Class + Sex + Hair.pattern, family=binomial(link="logit"), data=Dataset.mature)

ToolUserfull1<-glm(Tool.user ~ RankStand + Age.Class + Sex + Hair.pattern, family=binomial(link="logit"), data=Dataset.mature)
vif(ToolUserfull1, type ='terms') #<5 ok!

##Full/Null model comparison
mod0<-glm(Tool.user ~ 1, family=binomial(link="logit"), data = Dataset.mature)
#LRT
lrt<-anova(mod0,ToolUserfull0,test="Chisq")
lrt #S OK the full model is meaningful

mod1<-glm(Tool.user ~ Age.Class + Hair.pattern, family=binomial(link="logit"), data = Dataset.mature)
lrt<-anova(mod1,ToolUserfull0,test="Chisq")
lrt #S OK Test predictors (Rank and Sex) affect the model

##Check normality of model residuals
hist(residuals(ToolUserfull0))
qqnorm(residuals(ToolUserfull0))
qqline(residuals(ToolUserfull0))
Obj<-simulateResiduals(ToolUserfull0)
plot(Obj) #NS! OK

##Model stability
my.dfbeta<-function(m){
  xx=cbind(coef(m), coef(m)+ t(apply(X=dfbeta(m), MARGIN=2, FUN=range)))
  colnames(xx)=c("orig", "min", "max")
  return(xx)
}
my.dfbeta(m=ToolUserfull0) #Looks good (little variation)
source("C:/Users/Gwenn/Desktop/Revisions/GitHub/diagnostic_fcns.r")
m.stab.plot(my.dfbeta(m=ToolUserfull0)) 
summary(ToolUserfull0) #Age is not significant + Age is a control variable so it is OK

#Better plot
m.stab.plot<-function(est, lower=NULL, upper=NULL, xnames=NULL, col="black", center.at.null=F){
  if(ncol(est)==3){
    lower=est[, 2]
    upper=est[, 3]
    xnames=rownames(est)
    est=est[, 1]
    x.at=est
  }
  old.par = par(no.readonly = TRUE)
  par(mar=c(5, 0.5, 0.5, 0.5), mgp=c(3.5, 1.4, 0), tcl=-0.2)
  plot(x=est, y=1:length(est), pch=18, xlab="Estimate", ylab="", yaxt="n", 
       xlim=range(c(lower, upper)), type="n", ylim=c(1, length(est)+1), 
       cex.lab = 2.4, cex.axis = 2)
  abline(v=0, lty=3)
  if(center.at.null){
    x.at=rep(0, length(est))
    text(labels=xnames, x=x.at, y=(1:length(est))+0.3, cex=2) 
  }else{
    text(labels=xnames, x=x.at, y=(1:length(est))+0.3, cex=2, pos=c(2, 4)[1+as.numeric(est<0)]) 
  }
  points(x=est, y=1:length(est), pch=18, col=col, cex=2) 
  segments(x0=lower, x1=upper, y0=1:length(est), y1=1:length(est), col=col, lwd=2) 
  par(old.par)
}
Dataset.mature2<-Dataset.mature %>% 
  rename(
    Phenotype = Hair.pattern,
    Rank.brut = Rank,
    Rank = RankStand,
    Age=Age.Class
  )
ToolUserfullPLOT<-glm(Tool.user ~ Rank*Age + Sex + Phenotype, family=binomial(link="logit"), data=Dataset.mature2)
m.stab.plot(my.dfbeta(m=ToolUserfullPLOT))

###Check for influential cases
max(abs(dffits(ToolUserfull0))) #1.27 (not bigger than  2, so OK)
hist(dffits(ToolUserfull0))
## cooks distance
max(cooks.distance(ToolUserfull0))  ## (0.26 = OK, threshold = 1)
## leverage
l.thresh=lev.thresh(ToolUserfull0)    
l.thresh ##  Give a value which will be the threshold: Here 0.33.
max(as.vector(influence(ToolUserfull0)$hat)) ## Here 0.43, goes over the threshold = influential cases.
#But no specific reasons to remove individuals here.

###Reduced model
##With backwards stepwise elimination
#Effect of interactions
ToolUserfull0red<-glm(Tool.user ~ RankStand + Age.Class + Sex + Hair.pattern, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova1<-anova (ToolUserfull0red, ToolUserfull0, test="Chisq") #NS, then the interaction do not bring much to the model.
ResultsAnova1

#Effect of rank
ToolUserfull0red2<-glm(Tool.user ~ Age.Class + Sex + Hair.pattern, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova2<-anova (ToolUserfull0red2, ToolUserfull0, test="Chisq") #NS, then Rank and interaction do not bring much to the model.
ResultsAnova2

#Effect of Age
ToolUserfull0red3<-glm(Tool.user ~ RankStand + Sex + Hair.pattern, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova3<-anova (ToolUserfull0red3, ToolUserfull0, test="Chisq")  #NS, then Age and interaction do not bring much to the model.
ResultsAnova3

#Effect of Sex
ToolUserfull0red4<-glm(Tool.user ~ RankStand*Age.Class + Hair.pattern, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova4<-anova (ToolUserfull0red4, ToolUserfull0, test="Chisq") #S, then Sex is meaningful for the model.
ResultsAnova4

#Effect of phenotype
ToolUserfull0red5<-glm(Tool.user ~ RankStand*Age.Class + Sex, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova5<-anova (ToolUserfull0red5, ToolUserfull0, test="Chisq") #S, then Phenotype is meaningful for the model.
ResultsAnova5

#Results model
ResultsTU<-summary(ToolUserfull0)
#Reverse logit function
Reverselogit.Estm<-plogis(coef(ResultsTU)[,1])
Reverselogit.SE<-plogis(coef(ResultsTU)[,2])
CoefResultsTU<-summary(ToolUserfull0)$coefficients
ResultsTU<-cbind(CoefResultsTU, Reverselogit.Estm)
ResultsTU<-cbind(ResultsTU, Reverselogit.SE)
ResultsTU 


###Interpretation and plots
#Phenotype: Hybrid-like phenotype are more likely tool users than Common-like phenotype
Plot = ggplot(Dataset.mature, aes(x = Tool.user, fill=Hair.pattern)) + geom_bar(width = 0.6)
Plot = Plot+labs(x= 'Tool user', y = 'Number of individuals', fill="Phenotype")
#Plot = Plot+ggtitle("Hybrid-like phenotypes are \n more likely to be tool users")+theme(plot.title = element_text(size = 32))
Plot = Plot+theme(axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30))
Plot = Plot+theme(axis.text.x = element_text(size=26), axis.text.y = element_text(size=24))
Plot = Plot+theme(legend.title=element_text(size=26), legend.text = element_text(size=24))
Plot = Plot+scale_fill_manual(values = c("Common-like" = "#969696", "Hybrid-like" = "#ffd92f"))
Plot

#Sex: Males are more likely tool users than females
#Rename categories for plots
Dataset.mature$Sex<-as.character(Dataset.mature$Sex)
Dataset.mature$Sex=ifelse(Dataset.mature$Sex %in% "Female","Females", Dataset.mature$Sex)
Dataset.mature$Sex=ifelse(Dataset.mature$Sex %in% "Male","Males", Dataset.mature$Sex)
Dataset.mature$Sex<-as.factor(Dataset.mature$Sex)

Plot2 = ggplot(Dataset.mature, aes(x = Tool.user, fill=Sex)) + geom_bar(width = 0.6)
Plot2 = Plot2+labs(x= 'Tool user', y = 'Number of individuals', fill="Sex category")
#Plot2 = Plot2+ggtitle("Males are more likely to be tool users")+theme(plot.title = element_text(size = 22))
Plot2 = Plot2+theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
Plot2 = Plot2+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
Plot2 = Plot2+theme(legend.title=element_text(size=16), legend.text = element_text(size=14))
Plot2 = Plot2+scale_fill_manual(values = c("Females" = "#984ea3", "Males" = "#4daf4a"))
Plot2

#Rank*Age: NS
Plot3<-plot_model(ToolUserfull0, type = "pred", terms = c("RankStand[all]", "Age.Class"), jitter=T, color=c("blue","red","green"), show.data = TRUE)
Plot3<-Plot3+labs(title = "Tool user status according to rank and age")
Plot3<-Plot3+labs(x= 'Hierarchical rank', y = 'Tool user status')+theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
Plot3<-Plot3+theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
Plot3<-Plot3+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
Plot3



###Tool users as response variable----------------
#Full model WITHOUT interaction Rank*Age, to test rank effect alone
##Check for multi-collinearity
ToolUser0<-glm(Tool.user ~ RankStand + Sex + Hair.pattern + Age.Class, family=binomial(link="logit"), data=Dataset.mature)
vif(ToolUser0, type ='terms') #<5 ok!

##Full/Null model comparison
mod0<-glm(Tool.user ~ 1, family=binomial(link="logit"), data = Dataset.mature)
#LRT
lrt<-anova(mod0,ToolUser0,test="Chisq")
lrt #S OK the full model is meaningful

mod1<-glm(Tool.user ~ Hair.pattern + Age.Class, family=binomial(link="logit"), data = Dataset.mature)
lrt<-anova(mod1,ToolUser0,test="Chisq")
lrt #S OK Test predictors (Rank and Sex) affect the model

##Check normality of model residuals
hist(residuals(ToolUser0))
qqnorm(residuals(ToolUser0))
qqline(residuals(ToolUser0))
Obj<-simulateResiduals(ToolUser0)
plot(Obj) #NS! OK

###Model stability
my.dfbeta<-function(m){
  xx=cbind(coef(m), coef(m)+ t(apply(X=dfbeta(m), MARGIN=2, FUN=range)))
  colnames(xx)=c("orig", "min", "max")
  return(xx)
}
my.dfbeta(m=ToolUser0) #Looks good (little variation)
source("C:/Users/Gwenn/Desktop/Revisions/GitHub/diagnostic_fcns.r")
m.stab.plot(my.dfbeta(m=ToolUser0)) #OK
summary(ToolUser0)

#Better plot
ToolUserfullPLOT2<-glm(Tool.user ~ Rank + Sex + Phenotype + Age, family=binomial(link="logit"), data=Dataset.mature2)
m.stab.plot(my.dfbeta(m=ToolUserfullPLOT2))

###Check for influential cases
## max absolute change (full model): dffits 
max(abs(dffits(ToolUser0))) #0.69 (not bigger than  2)
hist(dffits(ToolUser0)) 
## cooks distance
max(cooks.distance(ToolUser0))  ## (0.12 = OK, threshold = 1)
## leverage 
l.thresh=lev.thresh(ToolUser0)    
l.thresh ##  Give a value which will be the threshold: Here 0.29.
max(as.vector(influence(ToolUser0)$hat)) ## Here 0.24, does not go over the threshold = no influential cases.

###Reduced model
##With backwards stepwise elimination
#Effect of rank
ToolUser0red1<-glm(Tool.user ~ Sex + Hair.pattern + Age.Class, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova1<-anova (ToolUser0red1, ToolUser0, test="Chisq") #NS, then Rank do not bring much to the model.
ResultsAnova1

#Effect of Age
ToolUser0red2<-glm(Tool.user ~ RankStand + Sex + Hair.pattern, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova2<-anova (ToolUser0red2, ToolUser0, test="Chisq")  #NS, then Age do not bring much to the model.
ResultsAnova2

#Effect of Sex
ToolUser0red3<-glm(Tool.user ~ RankStand + Hair.pattern + Age.Class, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova3<-anova (ToolUser0red3, ToolUser0, test="Chisq") #S, then Sex is meaningful for the model.
ResultsAnova3

#Effect of phenotype
ToolUser0red4<-glm(Tool.user ~ RankStand + Sex + Age.Class, family=binomial(link="logit"), data=Dataset.mature)
ResultsAnova4<-anova (ToolUser0red4, ToolUser0, test="Chisq") #S, then Phenotype is meaningful for the model.
ResultsAnova4

#Results model
ResultsTU2<-summary(ToolUser0)
#Reverse logit function
Reverselogit.Estm<-plogis(coef(ResultsTU2)[,1])
Reverselogit.SE<-plogis(coef(ResultsTU2)[,2])
CoefResultsTU2<-summary(ToolUser0)$coefficients
ResultsTU2<-cbind(CoefResultsTU2, Reverselogit.Estm)
ResultsTU2<-cbind(ResultsTU2, Reverselogit.SE)
ResultsTU2

#Rank*Age: NS
Plot4<-plot_model(ToolUser0, type = "pred", terms = "RankStand[all]", jitter=T, show.data = TRUE)
Plot4<-Plot4+labs(title = "Tool user status according to rank")
Plot4<-Plot4+labs(x= 'Hierarchical rank', y = 'Tool user status')+theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
Plot4<-Plot4+theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
Plot4<-Plot4+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
Plot4



###Rank as response variable----------------Just for sanity sake
#Full model
##Check for multi-collinearity
Rankfull<-glm(RankStand ~ Tool.user + Hair.pattern + Age.Class, family=poisson(link="log"), data=Dataset.mature)
vif(Rankfull) #<5 ok!

##Full/Null model comparison
mod0<-glm(RankStand ~ 1, family=poisson(link="log"), data = Dataset.mature)
#LRT
lrt<-anova(mod0,Rankfull,test="Chisq")
lrt #S the full model is meaningful

mod1<-glm(RankStand ~ Age.Class, family=poisson(link="log"), data = Dataset.mature)
lrt<-anova(mod1,Rankfull,test="Chisq")
lrt #NS Test predictors (Tool.user and Hair.pattern) do not affect the model

##Check normality of model residuals
hist(residuals(Rankfull)) 
qqnorm(residuals(Rankfull))
qqline(residuals(Rankfull)) 
Obj<-simulateResiduals(Rankfull)
plot(Obj) #NS! OK

###Model stability
my.dfbeta<-function(m){
  xx=cbind(coef(m), coef(m)+ t(apply(X=dfbeta(m), MARGIN=2, FUN=range)))
  colnames(xx)=c("orig", "min", "max")
  return(xx)
}
my.dfbeta(m=Rankfull) #Looks good (little variation)
source("C:/Users/Gwenn/Desktop/Revisions/GitHub/diagnostic_fcns.r")
m.stab.plot(my.dfbeta(m=Rankfull))

###Check for influential cases
## max absolute change (full model): dffits 
max(abs(dffits(Rankfull))) #0.79 (not bigger than  2)
hist(dffits(Rankfull))
## cooks distance
max(cooks.distance(Rankfull))  ## (0.27 = OK, threshold = 1)
## leverage 
l.thresh=lev.thresh(Rankfull)    
l.thresh ##  Give a value which will be the threshold: Here 0.29.
max(as.vector(influence(Rankfull)$hat)) ## Here 0.21, does not go over the threshold = no influential cases.

#Model
GLMRank = glm(RankStand ~ Tool.user + Hair.pattern + Age.Class, family=poisson(link="log"), data=Dataset.mature)
summary(GLMRank) #S for Age.

#Plot
table(Dataset.mature$Tool.user, Dataset.mature$Age.Class)
plot(RankStand~Age.Class, data=Dataset.mature)



################Grooming metrics as response variable################
#Are grooming metrics associated with tool user status?
#Literature: Age has an effect on grooming interactions. Here we use the dataset with all age categories.

###Social Position model on grooming network: EIGENVECTOR------------
#Full model
EigenGroomfull0<-glm(eigenG ~ Tool.user + Hair.pattern + RankStand*Sex + Age.Class, family=Gamma(link="log"), data=Dataset)

##Check for multi-collinearity without interactions
EigenGroomfull<-glm(eigenG ~ Tool.user + Hair.pattern + RankStand + Sex + Age.Class, family=Gamma(link="log"), data=Dataset)
vif(EigenGroomfull) #<5 OK!

##Full/Null model comparison
mod0<-glm(eigenG ~ 1, family=Gamma(link="log"), data = Dataset)
#LRT
lrt<-anova(mod0,EigenGroomfull0,test="Chisq")
lrt #S OK Model is meaningful

mod1<-glm(eigenG ~ RankStand*Sex + Age.Class, family=Gamma(link="log"), data = Dataset)
#LRT
lrt<-anova(mod1,EigenGroomfull0,test="Chisq")
lrt #NS Test predictors do not affect the model

##Check normality of model residuals
hist(residuals(EigenGroomfull0)) 
qqnorm(residuals(EigenGroomfull0))
qqline(residuals(EigenGroomfull0))
Obj<-simulateResiduals(EigenGroomfull0)
plot(Obj) #Outlier S pb but KS and Dispersion tests NS, need to check outliers

###Model stability
my.dfbeta<-function(m){
  xx=cbind(coef(m), coef(m)+ t(apply(X=dfbeta(m), MARGIN=2, FUN=range)))
  colnames(xx)=c("orig", "min", "max")
  return(xx)
}
my.dfbeta(m=EigenGroomfull0) #Looks good (little variation)
source("C:/Users/Gwenn/Desktop/Revisions/GitHub/diagnostic_fcns.r")
m.stab.plot(my.dfbeta(m=EigenGroomfull0))
summary(EigenGroomfull0) #But Tool use and Age are not significant + Age is a control variable so it is OK

###Check for influential cases
## max absolute change (full model): dffits 
max(abs(dffits(EigenGroomfull0))) #1.56 (not bigger than  2)
hist(dffits(EigenGroomfull0))
## cooks distance
max(cooks.distance(EigenGroomfull0))  ## (0.25 = OK, threshold = 1)
## leverage 
l.thresh=lev.thresh(EigenGroomfull0)    
l.thresh ##  Give a value which will be the threshold: Here 0.43.
max(as.vector(influence(EigenGroomfull0)$hat)) ## Here 0.38, does not go over the threshold = no influential cases.

#Permutations on the full model
perm = perm.net.nl(Dataset, labels = c("eigenG"), nperm = 10000)
GLMEigenGroom = stat.glm(perm, formula = eigenG ~ Tool.user + Hair.pattern + RankStand*Sex + Age.Class, family=Gamma(link="log"))
ant(GLMEigenGroom)$model.diagnostic
glmEigenGroom = ant(GLMEigenGroom)
glmEigenGroom

#Reverse log function
Reverselog.Estm<-exp(glmEigenGroom$model$coefficients$Estimate)
Reverselog.SE<-exp(glmEigenGroom$model$coefficients$`Std. Error`)
ResultsEigenGroom<-glmEigenGroom$model$coefficients
ResultsEigenGroom<-cbind(ResultsEigenGroom, Reverselog.Estm)
ResultsEigenGroom<-cbind(ResultsEigenGroom, Reverselog.SE)
ResultsEigenGroom


###Plots
#Rename categories for plots
Dataset2<-Dataset
Dataset2$Sex<-as.character(Dataset2$Sex)
Dataset2$Sex=ifelse(Dataset2$Sex %in% "Female","Females", Dataset2$Sex)
Dataset2$Sex=ifelse(Dataset2$Sex %in% "Male","Males", Dataset2$Sex)
Dataset2$Sex<-as.factor(Dataset2$Sex)

EigenGroomfull<-glm(eigenG ~ Tool.user + Hair.pattern + RankStand*Sex + Age.Class, family=Gamma(link="log"), data=Dataset2)
Plot5<-plot_model(EigenGroomfull, type = "pred", terms = c("RankStand [all]", "Tool.user", "Sex"), jitter=T, show.data = TRUE) +
  scale_color_manual(name = "Tool User", values = c("Yes" = "#FFA500", "No" = "blue")) + 
  scale_fill_manual(name = "Tool User", values = c("Yes" = "#FFA500", "No" = "blue")) +
  labs(title = "Predicted probabilities of being in a central \n position within the grooming network", 
       x= 'Hierarchical rank', y = 'Centrality position \n in grooming network')+
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by=2))
Plot5

##Plot Grooming Network Sociogram
#See R_Script_Grooming_Sociogram_Publication.



###Grooming received model: INStrength------------
#Full model
INStrGroomfull0<-glmmTMB(instrengthG ~ Tool.user + Hair.pattern + RankStand*Sex + Age.Class, family = tweedie(link = "log"), data=Dataset)

##Check for multi-collinearity without interactions
INStrGroomfull<-glmmTMB(instrengthG ~ Tool.user + Hair.pattern + RankStand + Sex + Age.Class, family = tweedie(link = "log"), data=Dataset)
check_collinearity(INStrGroomfull) #<5 OK!

##Full/Null model comparison
mod0<-glmmTMB(instrengthG ~ 1, family = tweedie(link = "log"), data = Dataset)
#LRT
lrt<-anova(mod0,INStrGroomfull0,test="Chisq")
lrt #NS!

#Cannot conclude



###Grooming Partners model: INDegree------------
#Full model
INDegGroomfull0<-glm(indegreeG ~ Tool.user + Hair.pattern + RankStand*Sex + Age.Class, family=poisson(link="log"), data=Dataset)

##Check for multi-collinearity without interactions
INDegGroomfull<-glm(indegreeG ~ Tool.user + Hair.pattern + RankStand + Sex + Age.Class, family=poisson(link="log"), data=Dataset)
vif(INDegGroomfull) #<5 OK!

##Full/Null model comparison
mod0<-glm(indegreeG ~ 1, family=poisson(link="log"), data = Dataset)
#LRT
lrt<-anova(mod0,INDegGroomfull0,test="Chisq")
lrt #S model meaningful

mod1<-glm(indegreeG ~ RankStand*Sex + Age.Class, family=poisson(link="log"), data = Dataset)
#LRT
lrt<-anova(mod1,INDegGroomfull0,test="Chisq")
lrt #NS Test predictors do not affect the model

##Check normality of model residuals
hist(residuals(INDegGroomfull0))
qqnorm(residuals(INDegGroomfull0))
qqline(residuals(INDegGroomfull0))
Obj<-simulateResiduals(INDegGroomfull0)
plot(Obj) #NS OK!

###Model stability
my.dfbeta<-function(m){
  xx=cbind(coef(m), coef(m)+ t(apply(X=dfbeta(m), MARGIN=2, FUN=range)))
  colnames(xx)=c("orig", "min", "max")
  return(xx)
}
my.dfbeta(m=INDegGroomfull0) #Looks good (little variation)
source("C:/Users/Gwenn/Desktop/Revisions/GitHub/diagnostic_fcns.r")
m.stab.plot(my.dfbeta(m=INDegGroomfull0))

###Check for influential cases
## max absolute change (full model): dffits
max(abs(dffits(INDegGroomfull0))) #1.20 (not bigger than  2)
hist(dffits(INDegGroomfull0))
## cooks distance
max(cooks.distance(INDegGroomfull0))  ## (0.28 = OK, threshold = 1)
## leverage
l.thresh=lev.thresh(INDegGroomfull0)    
l.thresh ##  Give a value which will be the threshold: Here 0.43.
max(as.vector(influence(INDegGroomfull0)$hat)) ## Here 0.60, goes over the threshold = influential cases.
#But no specific reasons to remove individuals.

#Permutations on the full model
perm = perm.net.nl(Dataset, labels = c("indegreeG"), nperm = 10000)
GLMIndegGroom = stat.glm(perm, formula = indegreeG ~ Tool.user + Hair.pattern + RankStand*Sex + Age.Class, family=poisson(link="log"))
ant(GLMIndegGroom)$model.diagnostic
glmIndegGroom = ant(GLMIndegGroom)
glmIndegGroom

#Reverse log function
Reverselog.Estm<-exp(glmIndegGroom$model$coefficients$Estimate)
Reverselog.SE<-exp(glmIndegGroom$model$coefficients$`Std. Error`)
ResultsIndegGroom<-glmIndegGroom$model$coefficients
ResultsIndegGroom<-cbind(ResultsIndegGroom, Reverselog.Estm)
ResultsIndegGroom<-cbind(ResultsIndegGroom, Reverselog.SE)
ResultsIndegGroom


###Interpretation and plots
#Males received grooming from less partners than females.
INDegGroomfull<-glm(indegreeG ~ Tool.user + Hair.pattern + RankStand*Sex + Age.Class, family=poisson(link="log"), data=Dataset2)
Plot6 = plot_model(INDegGroomfull, type = "pred", terms = c("RankStand [all]", "Tool.user", "Sex"), jitter=T, show.data = TRUE) +
  scale_color_manual(name = "Tool User", values = c("Yes" = "#FFA500", "No" = "blue")) + 
  scale_fill_manual(name = "Tool User", values = c("Yes" = "#FFA500", "No" = "blue")) + 
  labs(title = "Hierarchical rank had an opposite \n influence on females and males", 
       x= 'Hierarchical rank', y = 'Predicted number of individuals \n emitting grooming')+
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by=2))
Plot6



#######HOMOPHILY ACCORDING TO TOOL USE STATUS IN GROOMING NETWORK##########
#Is there a phenomenon of homophily when it comes to the use of tools?

###GROUP LEVEL PATTERNS OF GROOMING--------------------

###SEE R_Script_to_prepare_for_Python

###Then Python_Script_for_homophily_correction
#Corrections based on paper of Karimi & Oliveira 2023

###Results
##Mixing matrix (normalized)
#         NTU         TU
#NTU [[0.11242604 0.22781065] 
#TU  [0.22781065 0.43195266]]

##Adjusted mixing matrix
#         NTU         TU
#NTU [[0.05749041 0.12610909]
#TU  [0.10687838 0.22088422]]

##Adjusted assortativity with groups of unequal sizes (correction 1)
Adj_assort <- 0.010066249673468076 #=> close to 0; so no homophily
##Assessing asymmetric mixing patterns in networks (correction 2)
Asym_patt <- -0.013260728235921436 #=> close to 0; so no asymmetric mixing patterns


###SUBGROUP LEVEL PATTERNS OF GROOMING--------------------
###ASSOCIATIONS DEPENDING ON TOOL USE STATUS (interaction type) and other factors (rank, sex, age)

###Prepare data for ego networks creation
TU = Dataset$Tool.user == 'Yes' #Filter for tool users
NTU = Dataset$Tool.user == 'No' #Filter for non-tool users
d = ifelse(Dataset$Tool.user == 'Yes', T, F) #Convert to logical
f = ifelse(Dataset$Tool.user == 'No', T, F) 
N_TU = sum(d) #Number of tool users
N_NTU = nrow(Dataset) - N_TU #Number of non-tool users
alters_info = d #Logical vector for tool user status
alters_info2 = f #Logical vector for non-tool user status

###Filter Dataset to get only tool users ids and loop only through them
dfTU = Dataset[TU,] #Subset Dataset for tool users only
dfNTU = Dataset[NTU,] #Subset Dataset for non-tool users only
idsTU = dfTU$id #Get ids of tool users
idsNTU = dfNTU$id #Get ids of non-tool users

#For TU
rTU = NULL #Initialize result data frame
for (i in 1:length(idsTU)) { #Loop through each tool user
  egoTU = dfTU$id[i] #Get the current tool user's id
  
  ##In matrix mGroom
  # Get out- and in-strength (number or weight of interactions)
  outSt = mGroom[colnames(mGroom) == egoTU, ] #Get the out-strength for the current tool user
  inSt = mGroom[,colnames(mGroom) == egoTU] #Get the in-strength for the current tool user
  St = outSt + inSt #Calculate the strength of ties for the current tool user
  
  # Get degrees as binary presence/absence
  outDeg = (outSt > 0) * 1 #Get the out-degree for the current tool user
  inDeg = (inSt > 0) * 1 #Get the in-degree for the current tool user
  Deg = inDeg + outDeg #Calculate the total degree for the current tool user
  
  # Sum of ties/degree with TU vs NTU alters
  St_with_TU = sum(St[alters_info]) #Sum the strength of ties with tool users
  St_with_NTU = sum(St[!alters_info]) #Sum the strength of ties with non-tool users
  Deg_with_TU = sum(Deg[alters_info]) #Sum the degree with tool users
  Deg_with_NTU = sum(Deg[!alters_info]) #Sum the degree with non-tool users
  
  #Define interaction types and group sizes
  interaction_type = c('TU', 'NTU') #Define interaction types
  grp_size = c(N_TU, N_NTU) #Group sizes for tool users and non-tool users
  
  #Bind results into a data frame
  rTU = rbind(rTU, data.frame(egoTU, 'Strength' = c(St_with_TU, St_with_NTU), 
                              'Degree' = c(Deg_with_TU, Deg_with_NTU), 
                              interaction_type, grp_size))
}
rTU #Store results in a separate variable for tool users

#For NTU
rNTU = NULL #Initialize result data frame
for (i in 1:length(idsNTU)) { #Loop through each tool user
  egoNTU = dfNTU$id[i] #Get the current tool user's id
  
  ##In matrix mGroom
  # Get out- and in-strength (number or weight of interactions)
  outSt = mGroom[colnames(mGroom) == egoNTU, ] #Get the out-strength for the current tool user
  inSt = mGroom[,colnames(mGroom) == egoNTU] #Get the in-strength for the current tool user
  St = outSt + inSt #Calculate the strength of ties for the current tool user
  
  # Get degrees as binary presence/absence
  outDeg = (outSt > 0) * 1 #Get the out-degree for the current tool user
  inDeg = (inSt > 0) * 1 #Get the in-degree for the current tool user
  Deg = inDeg + outDeg #Calculate the total degree for the current tool user
  
  # Sum of ties/degree with TU vs NTU alters
  St_with_TU = sum(St[alters_info2]) #Sum the strength of ties with tool users
  St_with_NTU = sum(St[!alters_info2]) #Sum the strength of ties with non-tool users
  Deg_with_TU = sum(Deg[alters_info2]) #Sum the degree with tool users
  Deg_with_NTU = sum(Deg[!alters_info2]) #Sum the degree with non-tool users
  
  #Define interaction types and group sizes
  interaction_type = c('TU', 'NTU') #Define interaction types
  grp_size = c(N_TU, N_NTU) #Group sizes for tool users and non-tool users
  
  #Bind results into a data frame
  rNTU = rbind(rNTU, data.frame(egoNTU, 'Strength' = c(St_with_TU, St_with_NTU), 
                                'Degree' = c(Deg_with_TU, Deg_with_NTU), 
                                interaction_type, grp_size))
}
rNTU #Store results in a separate variable for non-tool users

###Add information for individuals
#For TU
rTU_merged <- merge(rTU, Dataset, by.x = "egoTU", by.y = "id")
rTU <- rTU_merged[,c("name", "Age.Class", "Sex", "Rank", "RankStand", "Hair.pattern", "Tool.user",
                     "egoTU", "Strength", "Degree", "interaction_type", "grp_size")]
rTU

#For NTU
rNTU_merged <- merge(rNTU, Dataset, by.x = "egoNTU", by.y = "id")
rNTU <- rNTU_merged[,c("name", "Age.Class", "Sex", "Rank", "RankStand", "Hair.pattern", "Tool.user",
                       "egoNTU", "Strength", "Degree", "interaction_type", "grp_size")]
rNTU

#Convert interaction_type to factor for modeling
rTU$interaction_type<-as.factor(rTU$interaction_type)
rNTU$interaction_type<-as.factor(rNTU$interaction_type)



###ASSOCIATION MODELS###

###Association Duration Model: STRENGTH for TU---
##Full model
model_Str_TU <- glmmTMB(
  Strength ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoTU),
  ziformula = ~1, 
  family = tweedie(link = "log"), 
  data = rTU
)

##Check for multi-collinearity without interactions
model_Str_TU2<-glmmTMB(
  Strength ~ interaction_type + RankStand + Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoTU),
  ziformula = ~1, 
  family = tweedie(link = "log"),  
  data = rTU
)
check_collinearity(model_Str_TU2) #<5 OK

##Full/Null model comparison
mod0<-glmmTMB(Strength ~ 1 + offset(log(grp_size)) + (1 | egoTU), ziformula = ~1, 
              family = tweedie(link = "log"),  data = rTU)
#LRT
lrt<-anova(mod0,model_Str_TU,test="Chisq")
lrt #S ok Model is meaningful

mod1<-glmmTMB(Strength ~ RankStand + Sex + Age.Class 
              + offset(log(grp_size)) + (1 | egoTU), ziformula = ~1, 
              family = tweedie(link = "log"),  data = rTU)
#LRT
lrt<-anova(mod1,model_Str_TU,test="Chisq")
lrt #NS Test predictor (interaction_type) does not affect the model

##Check normality of model residuals
hist(residuals(model_Str_TU))
qqnorm(residuals(model_Str_TU))
qqline(residuals(model_Str_TU))
Obj<-simulateResiduals(model_Str_TU)
plot(Obj) #NS OK!

###Model stability
leave_one_out_effects <- function(model, data) {
  original_coef <- fixef(model)$cond
  n <- nrow(data)
  delta <- matrix(NA, nrow = n, ncol = length(original_coef))
  colnames(delta) <- names(original_coef)
  
  cat("Start of the leave-one-out...\n")
  
  for (i in 1:n) {
    dat_i <- data[-i, ]
    
    # Offset recalculation if necessary
    offset_i <- log(dat_i$grp_size)
    
    fit_i <- tryCatch(
      glmmTMB(
        formula = formula(model),
        data = dat_i,
        family = eval(model$call$family),
        ziformula = eval(model$call$ziformula),
        offset = offset_i,
        REML = FALSE
      ),
      error = function(e) {
        message(sprintf("❌ EObservation  error %d : %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(fit_i)) {
      coef_i <- fixef(fit_i)$cond
      delta[i, ] <- coef_i - original_coef
    }
    
    if (i %% 10 == 0) cat(sprintf("… %d/%d Achived\n", i, n))
  }
  
  cat("Achived.\n")
  return(as.data.frame(delta))
}

# Apply function
delta <- leave_one_out_effects(model = model_Str_TU, data = rTU)
# Clean errors
delta_clean <- delta[complete.cases(delta), , drop = FALSE]
#Output similar to my.dfbeta:
# Coefficient original
coef_original <- fixef(model_Str_TU)$cond
# Min et max values calculation after leave-one-out
dfbeta_like <- data.frame(
  orig = coef_original,
  min = coef_original + apply(delta_clean, 2, min, na.rm = TRUE),
  max = coef_original + apply(delta_clean, 2, max, na.rm = TRUE)
)
round(dfbeta_like, 4)

#Recalculation of dfbetas
# 'delta_clean' already : estimate_leave_one_out - estimate_original
dfbeta_equiv <- delta_clean
#Zero centered for predictors
dfbeta_equiv_centered <- scale(dfbeta_equiv, center = TRUE, scale = FALSE)
dfbeta_equiv_centered <- as.data.frame(dfbeta_equiv_centered)
# Add id column
dfbeta_equiv_centered$ObsID <- 1:nrow(dfbeta_equiv_centered)
# Long format
dfbeta_long <- melt(dfbeta_equiv_centered, id.vars = "ObsID",
                    variable.name = "Predictor", value.name = "DFBETA")

# Plot visualization
ggplot(dfbeta_long, aes(x = DFBETA, y = Predictor)) +
  geom_point(alpha = 0.6, color = "skyblue") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", color = "black", fatten = 1.2) +
  theme_minimal() +
  labs(
    title = "Stability plot (DFBETA approach)",
    x = "Centered Leave-one-out effect",
    y = ""
  )

# Recall full model
model_Str_TU <- glmmTMB(
  Strength ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoTU),
  ziformula = ~1, 
  family = tweedie(link = "log"), 
  data = rTU
)

# Estimate extraction of interaction_typeTU
effect_TU <- summary(model_Str_TU)$coefficients$cond["interaction_typeTU", c("Estimate", "Std. Error")]
obs_estimate <- effect_TU["Estimate"] 
obs_se <- effect_TU["Std. Error"] 


# Function for data permutation and model adjustment
permute_and_fit <- function(data) {
  permuted_data <- do.call(rbind, lapply(split(data, data$egoTU), function(subdf) { # Split data by egoTU
    # Subset TU and NTU interaction types
    interaction_types <- subdf$interaction_type # Get unique interaction types for this egoTU
    subdf$interaction_type <- sample(interaction_types)  # permute within ego
    return(subdf) # Permuted subdf
  }))
  
  # Re-do model
  model_perm <- try(glmmTMB(
    Strength ~ interaction_type + RankStand * Sex + Age.Class + offset(log(grp_size)) + 
      (1 | egoTU),
    ziformula = ~1, family = tweedie(link = "log"),
    data = permuted_data
  ), silent = TRUE) # Permutation model fitting
  
  if (inherits(model_perm, "try-error")) return(NA) # If model fitting fails, return NA
  
  # Extract the estimate for interaction_typeTU
  est <- try(summary(model_perm)$coefficients$cond["interaction_typeTU", "Estimate"], silent = TRUE) 
  if (inherits(est, "try-error")) return(NA) # If extraction fails, return NA
  
  return(est) # Return the estimate
}

#Permutation
set.seed(123)
n_perm <- 10000
data_perm <- rTU

#perm_estimates <- replicate(n_perm, permute_and_fit(data_perm)) # Permutation without progress bar
#perm_estimates <- pbreplicate(n_perm, permute_and_fit(data_perm)) #with progression bar
perm_estimates <- pbapply::pbreplicate(n_perm, permute_and_fit(data_perm)) #same but faster

#Empiric P-value
perm_estimates <- perm_estimates[!is.na(perm_estimates)] # Suppress eventual errors (NA)

# bilateral p-value
p_empirical <- mean(abs(perm_estimates) >= abs(obs_estimate))
p_empirical 

# Results with interaction_typeTU
cat("Estimate =", round(obs_estimate, 3), "| SE =", round(obs_se, 3), 
    "| p-value =", round(p_empirical, 3), "\n")
Assoc_results_StrengthTU <- "Estimate = 0.24 | SE = 0.186 | p-value = 0.177"

ResultModel<-summary(model_Str_TU)
#Reverse log function
Reverselog.Estm<-exp(ResultModel$coefficients$cond[, "Estimate"])
Reverselog.SE<-exp(ResultModel$coefficients$cond[, "Std. Error"])
ResultStr<-ResultModel$coefficients$cond
ResultStr<-cbind(ResultStr, Reverselog.Estm)
ResultStr<-cbind(ResultStr, Reverselog.SE)
ResultStrTU<-ResultStr
ResultStrTU

# Visualisation
hist(perm_estimates, breaks = 50, main = "Permutate estimates distribution \n for interaction_typeTU")
abline(v = obs_estimate, col = "red", lwd = 2) # Add line of observed estimate



##Association Duration Model: STRENGTH for NTU---
# Full model
model_Str_NTU <- glmmTMB(
  Strength ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1, 
  family = tweedie(link = "log"), 
  data = rNTU
)

##Check for multi-collinearity without interactions
model_Str_NTU2<-glmmTMB(
  Strength ~ interaction_type + RankStand + Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1,
  family = tweedie(link = "log"), 
  data = rNTU
)

check_collinearity(model_Str_NTU2) #<5 OK

##Full/Null model comparison
mod0<-glmmTMB(Strength ~ 1 + offset(log(grp_size)) + (1 | egoNTU), ziformula = ~1, 
              family = tweedie(link = "log"), data = rNTU) 
#LRT
lrt<-anova(mod0,model_Str_NTU,test="Chisq")
lrt 

mod1<-glmmTMB(Strength ~ RankStand * Sex + Age.Class
              + offset(log(grp_size)) + (1 | egoNTU), ziformula = ~1, 
              family = tweedie(link = "log"), data = rNTU) 
#LRT
lrt<-anova(mod1,model_Str_NTU,test="Chisq")
lrt #S donc test predictor (interaction_type) affect the model

##Check normality of model residuals
hist(residuals(model_Str_NTU)) 
qqnorm(residuals(model_Str_NTU))
qqline(residuals(model_Str_NTU))
Obj<-simulateResiduals(model_Str_NTU)
plot(Obj) #NS OK!

###Model stability
leave_one_out_effects <- function(model, data) {
  original_coef <- fixef(model)$cond
  n <- nrow(data)
  delta <- matrix(NA, nrow = n, ncol = length(original_coef))
  colnames(delta) <- names(original_coef)
  
  cat("Starting of leave-one-out...\n")
  
  for (i in 1:n) {
    dat_i <- data[-i, ]
    
    # Offset recalculation if necessary
    offset_i <- log(dat_i$grp_size)
    
    fit_i <- tryCatch(
      glmmTMB(
        formula = formula(model),
        data = dat_i,
        family = eval(model$call$family),
        ziformula = eval(model$call$ziformula),
        offset = offset_i,
        REML = FALSE
      ),
      error = function(e) {
        message(sprintf("❌ EObservation error %d : %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(fit_i)) {
      coef_i <- fixef(fit_i)$cond
      delta[i, ] <- coef_i - original_coef
    }
    
    if (i %% 10 == 0) cat(sprintf("… %d/%d achived\n", i, n))
  }
  
  cat("Achived.\n")
  return(as.data.frame(delta))
}

# Apply function
delta <- leave_one_out_effects(model = model_Str_NTU, data = rNTU)
# Clean lines with errors
delta_clean <- delta[complete.cases(delta), , drop = FALSE]
# Original coefficient
coef_original <- fixef(model_Str_NTU)$cond
# Calculation of min et max values after leave-one-out
dfbeta_like <- data.frame(
  orig = coef_original,
  min = coef_original + apply(delta_clean, 2, min, na.rm = TRUE),
  max = coef_original + apply(delta_clean, 2, max, na.rm = TRUE)
)
round(dfbeta_like, 4)

#Recalculation of dfbetas
# 'delta_clean' already : estimate_leave_one_out - estimate_original
dfbeta_equiv <- delta_clean
#Zero centered around predictors
dfbeta_equiv_centered <- scale(dfbeta_equiv, center = TRUE, scale = FALSE)
dfbeta_equiv_centered <- as.data.frame(dfbeta_equiv_centered)

# Add ID column
dfbeta_equiv_centered$ObsID <- 1:nrow(dfbeta_equiv_centered)

# Long format transformation
dfbeta_long <- melt(dfbeta_equiv_centered, id.vars = "ObsID",
                    variable.name = "Predictor", value.name = "DFBETA")

# Plot
ggplot(dfbeta_long, aes(x = DFBETA, y = Predictor)) +
  geom_point(alpha = 0.6, color = "skyblue") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", color = "black", fatten = 1.2) +
  theme_minimal() +
  labs(
    title = "Stability plot (DFBETA approach)",
    x = "Centered leave-one-out effect",
    y = ""
  )

# Recall full model
model_Str_NTU <- glmmTMB(
  Strength ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1, family = tweedie(link = "log"),  
  data = rNTU
)

# Extraction estimate of interaction_typeTU
effet_NTU <- summary(model_Str_NTU)$coefficients$cond["interaction_typeTU", c("Estimate", "Std. Error")]
obs_estimate <- effet_NTU["Estimate"] 
obs_se <- effet_NTU["Std. Error"] 


# Function for data permutation and model adjustment
permute_and_fit <- function(data) {
  permuted_data <- do.call(rbind, lapply(split(data, data$egoNTU), function(subdf) { # Split data by egoTU
    # Subset TU and NTU interaction types
    interaction_types <- subdf$interaction_type # Get unique interaction types for this egoTU
    subdf$interaction_type <- sample(interaction_types)  # permute within ego
    return(subdf) # Permuted subdf
  }))
  
  # Re-do model
  model_perm <- try(glmmTMB(
    Strength ~ interaction_type + RankStand * Sex + Age.Class + offset(log(grp_size)) + 
      (1 | egoNTU),
    ziformula = ~1, family = tweedie(link = "log"),
    data = permuted_data
  ), silent = TRUE) # Permutation model fitting
  
  if (inherits(model_perm, "try-error")) return(NA) # If model fitting fails, return NA
  
  # Extract the estimate for interaction_typeTU
  est <- try(summary(model_perm)$coefficients$cond["interaction_typeTU", "Estimate"], silent = TRUE) 
  if (inherits(est, "try-error")) return(NA) # If extraction fails, return NA
  
  return(est) # Return the estimate
}

#Permutation
set.seed(123) 
n_perm <- 10000
data_perm <- rNTU

#perm_estimates <- replicate(n_perm, permute_and_fit(data_perm)) # Permutation without progress bar
#perm_estimates <- pbreplicate(n_perm, permute_and_fit(data_perm)) #with progression bar
perm_estimates <- pbapply::pbreplicate(n_perm, permute_and_fit(data_perm)) #same but faster

#Empiric p-value
perm_estimates <- perm_estimates_vec[!is.na(perm_estimates_vec)] # Suppress potential errors (NA)

# Bilateral p-value
p_empirical <- mean(abs(perm_estimates_vec) >= abs(obs_estimate))
p_empirical 

# Results of interaction_typeTU
cat("Estimate =", round(obs_estimate, 3), "| SE =", round(obs_se, 3), 
    "| p-value =", round(p_empirical, 3), "\n")
Assoc_results_StrengthNTU <- "Estimate = -1.096 | SE = 0.294 | p-value < 0.0001"

ResultModel<-summary(model_Str_NTU)
#Reverse log function
Reverselog.Estm<-exp(ResultModel$coefficients$cond[, "Estimate"])
Reverselog.SE<-exp(ResultModel$coefficients$cond[, "Std. Error"])
ResultStr<-ResultModel$coefficients$cond
ResultStr<-cbind(ResultStr, Reverselog.Estm)
ResultStr<-cbind(ResultStr, Reverselog.SE)
ResultStrNTU<-ResultStr
ResultStrNTU

# Visualisation
hist(perm_estimates, breaks = 50, main = "Permuted estimate distribution \n for interaction_typeTU")
abline(v = obs_estimate, col = "red", lwd = 2) # Add observed estimate line



##Association Partners Model: DEGREE for TU---
# Full model
model_Deg_TU <- glmmTMB(
  Degree ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoTU),
  ziformula = ~1, 
  family = poisson(link="log"),
  data = rTU)

hist(residuals(model_Deg_TU))

VarCorr(model_Deg_TU) # 0.27599: The random variance associated with egoTU is not zero, 
#which justifies keeping it in the model (it is not a useless artifact). 
#It does indeed help capture interindividual variance.

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
overdisp_fun(model_Deg_TU)
#Ratio 1.10 No real overdispersion, P-value NS

##Check for multi-collinearity without interactions
model_Deg_TU2<-glmmTMB( 
  Degree ~ interaction_type + RankStand + Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoTU),
  ziformula = ~1, 
  family = poisson(link = "log"), 
  data = rTU)

check_collinearity(model_Deg_TU2) #OK

##Full/Null model comparison
mod0<-glmmTMB(Degree ~ 1 + offset(log(grp_size)) + (1 | egoTU), ziformula = ~1, 
              family = poisson(link = "log"), data = rTU)
#LRT
lrt<-anova(mod0,model_Deg_TU,test="Chisq")
lrt #S ok Model meaningful

mod1<-glmmTMB(Degree ~ RankStand + Sex + Age.Class
              + offset(log(grp_size)) + (1 | egoTU), ziformula = ~1, 
              family = poisson(link = "log"), data = rTU)
#LRT
lrt<-anova(mod1,model_Deg_TU,test="Chisq")
lrt #NS Test predictor does not affect the model

##Check normality of model residuals
hist(residuals(model_Deg_TU))
qqnorm(residuals(model_Deg_TU))
qqline(residuals(model_Deg_TU)) 
Obj<-simulateResiduals(model_Deg_TU)
plot(Obj) #NS OK!

###Model stability
leave_one_out_effects <- function(model, data) {
  original_coef <- fixef(model)$cond
  n <- nrow(data)
  delta <- matrix(NA, nrow = n, ncol = length(original_coef))
  colnames(delta) <- names(original_coef)
  
  cat("Starting of leave-one-out...\n")
  
  for (i in 1:n) {
    dat_i <- data[-i, ]
    
    # Model refit without observation i
    fit_i <- tryCatch(
      glmmTMB(
        formula = formula(model),
        data = dat_i,
        family = eval(model$call$family),
        ziformula = if (!is.null(model$call$ziformula)) eval(model$call$ziformula) else ~0
      ),
      error = function(e) {
        message(sprintf("❌ EObservation error %d : %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(fit_i)) {
      coef_i <- fixef(fit_i)$cond
      delta[i, ] <- coef_i - original_coef
    }
    
    if (i %% 10 == 0) cat(sprintf("… %d/%d achived\n", i, n))
  }
  
  cat("Achived.\n")
  return(as.data.frame(delta))
}

delta <- leave_one_out_effects(model = model_Deg_TU, data = rTU)
delta_clean <- delta[complete.cases(delta), , drop = FALSE]

# Original coefficient
coef_original <- fixef(model_Deg_TU)$cond
# Calculation of min et max values after leave-one-out
dfbeta_like <- data.frame(
  orig = coef_original,
  min = coef_original + apply(delta_clean, 2, min, na.rm = TRUE),
  max = coef_original + apply(delta_clean, 2, max, na.rm = TRUE)
)
round(dfbeta_like, 4)

#dfbetas recalculation
# 'delta_clean' already : estimate_leave_one_out - estimate_original
dfbeta_equiv <- delta_clean
#Zero centered around predictors
dfbeta_equiv_centered <- scale(dfbeta_equiv, center = TRUE, scale = FALSE)
dfbeta_equiv_centered <- as.data.frame(dfbeta_equiv_centered)

# Add ID column
dfbeta_equiv_centered$ObsID <- 1:nrow(dfbeta_equiv_centered)

# Long format
dfbeta_long <- melt(dfbeta_equiv_centered, id.vars = "ObsID",
                    variable.name = "Predictor", value.name = "DFBETA")

# Plot
ggplot(dfbeta_long, aes(x = DFBETA, y = Predictor)) +
  geom_point(alpha = 0.6, color = "skyblue") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", color = "black", fatten = 1.2) +
  theme_minimal() +
  labs(
    title = "Stability plot (DFBETA approach)",
    x = "Centered leave-one-out effect",
    y = ""
  )


# Recall full model
model_Deg_TU <- glmmTMB(
  Degree ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoTU),
  ziformula = ~1, 
  family = poisson(link="log"),
  data = rTU)


# Extraction estimate of interaction_typeTU
effect_TU <- summary(model_Deg_TU)$coefficients$cond["interaction_typeTU", c("Estimate", "Std. Error")]
obs_estimate <- effect_TU["Estimate"]  
obs_se <- effect_TU["Std. Error"] 


# Function to permute data and adjust model
permute_and_fit <- function(data) {
  permuted_data <- do.call(rbind, lapply(split(data, data$egoTU), function(subdf) { # Split data by egoTU
    # Subset TU and NTU interaction types
    interaction_types <- subdf$interaction_type # Get unique interaction types for this egoTU
    subdf$interaction_type <- sample(interaction_types)  # permute within ego
    return(subdf) # Permuted subdf
  }))
  
  # Model refit
  model_perm <- try(glmmTMB(
    Degree ~ interaction_type + RankStand * Sex + Age.Class + offset(log(grp_size)) + 
      (1 | egoTU),
    ziformula = ~1,
    family = poisson(link = "log"),
    data = permuted_data
  ), silent = TRUE) # Permutation model fitting
  
  if (inherits(model_perm, "try-error")) return(NA) # If model fitting fails, return NA
  
  # Extract the estimate for interaction_typeTU
  est <- try(summary(model_perm)$coefficients$cond["interaction_typeTU", "Estimate"], silent = TRUE) 
  if (inherits(est, "try-error")) return(NA) # If extraction fails, return NA
  
  return(est) # Return the estimate
}

#Lancer la permutation
set.seed(123) 
n_perm <- 10000
data_perm <- rTU

#perm_estimates <- replicate(n_perm, permute_and_fit(data_perm)) # Permutation without progress bar
#perm_estimates <- pbreplicate(n_perm, permute_and_fit(data_perm)) #with progression bar
perm_estimates <- pbapply::pbreplicate(n_perm, permute_and_fit(data_perm)) #same but faster

#Empiric P-value
perm_estimates <- perm_estimates[!is.na(perm_estimates)] # Suppress potential errors (NA)

# Bilateral p-value
p_empirical <- mean(abs(perm_estimates) >= abs(obs_estimate))
p_empirical

# Results of interaction_typeTU
cat("Estimate =", round(obs_estimate, 3), "| SE =", round(obs_se, 3), 
    "| p-value =", round(p_empirical, 3), "\n")
Assoc_results_DegreeTU <- "Estimate = 0.089 | SE = 0.124 | p-value = 0.432"

ResultModel<-summary(model_Deg_TU)
#Reverse log function
Reverselog.Estm<-exp(ResultModel$coefficients$cond[, "Estimate"])
Reverselog.SE<-exp(ResultModel$coefficients$cond[, "Std. Error"])
ResultDeg<-ResultModel$coefficients$cond
ResultDeg<-cbind(ResultDeg, Reverselog.Estm)
ResultDeg<-cbind(ResultDeg, Reverselog.SE)
ResultDegTU<-ResultDeg
ResultDegTU

# Visualisation
hist(perm_estimates, breaks = 50, main = "Permutated estimate Distribution \n for interaction_typeTU")
abline(v = obs_estimate, col = "red", lwd = 2) # For obsered estimate



##Association Partners Model: DEGREE for NTU---
# Full model
model_Deg_NTU <- glmmTMB(
  Degree ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1, 
  family = poisson(link="log"),
  data = rNTU)

hist(residuals(model_Deg_NTU))

VarCorr(model_Deg_NTU) # 1.0417e-05: The random variance associated with egoNTU is virtually zero,
#so there is no reason to keep it in the model (it is not a useful artifact). 
#It does not contribute significantly to capturing inter-individual variance. 
#However, it remains important to control for it in repeated measures.

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
overdisp_fun(model_Deg_NTU)
#Ratio 0.89 No serious overdispersion, NS P-value

##Check for multi-collinearity without interactions
model_Deg_NTU2<-glmmTMB(
  Degree ~ interaction_type + RankStand + Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1, 
  family = poisson(link = "log"), 
  data = rNTU)

check_collinearity(model_Deg_NTU2)

##Full/Null model comparison
mod0<-glmmTMB(Degree ~ 1 + offset(log(grp_size)) + (1 | egoNTU), ziformula = ~1, family = poisson(link = "log"), 
              data = rNTU)
#LRT
lrt<-anova(mod0,model_Deg_NTU,test="Chisq")
lrt #S ok model meaningful

mod1<-glmmTMB(Degree ~ RankStand * Sex + Age.Class
              + offset(log(grp_size)) + (1 | egoNTU), ziformula = ~1, family = poisson(link = "log"), 
              data = rNTU)
#LRT
lrt<-anova(mod1,model_Deg_NTU,test="Chisq")
lrt #S ok test predictor (interaction_type) affect the model

##Check normality of model residuals
hist(residuals(model_Deg_NTU)) 
qqnorm(residuals(model_Deg_NTU))
qqline(residuals(model_Deg_NTU))
Obj<-simulateResiduals(model_Deg_NTU)
plot(Obj) #NS OK!

###Model stability
leave_one_out_effects <- function(model, data) {
  original_coef <- fixef(model)$cond
  n <- nrow(data)
  delta <- matrix(NA, nrow = n, ncol = length(original_coef))
  colnames(delta) <- names(original_coef)
  
  cat("Starting of leave-one-out...\n")
  
  for (i in 1:n) {
    dat_i <- data[-i, ]
    
    # Refit model without observation i
    fit_i <- tryCatch(
      glmmTMB(
        formula = formula(model),
        data = dat_i,
        family = eval(model$call$family),
        ziformula = if (!is.null(model$call$ziformula)) eval(model$call$ziformula) else ~0
      ),
      error = function(e) {
        message(sprintf("❌ EObservation error %d : %s", i, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(fit_i)) {
      coef_i <- fixef(fit_i)$cond
      delta[i, ] <- coef_i - original_coef
    }
    
    if (i %% 10 == 0) cat(sprintf("… %d/%d achived\n", i, n))
  }
  
  cat("Achived.\n")
  return(as.data.frame(delta))
}

delta <- leave_one_out_effects(model = model_Deg_NTU, data = rNTU)
delta_clean <- delta[complete.cases(delta), , drop = FALSE]

# Original Coefficient
coef_original <- fixef(model_Deg_NTU)$cond
# Min et max values calculation after leave-one-out
dfbeta_like <- data.frame(
  orig = coef_original,
  min = coef_original + apply(delta_clean, 2, min, na.rm = TRUE),
  max = coef_original + apply(delta_clean, 2, max, na.rm = TRUE)
)
round(dfbeta_like, 4)

#Recalculation of dfbetas
# 'delta_clean' already : estimate_leave_one_out - estimate_original
dfbeta_equiv <- delta_clean
#Zero centered around predictors
dfbeta_equiv_centered <- scale(dfbeta_equiv, center = TRUE, scale = FALSE)
dfbeta_equiv_centered <- as.data.frame(dfbeta_equiv_centered)

# Add id column
dfbeta_equiv_centered$ObsID <- 1:nrow(dfbeta_equiv_centered)

# Long format transformation
dfbeta_long <- melt(dfbeta_equiv_centered, id.vars = "ObsID",
                    variable.name = "Predictor", value.name = "DFBETA")

# Plot
ggplot(dfbeta_long, aes(x = DFBETA, y = Predictor)) +
  geom_point(alpha = 0.6, color = "skyblue") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", color = "black", fatten = 1.2) +
  theme_minimal() +
  labs(
    title = "Stability plot (DFBETA approach)",
    x = "Centered leave-one-out effect",
    y = ""
  )


# Recall full model
model_Deg_NTU <- glmmTMB(
  Degree ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1,
  family = poisson(link="log"),
  data = rNTU)

# Extraction of estimate for interaction_typeTU
effet_NTU <- summary(model_Deg_NTU)$coefficients$cond["interaction_typeTU", c("Estimate", "Std. Error")]
obs_estimate <- effet_NTU["Estimate"] 
obs_se <- effet_NTU["Std. Error"] 

# Function for data permutation and model adjustment
permute_and_fit <- function(data) {
  permuted_data <- do.call(rbind, lapply(split(data, data$egoNTU), function(subdf) { # Split data by egoTU
    # Subset TU and NTU interaction types
    interaction_types <- subdf$interaction_type # Get unique interaction types for this egoTU
    subdf$interaction_type <- sample(interaction_types)  # permute within ego
    return(subdf) # Permuted subdf
  }))
  
  # Re-do model
  model_perm <- try(glmmTMB(
    Degree ~ interaction_type + RankStand * Sex + Age.Class + offset(log(grp_size)) + 
      (1 | egoNTU),
    ziformula = ~1,
    family = poisson(link = "log"),
    data = permuted_data
  ), silent = TRUE) # Permutation model fitting
  
  if (inherits(model_perm, "try-error")) return(NA) # If model fitting fails, return NA
  
  # Extract the estimate for interaction_typeTU
  est <- try(summary(model_perm)$coefficients$cond["interaction_typeTU", "Estimate"], silent = TRUE) 
  if (inherits(est, "try-error")) return(NA) # If extraction fails, return NA
  
  return(est) # Return the estimate
}

#Permutation
set.seed(123)
n_perm <- 10000
data_perm <- rNTU

#perm_estimates <- replicate(n_perm, permute_and_fit(data_perm)) # Permutation without progress bar
#perm_estimates <- pbreplicate(n_perm, permute_and_fit(data_perm)) #with progression bar
perm_estimates <- pbapply::pbreplicate(n_perm, permute_and_fit(data_perm)) #same but faster

#Empiric p-value
perm_estimates <- perm_estimates[!is.na(perm_estimates)] # Suppress potential errors (NA)

#Bilateral p-value
p_empirical <- mean(abs(perm_estimates) >= abs(obs_estimate))
p_empirical

# Results for interaction_typeTU
cat("Estimate =", round(obs_estimate, 3), "| SE =", round(obs_se, 3), 
    "| p-value =", round(p_empirical, 3), "\n")
Assoc_results_DegreeNTU <- "Estimate = -1.275 | SE = 0.195 | p-value < 0.0001"

ResultModel<-summary(model_Deg_NTU)
#Reverse log function
Reverselog.Estm<-exp(ResultModel$coefficients$cond[, "Estimate"])
Reverselog.SE<-exp(ResultModel$coefficients$cond[, "Std. Error"])
ResultDeg<-ResultModel$coefficients$cond
ResultDeg<-cbind(ResultDeg, Reverselog.Estm)
ResultDeg<-cbind(ResultDeg, Reverselog.SE)
ResultDegNTU<-ResultDeg
ResultDegNTU

# Visualisation
hist(perm_estimates, breaks = 50, main = "Permutated estimate distribution \n for interaction_typeTU")
abline(v = obs_estimate, col = "red", lwd = 2) # Add line for observed estimate



###PLOT MODELS---------------------------------------

#Recall of models
model_Str_TU <- glmmTMB(
  Strength ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) +
    (1 | egoTU),
  ziformula = ~1,
  family = tweedie(link = "log"), 
  data = rTU
)

model_Str_NTU <- glmmTMB(
  Strength ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1,
  family = tweedie(link = "log"), 
  data = rNTU
)

model_Deg_TU <- glmmTMB(
  Degree ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoTU),
  ziformula = ~1, 
  family = poisson(link="log"),
  data = rTU)

model_Deg_NTU <- glmmTMB(
  Degree ~ interaction_type + RankStand * Sex + Age.Class + 
    offset(log(grp_size)) + 
    (1 | egoNTU),
  ziformula = ~1, 
  family = poisson(link="log"),
  data = rNTU)

#STRENGTH
Plot7<-plot_model(model_Str_TU, type = "pred", terms = c("interaction_type"), jitter=T, color="black", dot.size = 5, line.size = 1, show.data = TRUE)
Plot7<-Plot7+ggtitle("Predicted probabilities of \n interaction for a tool user") + theme(plot.title = element_text(size=20))
Plot7<-Plot7+labs(x= 'Interaction type of a \n tool user with', y = 'Grooming duration')+theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
Plot7<-Plot7+theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
Plot7<-Plot7+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
Plot7
Plot8<-plot_model(model_Str_NTU, type = "pred", terms = c("interaction_type"), jitter=T, color="black", dot.size = 5, line.size = 1, show.data = TRUE)
Plot8<-Plot8+ggtitle("Predicted probabilities of \n interaction for a non-tool user") + theme(plot.title = element_text(size=20))
Plot8<-Plot8+labs(x= 'Interaction type of a \n non-tool user with', y = 'Grooming duration')+theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
Plot8<-Plot8+theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
Plot8<-Plot8+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
Plot8
Plot7 + Plot8

#DEGREE
Plot9<-plot_model(model_Deg_TU, type = "pred", terms = c("interaction_type"), jitter=T, color="black", dot.size = 5, line.size = 1, show.data = TRUE)
Plot9<-Plot9+ggtitle("Predicted probabilities of \n interaction for a tool user") + theme(plot.title = element_text(size=20))
Plot9<-Plot9+labs(x= 'Interaction type of a \n tool user with', y = 'Number of partners')+theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
Plot9<-Plot9+theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
Plot9<-Plot9+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
Plot9
Plot10<-plot_model(model_Deg_NTU, type = "pred", terms = c("interaction_type"), jitter=T, color="black", dot.size = 5, line.size = 1, show.data = TRUE)
Plot10<-Plot10+ggtitle("Predicted probabilities of \n interaction for a non-tool user") + theme(plot.title = element_text(size=20))
Plot10<-Plot10+labs(x= 'Interaction type of a \n non-tool user with', y = 'Number of partners')+theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
Plot10<-Plot10+theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
Plot10<-Plot10+theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
Plot10
Plot9 + Plot10
