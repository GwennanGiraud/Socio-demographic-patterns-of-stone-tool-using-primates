#####################################Script R SNA Metrics calculation###############################

rm(list=ls())

#####Import libraries#####
library(ANTs)


##########Datasets########################
###Grooming dataset
#Data has been primarily cleaned from non-focal individuals and infants 
GData<- as.data.frame(read.csv("Grooming_Dataset.csv", header=TRUE, sep=",", dec=".")) 

###Attribute dataset
attributs <- as.data.frame(read.csv("Attributes_Table.csv", header=TRUE, sep=",", dec="."))



#####Transform the dataframe into matrix########
###Create an object to take into account the sampling effort for each individual
Tobs<-attributs$Focal.Effort

###Matrix construction
mGroom= df.to.mat(GData, actor = "Groomer", receiver = "Receiver", weighted = "Duration", sym = FALSE, tobs = Tobs, num.ids = FALSE) 



#############SNA metrics calculation######################
###Create an empty data frame the same size as the matrix "mGroom"
data.groom = df.create(mGroom)

###Metrics calculation for grooming network
data.groom = instrength = met.instrength(mGroom, df = data.groom, dfid=1) 
data.groom = outstrength = met.outstrength(mGroom, df = data.groom, dfid=1)
data.groom = indegree = met.indegree(mGroom, df = data.groom, dfid=1)
data.groom = outdegree = met.outdegree(mGroom, df = data.groom, dfid=1)
data.groom = eigen = met.eigen(mGroom, df = data.groom, dfid=1)
names(data.groom)[2:9] <- c("instrengthG","outstrengthG","indegreeG","outdegreeG", "eigenG")

###Combine the data
data.all = merge(attributs, data.groom, by = "id", all.x = T) # get information together
Data.all <- subset(data.all, select=c(id, name, Age.Class, Age.Group, Sex, Rank, Hair.pattern, Tool.user, n_Focals,
                                          Focal.Effort, instrengthG, outstrengthG, indegreeG, outdegreeG, eigenG))
