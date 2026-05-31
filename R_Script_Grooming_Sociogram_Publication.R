##############################R Script Grooming Sociogram#########################################

rm(list=ls())

#####Dataset#####
###Use the Environment_Analysis provided (with 42 individuals)
load("~/Environment_Analysis_Published.RData")

###Library
library(NetExplorer)

###Dataset with information per node
attributes<-data.frame(id=Dataset$id,
                       Sex=ifelse(as.character(Dataset$Sex)=="Male", "M", "F"),
                       Tool.user=ifelse(as.character(Dataset$Tool.user)=="Yes", "Yes", "No"),
                       instrengthG= Dataset$instrengthG,
                       indegreeG= Dataset$indegreeG,
                       outstrengthG= Dataset$outstrengthG,
                       outdegreeG= Dataset$outdegreeG,
                       eigenG= Dataset$eigenG)

#Matrix with interactions
m=mGroom
m <- as.data.frame(mGroom)

#Verifications
colnames(m)[!colnames(m) %in% attributes$id] 
m=m[colnames(m) %in% attributes$id, colnames(m) %in% attributes$id] 
all(colnames(m) %in% attributes$id) 

m<-as.matrix(m)

#Plot sociogram for connectivity during grooming interactions
NetExplorer::vis.net(
  df=attributes, 
  m=m, 
  col.id='id', 
  col.shape = 'Sex',
  shapes = c('triangle', 'circle'), # Shape for 'Male' and 'Female'
  col.size = 'eigenG',
  col.color = 'Tool.user',
  color = c("blue", "#FFA500"), #Color for 'No' and 'Yes' in Tool.user
  layers = 'Tool.user',
  background =  'white')

#Online, you have to draw the sociogram yourself according to what you want to show.
