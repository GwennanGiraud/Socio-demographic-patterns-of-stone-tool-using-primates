################R Script to prepare for homophily analysis with python code########################

#This is a R-Script for non-python users.
#You should have Git on your computer to use this code.

library(reticulate)

py_config()

py_install("pandas")
py_install("networkx")
py_install("graph_tools")
py_install("scipy")
py_install("matplotlib")

py_require("pandas")  
py_require("networkx")      
py_require("graph_tools") 
py_require("scipy") 
py_require("matplotlib")

#Get the path to src folder (originally issued from Karimi & Oliveira 2023 paper)
setwd("C:/Users/Gwenn/Desktop/Revisions/GitHub/src")

#Then open the Python_Script_for_homophily_correction.py file into R and run it. 
