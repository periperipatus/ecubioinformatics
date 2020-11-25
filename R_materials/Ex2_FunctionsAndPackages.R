#### Ex 2: Functions and Packages ####
#### Learning outcomes #####
# - Clearing R environment using 'rm(list=ls())'
# - Structure of functions.
# - installing packages
# - loading packages to the environment 
# - loading package functions
# - writing your own function


#We no longer need data from the previous exercise
rm(list=ls())

#load some data for the exercise.
dune<- read.csv("data/dune.csv")
str(dune) #a data frame of observations of 30 species at 20 sites in Dutch dune meadows.

#RESEARCH QUESTION: What is the species diversity at each site in these dunes?
#We could write our own function to calculate diversity.
#in R, writing your own function has the general format of: 
function_name <- function(argument1, argument2,...){                 
  statement(argument1)
  statement(argument2)
  return(output)                 
} 

### a function to calculate species diversity might look something like this

simpson_diversity<- function(data){
total<- sum(data)
D<- vector(mode="numeric")
for(i in 1:length(data)){
  species_count<- data[i]
  freq<- species_count/total
  D[i]<- freq^2
}
D<- na.omit(D)
return(1-sum(D))
}

simpson_diversity(dune[1,])
apply(dune, 1, simpson_diversity) #remember 'apply' repeats the application of a function over rows or columns.


### HOWEVER, if you then want to calculate shannon diversity and other common metrics, why re-invent the wheel?
### chances are, someone has already written a package with the function you want! 
### in this case, vegan is a package that contains the function(s)

install.packages("vegan") #installs a package of interest. Here we are using vegan which has a lot of functions designed for community ecology

library(vegan) # loads package into R environment.

diversity(dune,index = "simpson")

?diversity
#can also calculate shannon diversity
diversity(dune, index="shannon")


#huge time saver!! 

#package documentation includes all the different functions packages have: https://cran.r-project.org/web/packages/vegan/vegan.pdf
#some packages have many functions - others have one or two.


#### For the next two exercises you will need to install the following packages

install.packages("knitr")
install.packages("plyr")
install.packages("ggplot2")
install.packages("lubridate")
install.packages("oz")
install.packageS("sp")


# Here we will learn more about writing functions. 

# This example function that converts Fahrenheit to Kelvin. # Each line of code is annotated with comments so you can understand how it works

fahrenheit_to_kelvin <- function(temp_F) #fahrenheit_to_kelvin object is defined as a function with 1 argument(temp_f).
{ #the curly bracket indicates the body of the function which executes the function statements
  temp_K <- ((temp_F - 32) * (5 / 9)) + 273.15 #Take the value given as temp_F, and perform the arithmetic and deliver to new object
  return(temp_K) #show the answer
}

fahrenheit_to_kelvin(100) # apply your function.


hot_day_K<- fahrenheit_to_kelvin(100) #Function outputs can be saved as new objects

## EXERCISE: write a function called 'outside' that returns the first and last elements of a vector.
#Use the following vector to test your function.
vector_characters<- c("I'm", "hungry","I", "want","to", "eat","something","cheesy") #c function creates a vector composed of n elements.

#HINT: Function should be used like this: outside(vector_characters)
#HINT: if variable v is a vector then, v[1] is the first element and v[length(v)] is the last element.



#### SOURCE ####

#Dune data came as part of the vegan package.
#some lesson content has been adapted from software carpentry's "Creating Functions" module.
#http://swcarpentry.github.io/r-novice-inflammation/02-func-R/index.html