### Exercise 1 - R & R Studio Orientation ###

#### Learning outcomes #####
#### Ex 1.1 
# - Navigating R-Studio - File browser, Global environment
# - Variable asignment
# - Arithmetic
# - Numerical and character vectors
# - Vector indexing & subsetting
#### Ex 1.2 
# - specifying multiple values
# - viewing object contents
# - number of vector elements using 'length'
# - numeric and character data
# - interrogating objects
#### Ex 1.3
# - addressing vector elements
# - including and excluding vector elements.
# - replacing data.
#### Ex 1.4
# - Setting working directory.
# - Use commandline and 'base' code to import data "read.csv"
# - Use R-studio GUI "import" function with tibble package "read_csv"
# - data frames
# - Interrogate data frames using 'str'
#### Ex 1.5 subsetting arrays
# - using square brackets and indexing to access rows and columns
# - Addressing data.frame variables by column number or name
# - Subset data using logical operators.
#### Ex 1.6 Calculating Statistics
# - use of standard summary statistics functions, 'min' 'max' 'mean' 'sd' 'median'
# - peculiarities of input data.
#### Ex 1.7 Repeating operations
# - Operation automation through 'apply' function
# - Function arguments.
# - Accessing R documentation through R-studio.
#### Ex 1.8 Plotting Data
# - Plot scatter data using base 'plot'
# - Plot histogram using base 'hist'
# - Simple plot modifications
# - Plot export

##### Ex 1.1 - Layout of R-studio & Variable assignment #####



#The text editor is where this text is displayed, and can be used to direct code into the R commandline (below)
#Text preceeded by a '#' is a comment, and will not be read by R as executable code. 
#Use commenting to explain to your future self and other users what particular code does!

x<- 4 #create an object for the number four. To run this code press ALT+ENTER or have the caret(text cursor) in the relevant line and click 'Run'.
y<- 8 #code that has been run appears below in the R interface.

x+y #perform arithmetic with these R objects. The output of these appears in the R console below.
x*y
x/y
x^y
4/(x+y)

#can also perform arithmetic without assigning a variable first. e.g.
13 + 121



#### Ex 1.2 - Creating Vectors ####

a<- 1:10  #make an object called a with values 1 through 10.
a
b<- c(1,6,9) #make an object with non-sequential numbers. "c" is for combine.
b
d<- seq(from=0, to=100, by=10)
length(d) #how many elements in this vector

e<- c("I", "can", "not", "wait", "for", "tea") #create a character vector
class(a) #function class interrogates object to tell you how R is storing the information.
class(e) #class information can also be seen in the environment pane.

#### Ex 1.3 - Subsetting Vectors ####

e[-3] ## square brackets indicate that you want to extract something from within the object. This indicates what you want to exclude.
e[c(1,2,4:6)] #this indicates which vector elements you want to INCLUDE


#replace vector elements
e[2:3]<- c("have", "to") 
e



##### Ex 1.4 - Loading external data #####
# Set the working directory in R-Studio so R knows where to find your data. Can be done via R-Studio. But shows code in R Console.


#import via commandline
infl<- read.csv("./data/inflammation-01.csv", header=FALSE) 

#these data reflect measures of inflammation over days (columns) and each row represents observations from a single patient


#data can also be imported using the R-Studio GUI. Let's work together to import another file "aneurism.csv". Let's call the object 'anuerism'
#these data reflect case control studies with some aneurism measure as the response variable.

# To view the type of object it is run:
class(infl)


#Within dataframes, different columns can be stored with different variable types. 
str(dat) #gives you a summary of the data, showing the names of the variables (V1-40 here), and the type of variable (integer)
str(aneurism) 

# Data frames are essentially stored as a series of aggregated vectors (columns) with different properties (e.g. Numeric, Integer, Character, Factor, Boolean)
# different column variables can be accessed using the $ operator.
class(aneurism$Group)



##### Ex 1.5 Subsetting arrays of data ####

#data in the dataframe can be subsetted using the index provided in square brackets. 
infl[1,1] # The first number specifies the row, and the second (after the comma) the column.

##NOTE: Unlike Python, R is indexed from 1. 

infl[c(1,4),] #can select row 1 & 4 (patient 1 & 4)
infl[,10:12] # or multiple columns (days 10 through 12)

#columns can also be addressed by name.
infl[,c("V10","V11","V12")]


    ##### Q. 1.5.1 - Create a new dataframe object that only includes patient ID, Blood pressure and Age. ####
    names(aneurism) # use this to find out your column names

    
#Can also use logical operators to subset data
#e.g we want inflammation observations that are greater than 3 on day 5.
infl[infl$V5>3,]

aneurism[aneurism$Group=="Control",]


    ##### Q 1.5.2 - Create a new dataframe object from the aneurism data for all individuals greater than or equal to 16 years old#####
    # HINT: use <= or >= 


##### Ex 1.6 - Calculating statistics #####

#what is the maximum inflamation patient two experienced?
max(infl[2,])

# Across all patients, what's the minimum inflammation on day 7
min(infl[,7])

#What's the mean inlammation on day 7
mean(infl[,7])

#what's the mean inflammation for patient 7?
mean(infl[7,])

#these data return a warning message and an empty answer because the function mean does not convert in the dataframe rows
#into the vector it requires. Need to use an additional function to convert the rows into a numeric vector of data.
#see 
class(infl[7,])
class(infl[,7])

mean(as.numeric(infl[7,]))


    ##### Q 1.6.1 - calculate mean, median and standard deviation for treatment group 1 at variable aneurism q1 ####
    #HINT: Functions are median, sd.  


##### Ex 1.7 - Repeating operations #####

# for many analyses the same operation will need to be completed repeatedly on different subsets of the data.
# if you have done programming before, you will use for loops to perform this function.
# R has developed simple wrapper functions that can be used in this way.


#the function "apply" allows us to repeat a function (e.g. a mean) on all the rows, or columns of a data frame.
#R has an inbuilt manual for all functions, that tell you what arguments you need.

?apply

#apply requires 3 arguments, x= data, margin= column or row, fun=function to apply.


#to calculate the mean inflamation per patient

mean_per_patient<- apply(infl, 1, mean)
mean_per_day<- apply(infl,2,mean)

  ##### Q 1.7.1 - Calculate the standard deviation of inflamation per day ####
  
  

##there are many apply family functions  
#apply - apply over the margins of an array (e.g. the rows or columns of a matrix)
#lapply - apply over an object and return list
#sapply - apply over an object and return a simplified object (an array) if possible
#vapply - similar to sapply but you specify the type of object returned by the iterations

#these functions can be used to read in multiple datasets that have a similar naming structure to a list object

filenames <- list.files(path = "data", pattern = "inflammation-[0-9]{2}.csv", full.names = TRUE)
data<- lapply(filenames, FUN=read.csv, header=FALSE)


##### Q 1.7.2 repeat the above question using the second item in the 'data' list #####
#HINT: use [[]] to subset from a list.









##### Ex 1.8 Plotting the data #####

plot(mean_per_day) #basic plot, but could be easier to read.

plot(mean_per_day, xlab="Day", ylab="Mean Inlammation across subjects", pch=16, col="blue")

hist(aneurism$BloodPressure, breaks=10,ylim=c(0,20), xlab="Systolic blood pressure", main="Histogram of Blood Pressure in Aneurism study")

#### SOURCE ####

#the data and much of the lesson content was derived from Software Carpentry's R programming unit
# http://swcarpentry.github.io/r-novice-inflammation/01-starting-with-data/index.html

