##### Ex 3. Some control structures in R #####


##### Ex 3.1 If statements #####


#often in the process of data wrangling or analysis we need need to perform certain actions when a particular criterion(a) is met.
#Like excel's =IF(logical_test,[value_if_true],[value_if_false]), R gives you the capacity to write these statements


if (condition) {
  # do something
} else {
  # do something else
}


#for example: 

num <- 37
if (num > 100) {
  print("greater")
} else {
  print("not greater")
}

# Other useful operators are: & and, | or, == equals, != does not equal, >= greater than or equal, <= less than or equal...  


##### Ex 3.2 For Loops #####

# In many situations we want to repeat a task multiple times. For example, we want to print each word in this sentence

best_practice <- c("Let", "the", "computer", "do", "the", "work")

#it is inefficient to print each subset of the vector 6 times. 
#it's also not scaleable, so wont work  if you later want to use that code on a vector of 7 elements

#it is better to use a for loop, which has the general pattern of 

for(x in y){
  do thing to y x times.
}

#here we can build a function that uses a for-loop to print our words.

print_words <- function(sentence) {
  for (word in sentence) {
    print(word)
  }
}

print_words(best_practice)


#### Q.3.2 Write your own for loop #### 

#write a function called "total" that sums all the elements of a vector. DO NOT USE the inbuilt R function 'sum'.
ex_vec <- c(4, 8, 15, 16, 23, 42)





#some lesson content has been adapted from software carpentry's, "Creating Functions" "Loops in R" & "Making Choices" modules.
#http://swcarpentry.github.io/r-novice-inflammation/15-supp-loops-in-depth/index.html
#http://swcarpentry.github.io/r-novice-inflammation/04-cond/index.html
#http://swcarpentry.github.io/r-novice-inflammation/02-func-R/index.html