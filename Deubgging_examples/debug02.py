#%%
import numpy as np


#%% 
a = np.array([1, 2, 3])
b = np.array([[1], [2]])
result = a + b
print("result is:", result)

#TODO answer the following questions as comments below:
# What is the shape of a?
# What is the shape of b?
# What is the shape of result?
# Why did this happen? and what occurred? 

#######################################################################
# a is a 1-D array with 3 elements, b is a 2D array with 2 rows and 1 column.
# the result is a 2d array with 2 rows and 3 columns.
# this happened because of broadcasting. numpy tries to match the sizes of the arrays and effective turns a into an array: ([1,2,3],[1,2,3])
# it turns b into: ([1,1,1],[2,2,2]) then adds them.
#######################################################################



#%% 
a = np.array([1, 2, 3])
b = np.array([4])
result = a * b 
print("result is:", result)

#TODO answer the following questions as comments below:
# Is the output what you expected? Can you explain how it happened?
#####################################################################
# The output is what I expected. Numpy uses * to multiply arrays element-wise
######################################################################


#%% 
a = np.ones((2, 3))
b = np.ones((2, 3))
result = a + b  
print("result:", result)

# TODO answer the following questions as comments below:
# Why doesn't this work? If the 2D arrays are meant to be the same size, fix 
# the code so that you can successfully add them. 
######################################################################
# The arrays are not the same size. I've fixed them to both be 2 x 3 and now the output is correct
######################################################################

#%% 

a = np.array([1, 2, 3])
b = np.array([[1], [2], [3]])
result = a + b
print("result:", result)

#TODO - check the output of this addition, is it what you expected? Answer the 
# following questions as comments below:
#  
#  What is the shape of a and the shape of b? 
#  Why is the result 3x3? and what does the addition do in this case?
###########################################################################
# a is of size 1x3. b is 3x1 so when they are added they are broadcast together.
# The output is to be expected but it is somewhat confusing. The result is sa 3x3 because each array was extended
# to match the size of the other
############################################################################



#%%
a = np.array([1, 2, 3])
b = np.array(np.eye(3))
result1 = b*a
result2 = np.dot(b,a)
result3 = b@a

print("result1:", result1)
print("result2", result2)
print("result3", result3)

# which of these is the intended result? 
# TODO put your answer here as a comment and justify it. 
################################################################
# I would say that the third output makes the most sense because it uses matrix multiplication in the traditional sense
################################################################

#%%
a_column = a.reshape((3,1)) 
result1 = b*a_column
result2 = np.dot(b,a_column)
result3 = b@a_column 

print(a_column)
print("result1:", result1)
print("result2", result2)
print("result3", result3)

# which of these is the intended result? 
# TODO put your answer here as a comment and justify it. 
# None of these really make sense. It also depends on what you're looking for. If the intent is to matrix multipy these, that wouldn't normally be possible. 
# In all of these instances, numpy is "extending" (broadcasting) a_column
##########

#%%
result1 = a_column.T*a
result2 = a_column*a 
result3 = a*a_column
result4 = np.dot(a_column,a)
result5 = np.dot(a_column,a_column)
result6 = np.dot(a_column.T, a_column)
result7 = a_column.T@a_column

#TODO: answer these questions as comments below: 
#  Why does np.dot(a, a) work? but np.dot(a_column, a_column) doesn't? 
#  If I'm trying to do find the sum of the squares of the elements of a, which of these should I use?

# np.dot(a_column, a column) is impossible because in numpy they're now 2 dimensional arrays.
# you should use np.dot(a,a)



# if you still have questions about broadcasting, and array dimensions, please see the following:
# https://numpy.org/doc/stable/user/basics.broadcasting.html 
# https://www.youtube.com/watch?v=oG1t3qlzq14 