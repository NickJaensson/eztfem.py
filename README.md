# eztfem.py
## Running tests
executing "python -m unittest" runs all the tests in the "test/" folder that
have a filename matching "test*.py" (empty __init__.py in the test directory
is needed for this)


## Building documentation
From the sphinx folder, run the following command 
```
sphinx-build -b html ./source PATH_TO_HTML_BUILDDIR
```

## Steps for translating a Matlab loop:
* Copy the Matlab code to the Python file
* Translate to Python *as literal as possible*:
  - let indices run from range(1,n+1)
  - fix notation for loop
  - () to [], but use the same index values
  - remove end at the end of loops
  - remove ;
* Now add -1 to all indices into arrays and vectors A[i] becomes A[i-1]
* The code should run now
* Now replace the loops with the standard Python notation range(n) and add +1 
to all loop variables
* DONE!

## Checking for unfilled values using zeros
- In EZTFEM (e.g., in mesh_merge.m), sometimes an array filled with zeros is 
created which is supposed to hold node numbers. A check is performed to see if 
the value equals zero to see if a node has been added already. Care should be 
taken when translating this code to Python, since a node can have number zero, 
thus it might errounously think that no node was added! Where possible, it
is best to initialize arrays like this to -1 (and the modify the checks 
acoordingly)

## Neglecting output of a function
When assigning the output of a function to a single variable, in Python that
variable will be a tuple containing all output of the function:

def func():
    return a, b

aa = func() # aa will be a tuple containing both a and b
aa, _ = func() # only get the argument a from the function


In Matlab is only takes the value of the first variable being returned:

function [a,b] = func()

aa = func % aa will only contain the value of a
[aa,~] = func % same as the line above (but clearer)

