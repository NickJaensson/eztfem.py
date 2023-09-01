# eztfem.py

Steps for translating a Matlab loop:
* Copy the Matlab code to the Python file
* Translate to Python *as literal as possible*:
  - let indices run from range(1,n+1)
  - fix notation for loop
  - () to [], but use the same index values
  - remove end at the end of loops
  - remove ;
* Now add -1 to all indices into arrays and vectors A[i] becomes A[i-1]
* The code should run now
* Now replace the loops with the standard Python notation range(n) and add +1 to all loop variables
* DONE!
