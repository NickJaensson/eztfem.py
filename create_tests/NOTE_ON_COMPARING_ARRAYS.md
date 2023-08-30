**Test**  
Arrays (and scalars) in Matlab have a minimal number of dimensions of 2, 
which can be shown using the following code:
```
sca=2; disp(size(sca))            % prints: 1  1
vec=[0,0,0]'; disp(size(vec))     % prints: 3  1
```
In contrast, in Python, we do have *true* 1D arrays:
```
sca=2; print(sca.shape)                 # prints an error
vec=np.array([0,0,0]); print(vec.shape) # prints:   (3,)
```

A problem arrises when comparing arrays from Matlab to arrays in Python. The
reason for this is that without additional information, we don't know if a 
N x 1 Matlab array is intended as *true* 1D, or just happens to have one 
dimenions equal to 1. Since Python does make this distinction (an array 
of size N cannot be equal to an array of size N x 1), we might get
false negatives when comparing.  
  
A number of steps are taken in this code to avoid this problem.
  
*On the Matlab side:*
- write scalars as scalars (using isscalar)
- write 2D arrays (so can be 1D or "true" 2D) as 1D when one of the 
  dims == 1
  
*And on the Python side:*
- use numpy.squeeze on all 2D arrays to get rid of dim == 1

Now it does not matter if the arrays are of size N or N x 1, since they end
results will be the same on both the Matlab side and the Python side.  
  
NOTE: this problem also arises with matrices of dimension 3 (or higher), as the
following code demonstrates:
```
mat=zeros(1,3,3); disp(size(mat))            % prints: 1  3  3
mat=zeros(3,1,3); disp(size(mat))            % prints: 3  1  3
mat=zeros(3,3,1); disp(size(mat))            % prints: 3  3
```
A solution for these cases is to always use squeeze on the Matlab size, as 
well as squeeze on the Python side.