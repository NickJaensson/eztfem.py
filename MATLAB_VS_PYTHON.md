**Data types**  
All variables stored in a computed program are of a certain *data-type*. 
Integers from mathematics [...,-2,-1,0,1,2,...] are also in computer science  
known as *integers* (or *ints*). The real numbers from mathematics [1.2, 0.973, 
$\pi$, $e$, ...] are known in computer science as *floating point numbers* (or 
*floats*).  
  
Both ints and floats are stored using a certain number of *bits*. For example,
an 8 bit integer can hold 256 values (2^8), or 128 if the sign is stored as
well. We typically need integers of at least 32 bits (slightly more than 4 
billion values). Floats are typically stored as *single* (32 bits) or *double* (64 bits).

|           | Matlab         | Python      |
|-----------| ---------------| ----------- |
| Data types| double, single | float       |
|           | (u)int8, (u)int16, (u)int32, (u)int64 | int |