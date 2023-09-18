# Matlab-Python Rosetta Stone
## General
- What is it exactly?
  - *Matlab* is a software product developed by the company Mathworks. Matlab includes an
integrated development enverinment (IDE), the Matlab language, toolboxes, official support, and more.
  - *Python* is an open-source programming language. The Python comminuty drives the 
  development of interpreters/compilers, IDEs, packages, support, and more
- Documentation: 
  - *Matlab*: https://www.mathworks.com/help/matlab/
  - *Python*: https://docs.python.org/3/
- Commenting:
  - *Matlab*: use % for comments.
  - *Python*: uses # for comments.
- In Matlab, a line that does not end with a semicolon (;), will generally write something to the command-line. In Python this is not the case.
- Both Matlab and Python can be used *interactively*
- Indentation in Matlab is optional for readibility, and the keyword ```end```
is used the control flow. In Python, indentation is not optional as it plays a
crucial role in control flow.
- line continutation
- something on parentheses? () vs [] vs {}
- row-major (Python) vs column major (Matlab) (Reshape!)
- Structuring a project:
  - Matlab: functions in same folder or use addpath, 
  - Python: modules in same folder and use import, for different 
  folder not in current folder use sys.path.append (relative vs absolute imports, __init__.py file etc.)
- import math in Python

## Working with variables
### Data-types
All variables stored in a computed program are of a certain *data-type*. 
Integers from mathematics [..., -2, -1, 0, 1, 2, ...] are also in computer science 
known as *integers* (or *ints*). The real numbers from mathematics [1.2, 0.973, 
$\pi$, $e$, ...] are known in computer science as *floating point numbers* (or 
*floats*).  
  
Both ints and floats are stored using a certain number of *bits*. For example,
an 8 bit integer can hold 256 values ($2^8$), or 128 if the sign is stored as
well. We typically need integers of at least 32 bits (slightly more than $4 \times 
10^9$ values). Python *ints* have a variable size, which can grow 
as needed (e.g, as it gets assigned a value larger than it can currently hold).

Floats are typically expressed in base-10 as $A \times 10^B$, where a certain
number of bits is reserved for $A$ (the significand) and for $B$ (the exponent):
- *single precision* occupies 32 bits in total and its significand has a 
  precision of 24 bits (about 7 digits).
- *double precision* occupies 64 bits in total and its significand has a precision 
  of 53 bits (about 16 digits).
  
Matlab does not deal with ints by default, but makes everything floats of double precision.

|            | Matlab           | Python          |
|------------| -----------------| --------------- |
| floats     | double (default) | double (default)|
| ints       | 8 to 64 bits     | grows as needed |
| text       | string and character| string       |
| true/false | logical          | bool            |
| composites | structs and cell-arrays  | lists, dicts and tuples |

### Indexing
<table width="100%">
<tr>
<th style="border-bottom: none;">Matlab</th>
<th style="border-bottom: none;">Python</th>
</tr>
<tr>
<td>

```
a = [4,8,3];
disp(a(1)); % will print the value "4"
```
</td>
<td>

```
a = [4,8,3]
print(a[0]) # will print the value "4"
```
</td>
</tr>
</table>

- slices left and right-inclusive (Matlab) or only left inclusive (Python)
- slicing to end (left and right?)
- reverse?
### Assignment
- Python: variables are labels
- Matlab: variables are boxes
- Assignment statements in Python do not copy objects, they create bindings 
between a target and an object. Use .copy() of .deepcopy() for real copies (Python)
https://docs.python.org/3/library/copy.html
- mutable vs immutable (Python)
- Values copied when assigned to a new variable (Matlab)
- Everything is an array in Matlab, everything is an object in Python

## Control flow

### For-loops
<table width="100%">
<tr>
<th style="border-bottom: none;">Matlab</th>
<th style="border-bottom: none;">Python</th>
</tr>
<tr>
<td>

```
for i = 1:10
    % do stuff (i goes from 1 to 10)
end
```
</td>
<td>

```
for i in range(10):
    # do stuff (i goes from 0 to 9)

```
</td>
</tr>
</table>

### If-else statement

<table width="100%">
<tr>
<th style="border-bottom: none;">Matlab</th>
<th style="border-bottom: none;">Python</th>
</tr>
<tr>
<td>

```
if x < 3
    % do stuff
elseif x < 5
    % do stuff
else
    % do stuff
end
```
</td>
<td>

```
if x < 3:
    # do stuff
elif x < 5:
    # do stuff
else
    # do stuff

```
</td>
</tr>
</table>

### While-loops
<table width="100%">
<tr>
<th style="border-bottom: none;">Matlab</th>
<th style="border-bottom: none;">Python</th>
</tr>
<tr>
<td>

```
while x < 3
    % do stuff
end
```
</td>
<td>

```
while x < 3:
    # do stuff

```
</td>
</tr>
</table>

### Case-matching
<table width="100%">
<tr>
<th style="border-bottom: none;">Matlab</th>
<th style="border-bottom: none;">Python</th>
</tr>
<tr>
<td>

```
switch x
    case 1
        % do stuff
    case 2
        % do stuff
    otherwise
        % do stuff (none of the cases match)
end
```
</td>
<td>

```
match x:
    case 1:
        # do stuff
    case 2:
        # do stuff
    other:
        # do stuff (none of the cases match)
    
```
</td>
</tr>
</table>

### Try and catch exceptions
<table width="100%">
<tr>
<th style="border-bottom: none;">Matlab</th>
<th style="border-bottom: none;">Python</th>
</tr>
<tr>
<td>

```
try
    % try something
catch
    disp("An exception occurred")
end
```
</td>
<td>

```
try:
    print(x)
except:
    print("An exception occurred")

    
```
</td>
</tr>
</table>

## Functions
- Pass by value: Values copied when entering a function (Matlab)
- Pass by assignment: For some datatypes value copied, and for some datatypes
                      reference passed
- Matlab: pass-by-copy https://www.mathworks.com/help/matlab/matlab_prog/avoid-unnecessary-copies-of-data.html
- Python: pass-by-assignment https://docs.python.org/3/faq/programming.html#how-do-i-write-a-function-with-output-parameters-call-by-reference
- neglecting output of function Matlab vs Python
- function in Matlab at the end or in different file (main variables not known
  inside function)

## Object-oriented programming (advanced topic)

## Numerical computing in Python (advanced topic)
