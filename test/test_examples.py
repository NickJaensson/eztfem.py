from math import isclose
import sys
sys.path.append('..')
import examples.poisson.poisson1
import examples.poisson.poisson2
import examples.poisson.poisson3
import examples.poisson.poisson4
import examples.stokes.stokes1
import examples.stokes.stokes2
import examples.stokes.streamfunction1

# tests for poisson
result1 = examples.poisson.poisson1.main()
result2 = examples.poisson.poisson2.main()
result3 = examples.poisson.poisson3.main()
result4 = examples.poisson.poisson4.main()

# tests for stokes
result5, _, _, _ = examples.stokes.stokes1.main()
result6, _, _, _ = examples.stokes.stokes2.main()

# test for streamfunction
_, mesh, problem, u = examples.stokes.stokes1.main()
result7 = examples.stokes.streamfunction1.main(mesh,problem,u)

eps = 1e-15

print('='*20)
print('Test result:')
print(isclose(result1, 4.071378136849546e-07, abs_tol=eps))
                       # 4.071378150172222e-07 (Matlab)
print(isclose(result2, 0.29483065971774913, abs_tol=eps))
                     # 0.294830659717740 (Matlab)
print(isclose(result3, 0.14934372115313188, abs_tol=eps))
                     # 0.149343721153131 (Matlab)
print(isclose(result4, 6.522070278291991e-07, abs_tol=eps))
                     # 6.522070279402215e-07 (Matlab)
print(isclose(result5, 25.865650243176205, abs_tol=eps))
                     # 25.865650243176205 (Matlab)
print(isclose(result6, 0.0833333333333339, abs_tol=eps))
                     # 0.083333333333333 (Matlab)    
print(isclose(result7, 0.008522786557203416, abs_tol=eps))
                     # 0.008522786557203416 (Matlab)                   