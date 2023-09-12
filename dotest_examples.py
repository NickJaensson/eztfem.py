from math import isclose
import examples.poisson.poisson1
import examples.poisson.poisson2
import examples.poisson.poisson3
import examples.poisson.poisson4
import examples.stokes.stokes1

result1 = examples.poisson.poisson1.main()
result2 = examples.poisson.poisson2.main()
result3 = examples.poisson.poisson3.main()
result4 = examples.poisson.poisson4.main()
result5, _, _, _ = examples.stokes.stokes1.main()

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