import examples.poisson.poisson1
import examples.poisson.poisson2
import examples.poisson.poisson3
import examples.poisson.poisson4
import examples.stokes.stokes1

result1 = examples.poisson.poisson1.main()
result2 = examples.poisson.poisson2.main()
result3 = examples.poisson.poisson3.main()
result4 = examples.poisson.poisson4.main()
result5 = examples.stokes.stokes1.main()

print('='*20)
print('Test result:')
print(result1==4.071378136849546e-07)
            #  4.071378150172222e-07 (Matlab)
print(result2==0.29483065971774913)
            #  0.294830659717740 (Matlab)
print(result3==0.14934372115313188)
            #  0.149343721153131 (Matlab)
print(result4==6.522070278291991e-07)
            #  6.522070279402215e-07 (Matlab)
print(result5==25.865650243176205)
            #  25.865650243176205 (Matlab)