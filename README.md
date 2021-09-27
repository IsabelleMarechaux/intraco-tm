# intraco-tm
INTRACO theoretical model



Scenarios:

A scenario = a shape of attribute A + a function Performance=f(A) + a definition of Performance + total variance fixed or not

* Attribute:
1. Species differences, no IV: A(i,j,k) = a_i
2. Species differences + Endogenous IV: A(i,j,k) = a_i + c_j
3. Species differences + exogenous IV: A(i,j,k) = a_i +E_k*beta_i
4. Species differences + exogenous IV + endogenous IV: A(i,j,k) = a_i + E_k*beta_i + c_j

The different terms of A could follow the following distributions:
a_i: N(0, var_inter)
beta_i: N(0, var_inter_env) (the beta_i depicts the species-specific responses to environmental variables).
c(k): N(0, var_intra_i)
E: each environmental variables are drawn from N(0, var_site)

The dimension of E (and hence of the beta_i) can vary.

* Performance= f(A):

f(A)=(1-nonlin)*A+ nonlin*A^2

or

f(A)=ax^b

* Performance:

1. ability to preempt a cell
2. fecundity (number of seeds per individual)
3. mortality rate 


* Total variance:
1. fixed to a low value
2. fixed to a high value
3. not fixed (increases when you add a new level of attribute variation).
