# This script is used to demonstrate that the rlwe parameters described in TSBS23 are not
# secure against arora-gb.

# Furthermore, even the previous is not true, all parameter sets (Table 1) in TSBS23 are
# vulnerable to the attack against exact FHE schemes [CCP+ 24] (https://eprint.iacr.org/2024/127)
# Due to their high failure rate

from estimator import *

N=2**11
n=2**10

trlwe_alpha=9.6e-11
tlwe_alpha=6.5e-8

# The authors indicate in the work that they use tfhe.github.io for the
# implementation. However, this means that torus itself is discretized over a 32bit
# integer (see polynomial.h in the library)
# hence:
q = 2**32

lwe_params = LWE.Parameters(n=n,q=q,Xs=ND.UniformMod(2),Xe=ND.DiscreteGaussianAlpha(alpha=tlwe_alpha,q=q))
rlwe_params = LWE.Parameters(n=N,q=q,Xs=ND.UniformMod(2),Xe=ND.DiscreteGaussianAlpha(alpha=trlwe_alpha,q=q))

print("Estimating RLWE security")
LWE.estimate(rlwe_params)
print("Estimating LWE security")
LWE.estimate(lwe_params)


