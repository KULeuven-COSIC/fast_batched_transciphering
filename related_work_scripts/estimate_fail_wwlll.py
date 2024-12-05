# This script estimates the failure probability of the method described in 
# "Fregata: Faster Homomorphic Evaluation of AES via TFHE"
# https://link.springer.com/chapter/10.1007/978-3-031-49187-0_20
#
# The following related works were used as formulas and definitions were not included in the actual work:
#
# https://eprint.iacr.org/2017/430.pdf
# https://eprint.iacr.org/2014/1018.pdf
# https://eprint.iacr.org/2022/198.pdf


import math

nbar = 635
N = 1024
Nbar = 2048
abar = 2**(-15)
a = 2**(-25)
abbar = 2**(-44)

Bg = 2**6
Bgbar = 2**9
Bks = 2**2
Bksbar = 2**3
l=3
lbar = 4
t=7
tbar=10
# https://eprint.iacr.org/2017/430.pdf
beta = Bgbar / 2

# not specified in paper, using value in parameter set from library
k = 1

q_lv0 = 2**32
q_lv1 = 2**32
q_lv2 = 2**64

# https://eprint.iacr.org/2014/1018.pdf
var_ks = q_lv1 * abar / math.sqrt(2 * math.pi)
var_bsk = q_lv2 * abbar / math.sqrt(2 * math.pi)
var_privks = var_bsk

print("Variance of keyswitching key (lv1->lv0): ", var_ks, ". Modulus: ", q_lv0)
print("Variance of bootstrapping key (lv2): ", var_bsk, ". Modulus: ", q_lv2)
print("Variance of private keyswitching (lv2 -> lv1)", var_privks, ". Moduli: ", q_lv2, " (Source) ", q_lv1, " (Target) ")

var_pbs_many = nbar * lbar * (k + 1) * Nbar * (Bgbar**2 + 2) / 12 * var_bsk + nbar * (q_lv2**2 - Bgbar**(2 * lbar))/(24 * Bgbar**(2 * lbar)) * (1 + k * Nbar / 2) + nbar * k * Nbar / 32 + nbar/16 * (1 - k * Nbar / 2)**2

# Note, the authors of Fregata do not cite the proper theorem and the bound for Var_trgsw in the paper for 
# private key switching does not match Theorem 2.16 in the cited paper.
var_trgsw = var_pbs_many + tbar * Nbar * Bksbar**2 * var_privks

# Best case estimate
var_trgsw_best = tbar * Nbar * Bksbar**2 * var_privks

# Var_sbox = 0, they can be given in plaintext
var_tlwe = 8 * ((k + 1) * l * N * beta**2 * var_trgsw) + 0 + nbar * t * N * Bks**2 * var_ks
var_tlwe_best = 8 * ((k + 1) * l * N * beta**2 * var_trgsw_best) + 0 + nbar * t * N * Bks**2 * var_ks

sig_tlwe = math.sqrt(var_tlwe)
sig_tlwe_best = math.sqrt(var_tlwe_best)

plaintext_space = 2

# https://eprint.iacr.org/2022/198.pdf
p_fail = math.erfc(q_lv0 / (2 * plaintext_space * math.sqrt(2) * sig_tlwe))
p_fail_best = math.erfc(q_lv0 / (2 * plaintext_space * math.sqrt(2) * sig_tlwe_best))

print("Failure probability is: ", p_fail)
print("Failure probability in best case is: ", p_fail_best)
