r"""
Verify Kim's Conjecture for the thrice-punctured line over ℤ[1/3] for all primes p 
up to a given bound. The S-unit equation for S={3} has no solutions. 
If ``use_S3_symm`` is set to ``True``, the S_3-symmetrized polylogarithmic 
level 4 locus `X(ℤₚ)_{{3},PL,4}^{S_3}` is verified to be empty by checking that 
the polylogarithmic locus `X(ℤₚ)_{{3},PL,4}` equals {-1}. Otherwise, the full 
level 4 locus `X(ℤₚ)_{{3},4}` is computed and verified to be empty. 
The S_3-symmetrized variant is more efficient. This proves Theorem 5.9 of [CDCL26].

REFERENCES:

    [CDCL26] David Corwin, Ishai Dan-Cohen, Martin Lüdtke, "Mixed Tate Chabauty—Kim beyond 
        the polylog quotient", https://arxiv.org/abs/2610.25365v2

AUTHOR:

    - Martin Lüdtke (2026-11-01)

"""

load("BeyondPolylog.sage")

prime_bound = 1000
prec = 10
work_parallel = True
use_S3_symm = True

primes = [2] + prime_range(5,prime_bound)
expected_sols = [-1] if use_S3_symm else []

def check_prime(p):
    li_data =  polylog_data(4,p,prec)
    kim_functions = Z_one_third_kim_functions(p,prec,li_data)
    if use_S3_symm:
        # use only the two polylogarithmic functions if use_S3_symm is set
        kim_functions = kim_functions[:2]
    try:
        locus = compute_depth4_locus(kim_functions,li_data)
    except PrecisionError:
        print(f"p = {p}: PrecisionError")
    else:
        if len(locus) == len(expected_sols) or (p == 2 and len(locus) == 0):
            print(f"p = {p}: ✓")
        else:
            print(f"p = {p}: locus contains extra points:")
            for z in locus:
                if z not in expected_sols:
                    print(f"  {z}")

@parallel
def check_primes_parallel(p):
    check_prime(p)

if work_parallel:
    for _ in check_primes_parallel(primes):
        pass
else:
    for p in primes:
        check_prime(p)
