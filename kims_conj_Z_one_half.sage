r"""
Verify Kim's Conjecture for the thrice-punctured line over ℤ[1/2] for all primes p 
up to a given bound. The code computes the Chabauty—Kim locus in X(ℤₚ) for the 
full level 4 quotient of the fundamental group and verifies that it contains
only the 3 solutions {2, -1, 1/2} of the S-unit equation for S={2}. This proves
Theorem 5.7 of [CDCL26].

REFERENCES:

    [CDCL26] David Corwin, Ishai Dan-Cohen, Martin Lüdtke, "Mixed Tate Chabauty—Kim beyond 
        the polylog quotient", https://arxiv.org/abs/2610.25365v2

AUTHOR:

    - Martin Lüdtke (2026-11-01)

"""

load("BeyondPolylog.sage")

prime_bound = 50
prec = 10
work_parallel = True

s_unit_sols = [ 2, -1, 1/2 ]
primes = prime_range(3,prime_bound)

def check_prime(p):
    li_data =  polylog_data(4,p,prec)
    #kim_functions = Z_one_half_kim_functions(p,prec,li_data)
    # using only the first two of the four Kim functions is enough:
    kim_functions = Z_one_half_kim_functions(p,prec,li_data)[:2]
    try:
        locus = compute_depth4_locus(kim_functions,li_data)
    except PrecisionError:
        print(f"p = {p}: PrecisionError")
    else:
        if len(locus) < len(s_unit_sols):
            print(f"p = {p}: something is wrong, locus contains only {len(locus)} points!")
        elif len(locus) == len(s_unit_sols):
            print(f"p = {p}: ✓")
        else:
            print(f"p = {p}: locus contains extra points:")
            for z in locus:
                if z not in s_unit_sols:
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
