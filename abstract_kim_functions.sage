r"""
Compute two level 4 abstract Kim functions for X = ℙ¹ ⧵ {0,1,∞} and S a set of primes
with #S = 2. The procedure is described in §5.5 of [CDCL26].
The Kim functions are elements of A(Z)[Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211]
where A(Z) is the (weight-truncated) algebra of unipotent mixed Tate motivic periods 
over Z = Spec(ℤ[1/S]), represented in Sage as the polynomial ring 
    QQ[a1,a2,a12,a112,a122,aσ,a1112,a1122,a1222,a1σ,a2σ].
Here, 1 and 2 are shorthand for the generators τ_{ℓ_1} and τ_{ℓ_2} corresponding to
the two primes in S, and σ is shorthand for σ₃.
The list of the two resulting Kim functions is saved to "abstract_kim_functions.txt" 
in human-readable form and to "abstract_kim_functions.sobj" as a Sage object, so that it
can be loaded via ``load("abstract_kim_functions.sobj")``.

REFERENCES:

    [CDCL26] David Corwin, Ishai Dan-Cohen, Martin Lüdtke, "Mixed Tate Chabauty—Kim beyond 
        the polylog quotient", https://arxiv.org/abs/2610.25365v2

AUTHOR:

    Martin Lüdtke (2026-11-01)

"""


human_readable_output_file = "abstract_kim_functions.txt" # set to None if not wanted
sage_output_file = "abstract_kim_functions.sobj"  # dito

load("shuffle_algebra.sage")

alphabet = ["1","2","σ"]
weights = [1,1,3]
weight_bound = 4

# weight-truncated ring of motivic periods as a polynomial ring in Lyndon variables
A, a = shuffle_algebra_as_polynomial_ring(alphabet, weights, weight_bound)
K = FractionField(A)

selmer_coords = ["x1","y1","x2","y2","z3","z21"]
deRham_coords = ["Log","Li1","Li2","Li3","Li21","Li4","Li31","Li211"]

# (half-)weights of functions Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211
dR_weights = [1,1,2,3,3,4,4,4]

# polynomial rings and variables
R = PolynomialRing(K, names = selmer_coords+deRham_coords)
x1,y1,x2,y2,z3,z21, Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211 = R.gens()
R_dR_mot = PolynomialRing(A, deRham_coords)
R_dR_mot_X.<X> = PolynomialRing(R_dR_mot, 'X')

S = ["1", "2"]  # "1" and "2" represent the two primes in S
x = { "1": x1, "2": x2 }
y = { "1": y1, "2": y2 }

# z3 resp. z21 can be expressed in terms of x1,x2,y1,y2 and Li3 resp. Li21
z3_elim = (Li3 - sum(a(l1+l2+q)*x[l1]*x[l2]*y[q] for l1 in S for l2 in S for q in S))/a("σ")
z21_elim = (Li21 - sum(a(l+q1+q2)*x[l]*y[q1]*y[q2] for l in S for q1 in S for q2 in S))/a("σ")

# x2 can be expressed in terms of x1 and Log
x2_elim = (Log - a("1")*x1)/a("2")

# y2 can be expressed in terms of y1 and Li1
y2_elim = (Li1 - a("1")*y1)/a("2")

# Log*y1 can be expressed in terms of x1 and Log,Li1,Li2
b12 = 1/2*(a("12") - a("21"))
y1_elim =  (Li1*x1 - a("2")/b12*(Li2 - 1/2*Log*Li1))/Log


Li4_pullback = \
    sum(a(l1+l2+l3+q)*x[l1]*x[l2]*x[l3]*y[q] for l1 in S for l2 in S for l3 in S for q in S) \
    + sum (a(l+"σ")*x[l]*z3 for l in S)

Li31_pullback = \
    sum(a(l1+l2+q1+q2)*x[l1]*x[l2]*y[q1]*y[q2] for l1 in S for l2 in S for q1 in S for q2 in S) \
    + sum(a(l + "σ")*x[l]*z21 for l in S) + sum(a("σ" + l)*y[l]*(z3-1) for l in S)

Li211_pullback = \
    sum(a(l+q1+q2+q3)*x[l]*y[q1]*y[q2]*y[q3] for l in S for q1 in S for q2 in S for q3 in S) \
    + sum(a("σ" + l)*y[l]*(z21-1) for l in S)


P4 = Li4 - Li4_pullback
P31 = Li31 - Li31_pullback
P211 = Li211 - Li211_pullback

Ps = [P4, P31, P211]
suffix = ["4", "31", "211"]
print("computing ν4, ν31, ν211...")
nus = []
for i in range(3):
    P = Ps[i]

    # eliminate z3
    P = P.subs(z3 = z3_elim)

    # eliminate z21
    P = P.subs(z21 = z21_elim)

    # eliminate x2
    P = P.subs(x2 = x2_elim)

    # eliminate y2
    P = P.subs(y2 = y2_elim)

    # eliminate y1
    # multiply by power of Log to ensure each y1 can be paired with a Log
    y1_var_index = R.gens().index(y1)
    Log_var_index = R.gens().index(Log)
    needed_exponent = max(0, max(e[y1_var_index] - e[Log_var_index] for e in P.exponents()))
    P *= Log^needed_exponent
    # now substitute y1_elim for y1
    P = P.subs(y1 = y1_elim)
    # we have ensured that the substitution didn't introduce Log's in denominators
    P = R(P)
    
    # clear denominators and divide by gcd of coefficients
    P /= gcd(P.coefficients())

    # view x1 as a new variable X
    nu = R_dR_mot_X(P.subs(x1 = X))
    print(f"computed ν{suffix[i]} of X-degree {nu.degree()}")
    nus.append(nu)
print("")

print("computing abstract Kim functions...")
abstract_kim_functions = []

# now the de Rham variables live in R_dR_mot (without Selmer variables)
Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211 = R_dR_mot.gens()

# check that the functions ν_4, ν_{3,1}, ν_{2,1,1} have x1-degree 2
assert(all(nu.degree() == 2 for nu in nus))

# first Kim function via determinant
M = matrix(R_dR_mot, [nus[j].coefficients() for j in range(3)])
F = M.determinant()
abstract_kim_functions.append(F)

# second Kim function via resultant
F = nus[0].resultant(nus[1])
abstract_kim_functions.append(F)

for i in range(len(abstract_kim_functions)):
    F = abstract_kim_functions[i]

    # normalise coefficients
    F //= gcd(F.coefficients())
    F //= gcd(gcd(r for r in A(c).coefficients()) for c in F.coefficients())

    # both functions are divisible by Li2 - 1/2*Log*Li1, can remove this factor
    fac = Li2 - 1/2*Log*Li1
    assert fac.divides(F)
    F //= fac

    abstract_kim_functions[i] = F
    wt = F.weighted_degree(dR_weights)
    number_of_terms = F.number_of_terms()
    print(f"computed Kim function F{i+1} in weight {wt} with {number_of_terms} terms")
print("")

# we have removed the common factor, the polynomials are now coprime
assert gcd(abstract_kim_functions) == 1

if human_readable_output_file:
    with open(human_readable_output_file, "w") as out:
        for i in range(len(abstract_kim_functions)):
            F = abstract_kim_functions[i]
            print(f"F{i+1} = ", file=out)
            for e,c in F.iterator_exp_coeff():
                print(f"   + {c*parent(F).monomial(e)}", file=out)
            print("", file=out)
    print(f"saved Kim functions in {human_readable_output_file}")

if sage_output_file:
    save(abstract_kim_functions, sage_output_file)
    print(f"saved list of Kim functions as Sage object in {sage_output_file}")
