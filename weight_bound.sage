"""
Computes an upper bound for the dimension of the level 4 motivic Chabauty—Kim 
ideal I_{Z,4}^{mot} in half-weight ≤ W in the case Z = Spec ℤ[1/S] with #S = 2.

The motivic Chabauty—Kim ideal I_{Z,4}^{mot} is the kernel of the K(Z)-linear
universal cocycle evaluation pullback map
    ev^♯: K(Z)[Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211] → K(Z)[x1,y1,x2,y2,z3,z21]
given by Theorem 4.1 of [CDCL26]. The restriction to half-weight ≤ W is a map of
finite dimensional K(Z)-vector spaces. An upper bound for the dimension of the 
kernel is obtained by arbitrarily specializing the coefficients along a homomorphism
A(Z) → ℚ, and computing the kernel dimension of the corresponding matrix over ℚ.

With the weight bound W = 10, the upper bound for the dimension of 
  I_{Z,4}^{mot} ∩ (hwt ≤ 10)
is computed as 1, which matches the lower bound from Theorem 5.13 and proves 
Theorem 5.14 of [CDCL26].


REFERENCES:

    [CDCL26] David Corwin, Ishai Dan-Cohen, Martin Lüdtke, "Mixed Tate Chabauty—Kim beyond 
        the polylog quotient", https://arxiv.org/abs/2610.25365v2

AUTHOR:

    - Martin Lüdtke (2026-11-01)

"""

# half-weight bound
W = 10

# given a non-negative weight vector w_1,...,w_k and N ≥ 0, computes a list P 
# such that for 0 ≤ n ≤ N, P[n] is a list containing all tuples a = (a_1,...,a_k)
# with weight w(a) := ∑ w_i a_i equal to n
def weighted_partitions(N, weights):
    P = [[] for n in range(N+1)]
    P[0].append((0,)*len(weights)) # in weight 0 we have the partition (0,...,0)
    for i in range(len(weights)):
        w = weights[i]
        for n in range(w,N+1):
            for p in P[n-w]:
                # form new tuple with i-th entry increased by one
                new_p = tuple((p[i]+1 if j == i else p[j]) for j in range(len(weights)))
                P[n].append(new_p)
    return P


load("shuffle_algebra.sage")

alphabet = ["1","2","σ"]
weights = [1,1,3]
max_weight = 4
# weight-truncated ring of motivic periods as a polynomial ring in Lyndon variables
A, a = shuffle_algebra_as_polynomial_ring(alphabet, weights, max_weight)

selmer_coords = ["x1","y1","x2","y2","z3","z21"]
selmer_weights = [1,1,1,1,3,3]

deRham_coords = ["Log","Li1","Li2","Li3","Li21","Li4","Li31","Li211"]
dR_weights = [1,1,2,3,3,4,4,4]

RSel_mot = PolynomialRing(A,selmer_coords)
x1,y1,x2,y2,z3,z21 = RSel_mot.gens()

S = ["1", "2"]  # "1" and "2" represent the two primes in S
x = { "1": x1, "2": x2 }
y = { "1": y1, "2": y2 }

# evaluation map
Log_ev = sum(a(l)*x[l] for l in S)
Li1_ev = sum(a(l)*y[l] for l in S)
Li2_ev = sum(a(l+q)*x[l]*y[q] for l in S for q in S)
Li3_ev = sum(a(l1+l2+q)*x[l1]*x[l2]*y[q] for l1 in S for l2 in S for q in S) + a("σ")*z3
Li21_ev = sum(a(l+q1+q2)*x[l]*y[q1]*y[q2] for l in S for q1 in S for q2 in S) + a("σ")*z21
Li4_ev = sum(a(l1+l2+l3+q)*x[l1]*x[l2]*x[l3]*y[q] for l1 in S for l2 in S for l3 in S for q in S) \
    + sum (a(l+"σ")*x[l]*z3 for l in S)
Li31_ev = sum(a(l1+l2+q1+q2)*x[l1]*x[l2]*y[q1]*y[q2] for l1 in S for l2 in S for q1 in S for q2 in S) \
    + sum(a(l + "σ")*x[l]*z21 for l in S) + sum(a("σ" + l)*y[l]*(z3-1) for l in S)
Li211_ev = sum(a(l+q1+q2+q3)*x[l]*y[q1]*y[q2]*y[q3] for l in S for q1 in S for q2 in S for q3 in S) \
    + sum(a("σ" + l)*y[l]*(z21-1) for l in S)

evaluations = [Log_ev,Li1_ev,Li2_ev,Li3_ev,Li21_ev,Li4_ev,Li31_ev,Li211_ev]

# specialise the periods arbitrarily to some field
K = QQ
periods = 1,2,3,4,5,6,7,8,9,10,11

period_homom = A.hom(periods, K)
evaluations = [ f.change_ring(period_homom) for f in evaluations ]

RdR = PolynomialRing(K, names=deRham_coords)
RSel = PolynomialRing(K, names=selmer_coords)
ev_hom = RdR.hom(evaluations, RSel)


# Sel_bases[n] contains all exponent vectors for Selmer monomials of weight = n
Sel_bases = weighted_partitions(W,selmer_weights)
# similar for dR_bases[n]
dR_bases = weighted_partitions(W,dR_weights)

Sel_basis = [sm for w in range(W+1) for sm in Sel_bases[w]]
dR_basis = [dRm for w in range(W+1) for dRm in dR_bases[w]]

Sel_dim = len(Sel_basis)
dR_dim = len(dR_basis)

# compute matrix for the evaluation map
print(f"Computing {Sel_dim}×{dR_dim} matrix in half-weight ≤ {W}... ", flush=True, end="")
M = matrix(K,Sel_dim,dR_dim)
for j in range(dR_dim):
    ev_im = ev_hom(RdR.monomial(dR_basis[j]))
    for i in range(Sel_dim):
        monom = RSel.monomial(Sel_basis[i])
        M[i,j] = ev_im.monomial_coefficient(monom)
print("✓")

# compute kernel dimension
print("Computing kernel dimension... ", flush=True, end="")
kernel_dim = M.right_nullity()
print("✓")
print(f"Kernel dimension is {kernel_dim}")
