r"""
Code for computing level 4 Chabauty—Kim loci beyond the polylog quotient for 
the thrice-punctured line X = ℙ¹ ⧵ {0,1,∞}, as described in [CDCL26].

EXAMPLE:

    We compute the level 4 Chabauty—Kim locus for Z = Spec(ℤ[1/6]) using p = 37,
    and verify that Kim's Conjecture holds, i.e., the locus contains precisely
    the 21 solutions to the S-unit equation for S = {2,3}.

    sage: p = 37; prec = 10
    sage: abstract_kim_functions = load("abstract_kim_functions.sobj")
    sage: periods = Z_one_sixth_periods(p,prec)
    sage: p_adic_kim_functions = [ specialize_kim_function(f,periods) for f in abstract_kim_functions]
    sage: L = compute_depth4_locus(p_adic_kim_functions); L
    [2 + O(37^7),
     3 + O(37^7),
     4 + O(37^7),
     5 + 4*37 + 4*37^2 + 4*37^3 + 4*37^4 + 4*37^5 + 4*37^6 + O(37^7),
     9 + O(37^7),
     10 + 9*37 + 9*37^2 + 9*37^3 + 9*37^4 + 9*37^5 + 9*37^6 + O(37^7),
     12 + 12*37 + 12*37^2 + 12*37^3 + 12*37^4 + 12*37^5 + 12*37^6 + O(37^7),
     13 + 12*37 + 12*37^2 + 12*37^3 + 12*37^4 + 12*37^5 + 12*37^6 + O(37^7),
     15 + 32*37 + 13*37^2 + 32*37^3 + 13*37^4 + 32*37^5 + 13*37^6 + O(37^7),
     18 + 18*37 + 18*37^2 + 18*37^3 + 18*37^4 + 18*37^5 + 18*37^6 + O(37^7),
     19 + 18*37 + 18*37^2 + 18*37^3 + 18*37^4 + 18*37^5 + 18*37^6 + O(37^7),
     20 + 18*37 + 18*37^2 + 18*37^3 + 18*37^4 + 18*37^5 + 18*37^6 + O(37^7),
     23 + 4*37 + 23*37^2 + 4*37^3 + 23*37^4 + 4*37^5 + 23*37^6 + O(37^7),
     25 + 24*37 + 24*37^2 + 24*37^3 + 24*37^4 + 24*37^5 + 24*37^6 + O(37^7),
     26 + 24*37 + 24*37^2 + 24*37^3 + 24*37^4 + 24*37^5 + 24*37^6 + O(37^7),
     28 + 27*37 + 27*37^2 + 27*37^3 + 27*37^4 + 27*37^5 + 27*37^6 + O(37^7),
     29 + 36*37 + 36*37^2 + 36*37^3 + 36*37^4 + 36*37^5 + 36*37^6 + O(37^7),
     33 + 32*37 + 32*37^2 + 32*37^3 + 32*37^4 + 32*37^5 + 32*37^6 + O(37^7),
     34 + 36*37 + 36*37^2 + 36*37^3 + 36*37^4 + 36*37^5 + 36*37^6 + O(37^7),
     35 + 36*37 + 36*37^2 + 36*37^3 + 36*37^4 + 36*37^5 + 36*37^6 + O(37^7),
     36 + 36*37 + 36*37^2 + 36*37^3 + 36*37^4 + 36*37^5 + 36*37^6 + O(37^7)]
    sage: s_unit_sols = [
        2, -1, 1/2, 
        3, 1/3, 2/3, 3/2, -1/2, -2, 
        4, 1/4, 3/4, 4/3, -1/3, -3, 
        9, 1/9, 8/9, 9/8, -1/8, -8
    ]
    sage: len(L) == len(s_unit_sols)
    True
    sage: all(z in L for z in s_unit_sols)
    True


REFERENCES:

    [CDCL26] David Corwin, Ishai Dan-Cohen, Martin Lüdtke, "Mixed Tate Chabauty—Kim beyond 
        the polylog quotient", https://arxiv.org/abs/2610.25365v2

AUTHOR:

    Martin Lüdtke (2026-11-01)

"""


# use faster implementation of p-adic polylogarithms from https://github.com/martinluedtke/RefinedCK
load("RefinedCK.sage")
load("Zproots.sage")

load("shuffle_algebra.sage")


def compute_depth4_locus(kim_functions, li_data=None):
    r"""
    Compute the common zero locus in X(ℤₚ) of a list of polynomials 
    in ℚₚ[Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211], where X = ℙ¹ ⧵ {0,1,∞}.
    First the zero locus of ``kim_functions[0]`` is computed using ``Zproots``
    method, then the vanishing of the remaining Kim functions is checked.

    INPUT:

        - ``kim_functions`` --  a polynomial or a list of polynomials in ℚₚ[Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211]

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.

    OUTPUT:

        - the list of common roots in X(ℤₚ)
    
    EXAMPLES:

        Verify that the common vanishing locus in X(ℤ_{11}) of the functions Log and Li2 is {-1}:

        sage: K = Qp(11, prec=10)
        sage: R.<Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211> = PolynomialRing(K)
        sage: L = compute_depth4_locus([Log, Li2]); L
        [10 + 10*11 + 10*11^2 + 10*11^3 + 10*11^4 + 10*11^5 + 10*11^6 + 10*11^7 + 10*11^8 + 10*11^9 + O(11^10)]
        sage: L == [-1]
        True

        Verify Kim's Conjecture for S = {2} and p = 13:

        sage: p = 13; prec = 10
        sage: kim_functions = Z_one_half_kim_functions(p,prec)
        sage: L = compute_depth4_locus(kim_functions); L
        [2 + O(13^9),
        7 + 6*13 + 6*13^2 + 6*13^3 + 6*13^4 + 6*13^5 + 6*13^6 + 6*13^7 + 6*13^8 + O(13^9),
        12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + 12*13^7 + 12*13^8 + O(13^9)]
        sage: L == [2,1/2,-1]
        True

    """

    # if a single Kim function is provided, make it a singleton list
    if not isinstance(kim_functions, list):
        kim_functions = [kim_functions]
    
    if len(kim_functions) == 0:
        raise ValueError("Need at least one Kim function!")
    
    if kim_functions[0] == 0:
        raise ValueError("Kim function need to be non-zero!")

    K = kim_functions[0].base_ring()
    p = K.prime()
    prec = K.precision_cap()
    Kt.<t> = PolynomialRing(K)

    # compute polylog data if not provided
    if li_data is None:
        li_data = polylog_data(4,p,prec)

    solutions = []

    # loop over residue discs
    for residue in range(2,p):
        # root of unity in this residue disc
        ζ = K.teichmuller(residue)

        # power series expansions for the multiple polylog functions
        mpolylog_series = compute_level4_multiple_polylog_series(residue, p, prec, li_data)

        # compute in the Tate algebra to automatically keep track of valuations
        T = TateAlgebra(K, names='t')
        mpolylog_series = [ T(series).add_bigoh(prec) for series in mpolylog_series ]

        # compute power series expansions of Kim functions
        kim_series = []
        precisions = []
        for F in kim_functions:
            f = F(mpolylog_series)  # substitute power series for the multiple polylog functions
            f //= p^(f.valuation()) # normalise
            # lift Tate algebra element to polynomial and remember its p-adic precision
            f_poly = sum(term.coefficient() * t^term.exponent()[0] for term in f.terms())
            f_prec = f.precision_absolute()
            if f_prec <= 0:
                raise PrecisionError(f"Kim function has coefficients of precision <= 0")
            kim_series.append(f_poly)
            precisions.append(f_prec)

        # can skip residue disk if functions are coprime mod p
        fmodp = [f.change_ring(GF(p)) for f in kim_series]
        if any(gcd(fmodp[0], fmodp[i]) == 1 for i in range(1,len(fmodp))):
            continue

        # check which roots of the first Kim function are also roots of the other
        roots = Zproots(kim_series[0], precisions[0])
        for root in roots:
            if all(eval_poly(kim_series[i],root).add_bigoh(precisions[i]) == 0 for i in range(1,len(kim_series))):
                solutions.append(ζ + p*root)

    return solutions


def specialize_kim_function(kim_function, periods):
    r"""
    p-adic specialization of Kim functions. Given a polynomial over a polynomial ring Q[x_1,...,x_n]
    and a list of p-adic values a_1,...,a_n, substitutes those values for the variables x_1,...,x_n.

    """

    A = kim_function.base_ring()
    K = periods[0].parent()
    period_homom = A.hom(periods, K)
    F = kim_function.change_ring(period_homom)
    return F


def multiple_polylog_functions(p,prec,li_data=None):
    r"""
    Given a prime p and precision prec, returns a list of functions
        logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211
    which compute p-adic multiple polylogarithms up to level 4. 
    This uses the ``fast_polylog`` function from "RefinedCK.sage" and 
    exploits the fact that all multiple polylogarithms up to level 4 
    can be expressed in terms of single polylogarithms.

    INPUT:

        - ``p`` - a prime

        - ``prec`` - a positive integer specifying the precision

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.
    
    OUTPUT: a list of functions ``[logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211]``
        computing p-adic multiple polylogarithms
    

    EXAMPLE:

        sage: p = 7; prec = 10
        sage: logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211 = multiple_polylog_functions(p,prec)
        sage: Lip31(100)  # Li_{3,1}(100)
        3*7^4 + 5*7^5 + 5*7^6 + 4*7^7 + 3*7^8 + 4*7^9 + O(7^10)

    """

    # compute polylog data if not provided
    if li_data is None:
        li_data = polylog_data(4,p,prec)
    g,k0 = li_data

    K = Qp(p,prec)
    if p != 3:
        ζ3 = -4/3*fast_polylog(-1,3,p,prec,g,k0)  # ζ(3) = -4/3 Li_3(-1)
    else:
        # use alternative formula for ζ(3) when p = 3 to avoid precision loss
        ζ3 = 8/7*fast_polylog(2,3,p,prec,g,k0)  # ζ(3) = 8/7 Li_3(2)

    # exploit the fact that multiple polylogarithms up to level 4 can be expressed 
    # in terms of single polylogarithms
    logp = lambda x: K(x).log()
    Lip1 = lambda x: -K(1-x).log()
    Lip2 = lambda x: fast_polylog(x,2,p,prec,g,k0)
    Lip3 = lambda x: fast_polylog(x,3,p,prec,g,k0)
    Lip21 = lambda x: ζ3 + 1/2*logp(x)*Lip1(x)^2 - Lip1(x)*Lip2(1-x) - Lip3(1-x)
    Lip4 = lambda x: fast_polylog(x,4,p,prec,g,k0)
    Lip31 = lambda x: Lip4(x)  - Lip4(1-x) + Lip4(x/(x-1)) + Lip1(x)*Lip3(x) \
            + 1/24*Lip1(x)^4 + 1/6*logp(x)*Lip1(x)^3 - ζ3*Lip1(x)
    Lip211 = lambda x: 1/6*logp(x)*Lip1(x)^3 - Lip4(1-x) - Lip1(x)*Lip3(1-x) \
            - 1/2*Lip1(x)^2*Lip2(1-x)
    
    return logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211


def compute_level4_multiple_polylog_series(residue, p, prec, li_data=None):
    r"""
    Compute the power series expansions on a given residue disc for the multiple
    polylogarithm functions logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211.
    The residue disc is parametrized as ζ + pt with t ∈ ℤₚ and ζ the unique
    root of unity in the disc. The power series are returned as polynomials in ``t``
    which are a correct approximation of the power series modulo p^{prec}.

    INPUT:
    
        - ``residue`` -- a p-adic integer specifying the residue disc. It has to be
        different from 0 or 1 modulo p.

        - ``p`` -- a prime

        - ``prec`` -- a positive integer specifying the precision

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.

    """

    # compute polylog data if not provided
    if li_data is None:
        li_data = polylog_data(4,p,prec)
    g,k0 = li_data

    # delta = floor(log_p(k0-1))
    delta = 1
    while p^(delta+1) < k0:
        delta += 1
    M = max(prec, prec + (delta-1))
    K = Qp(p,M)  # store M p-adic digits for each coefficient
    R = PolynomialRing(K, names='t')
    ζ = K.teichmuller(residue)

    log_series = R([0] + [-(-p)^k/(k*ζ^k) for k in range(1,k0)])
    Li_series = compute_polylog_series(residue, 4, p, M, g, k0)
    logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211 = multiple_polylog_functions(p,M,li_data)

    # given series expansion for Li_w, compute the one for Li_{e₀w} but without constant coefficient
    def e0_integrate(series): 
        coeffs = [0]
        for k in range(1,k0):
            coeffs.append(-1/k*sum((-p/ζ)^(k-j)*series[j] for j in range(k)))
        return R(coeffs)
    
    # get Li_{2,1} = Li_{011} from Li_{11} = 1/2 Li_1^2
    Li11_series = 1/2 * Li_series[1]^2
    Li21_series = Lip21(ζ) + e0_integrate(Li11_series)

    # get Li_{3,1} = Li_{0011} from Li_{2,1} = Li_{011}
    Li31_series = Lip31(ζ) + e0_integrate(Li21_series)

    # get Li_{2,1,1} = Li_{0111} from Li_{111} = 1/6 Li_1^3
    Li111_series = 1/6 * Li_series[1]^3
    Li211_series = Lip211(ζ) + e0_integrate(Li111_series)

    return [
        log_series,
        Li_series[1],
        Li_series[2],
        Li_series[3],
        Li21_series,
        Li_series[4],
        Li31_series,
        Li211_series
    ]


def Z_one_sixth_periods(p,prec,li_data=None):
    r"""
    Values of the p-adic periods
      a3,a2,a32,a332,a322,aσ,a3332,a3322,a3222,a3σ,a2σ
    for Z[1/6] as listed in [CDCL26], Proposition 5.20.
    Here, "2" and "3" are shorthand for τ₂ and τ₃ and σ is shorthand for σ₃.

    INPUT:

        - ``p`` -- a prime ≥ 5

        - ``prec`` -- a positive integer specifying the p-adic precision

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.

    OUTPUT: the list of p-adic periods ``[a3,a2,a32,a332,a322,aσ,a3332,a3322,a3222,a3σ,a2σ]``

    """

    # compute polylog data if not provided
    if li_data is None:
        li_data = polylog_data(4,p,prec)

    logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211 = multiple_polylog_functions(p,prec,li_data)
    a2 = logp(2)
    a3 = logp(3)
    aσ = -4/3*Lip3(-1)  # ζ(3) = -4/3 Li_3(-1)
    a32 = -Lip2(3)
    a322 = -Lip3(-2) + logp(2)*Lip2(-2) + 1/2*logp(2)^2*logp(3)
    a332 = -Lip3(3)
    as2 = 4/7*Lip31(-1)
    a3σ = 18/13*Lip4(3) - 3/52*Lip4(9)
    a3332 = -Lip4(3)
    a3322 = - Lip4(2/3) - Lip4(-2) - (8/7)*(logp(2)^4/24 + Lip4(1/2)) \
        - (5/13)*Lip4(3) + (3/52)*Lip4(9) - (1/4)*logp(2)^2*logp(3)^2 \
        - logp(2)*Lip3(3) + (1/6)*logp(2)*logp(3)^3 - (1/24)*logp(3)^4
    a3222 = Lip4(-2) - logp(2)*Lip3(-2) + 1/2*logp(2)^2*Lip2(-2) \
        + 1/6*logp(2)^3*logp(3)
    a2σ = a2*aσ - as2  # shuffle relation

    return a3,a2,a32,a332,a322,aσ,a3332,a3322,a3222,a3σ,a2σ


def Z_one_sixth_kim_functions(p,prec,abstract_kim_functions, li_data=None):
    r"""
    The two p-adic Kim functions vanishing on the level 4 Chabauty—Kim locus X(ℤₚ)_{Z,4}
    for Z = Spec(ℤ[1/6]). It takes the abstract Kim functions for #S = 2 as input
    and computes the p-adic specialization for S = {2,3} using the periods from
    ``Z_one_sixth_periods``. Code for computing the abstract Kim functions is
    contained in ``abstract_kim_functions.sage``.

    INPUT:

        - ``p`` -- a prime ≥ 5

        - ``prec`` -- a positive integer specifying the p-adic precision

        - ``abstract_kim_functions`` - a list of polynomials in
            ℚ[a1,a2,a12,a112,a122,aσ,a1112,a1122,a1222,a1σ,a2σ][Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211]

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.
    
    OUTPUT: the list of polynomials in ℚₚ[Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211] obtained by specializing
        the coefficients with the p-adic periods from ``Z_one_sixth_periods``

    """

    if li_data is None:
        li_data = polylog_data(4,p,prec)
    periods = Z_one_sixth_periods(p,prec,li_data)
    return [ specialize_kim_function(F, periods) for F in abstract_kim_functions ]


def abstract_Z_one_over_ell_functions():
    r"""
    The four abstract Kim functions in level 4 for Z = Spec(ℤ[1/ℓ]),
    as derived in §5.1 of [CDCL26]. The functions are polynomials in
      ℚ[aτ, aσ, aτσ][Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211],
    with τ being shorthand for τ_ℓ, and σ being shorthand for σ₃.

    OUTPUT: the list of polynomials in ℚ[aτ, aσ, aτσ][Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211]
    from Theorem 5.2 of [CDCL26].

    """

    alphabet = ["τ","σ"]
    weights = [1,3]

    # weight-truncated ring of motivic periods as a polynomial ring in Lyndon variables
    A, a = shuffle_algebra_as_polynomial_ring(alphabet, weights, weight_bound=4)
    R_dR_mot.<Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211> = PolynomialRing(A)

    return [
        Li2 - 1/2 * Log*Li1,
        a("τ")*a("σ")*Li4 - a("τσ")*Log*Li3 - 1/24*(a("τ")*a("σ") - 4*a("τσ"))*Log^3*Li1,
        a("τ")*a("σ")*(Li31 + 1/8*Log^2*Li1^2) - a("τσ")*Log*Li21 - a("στ")*Li1*Li3 + a("σ")*a("στ")*Li1,
        a("τ")*a("σ")*Li211 - a("στ")*Li1*Li21 - 1/24*(a("τ")*a("σ") - 4*a("στ"))*Log*Li1^3 + a("σ")*a("στ")*Li1
    ]


def Z_one_half_periods(p,prec,li_data=None):
    r"""
    The values of the p-adic periods a_τ, a_σ, a_{τσ} ∈ A(Z) for Z = Spec(ℤ[1/2])
    as listed in Propositions 2.21 and 5.4 of [CDCL26]:
        a_τ = log(2),  a_σ = ζ(3),  a_{τσ} = 8/7 Li_4(2).
    
    INPUT:

        - ``p`` -- an odd prime

        - ``prec`` -- a positive integer specifying the p-adic precision

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.
    
    OUTPUT: the list of p-adic periods ``[a_τ, a_σ, a_{τσ}]`` for Z = Spec(ℤ[1/2])

    """

    if li_data is None:
        li_data = polylog_data(4,p,prec)
    g,k0 = li_data
    Kp = Qp(p,prec)
    aτ = log(Kp(2))
    aσ = -4/3*fast_polylog(-1,3,p,prec,g,k0) # ζ(3)
    aτσ = 8/7*fast_polylog(2,4,p,prec,g,k0)  # 8/7 Li_4(2)
    return [aτ, aσ, aτσ]


def Z_one_half_kim_functions(p,prec,li_data=None):
    r"""
    The four p-adic Kim functions vanishing on the level 4 Chabauty—Kim locus X(ℤₚ)_{Z,4}
    for Z = Spec(ℤ[1/2]). It takes the abstract Kim functions for #S = 1 from
    ``abstract_Z_one_over_ell_functions`` and computes their p-adic specialization 
    for S = {2} using the periods from ``Z_one_half_periods``. 

    INPUT:

        - ``p`` -- an odd prime

        - ``prec`` -- a positive integer specifying the p-adic precision

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.
    
    OUTPUT: the list of four polynomials in ℚₚ[Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211] from
        Theorem 5.2 of [CDCL26], specialized to the case Z = Spec(ℤ[1/2]).

    """

    if li_data is None:
        li_data = polylog_data(4,p,prec)
    periods = Z_one_half_periods(p,prec,li_data)
    abstract_kim_functions = abstract_Z_one_over_ell_functions()
    return [ specialize_kim_function(F_abs, periods) for F_abs in abstract_kim_functions ]


def Z_one_third_periods(p,prec,li_data=None):
    r"""
    The values of the p-adic periods a_τ, a_σ, a_{τσ} ∈ A(Z) for Z = Spec(ℤ[1/3])
    as listed in Propositions 2.21 and 5.8 of [CDCL26]:
        a_τ = log(3),  a_σ = ζ(3),  a_{τσ} = 18/13 Li_4(3) - 3/52 Li_4(9).
    
    INPUT:

        - ``p`` -- a prime ≠ 3

        - ``prec`` -- a positive integer specifying the p-adic precision

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.
    
    OUTPUT: the list of p-adic periods ``[a_τ, a_σ, a_{τσ}]`` for Z = Spec(ℤ[1/3])

    """

    Kp = Qp(p,prec)
    if p != 2:
        if li_data is None:
            li_data = polylog_data(4,p,prec)
        g,k0 = li_data
        Li4 = lambda z: fast_polylog(z,4,p,prec,g,k0)
        aτ = log(Kp(3))
        aσ = -4/3*fast_polylog(-1,3,p,prec,g,k0) # ζ(3)
        aτσ = 18/13*Li4(3) - 3/52*Li4(9)
    else:
        # fast_polylog only works for z ≢ 0,1 mod p; use built-in polylog method for p=2
        aτ = log(Kp(3))
        aσ = Kp(1).polylog(3)
        aτσ = 18/13*Kp(3).polylog(4) - 3/52*Kp(9).polylog(4)
    return [aτ, aσ, aτσ]


def Z_one_third_kim_functions(p,prec,li_data=None):
    r"""
    The four p-adic Kim functions vanishing on the level 4 Chabauty—Kim locus X(ℤₚ)_{Z,4}
    for Z = Spec(ℤ[1/3]). It takes the abstract Kim functions for #S = 1 from
    ``abstract_Z_one_over_ell_functions`` and computes their p-adic specialization 
    for S = {3} using the periods from ``Z_one_third_periods``. 

    INPUT:

        - ``p`` -- a prime ≠ 3

        - ``prec`` -- a positive integer specifying the p-adic precision

        - ``li_data`` (optional) -- precomputed output of ``polylog_data(4,p,prec)`` for efficiency.
    
    OUTPUT: the list of four polynomials in ℚₚ[Log,Li1,Li2,Li3,Li21,Li4,Li31,Li211] from
        Theorem 5.2 of [CDCL26], specialized to the case Z = Spec(ℤ[1/3]).

    """

    if li_data is None:
        li_data = polylog_data(4,p,prec)
    periods = Z_one_third_periods(p,prec,li_data)
    abstract_kim_functions = abstract_Z_one_over_ell_functions()
    return [ specialize_kim_function(F_abs, periods) for F_abs in abstract_kim_functions ]
