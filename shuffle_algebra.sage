def shuffle_algebra_as_polynomial_ring(alphabet, weights, weight_bound):
    r"""
    Given a finite ordered alphabet Σ = {σ₁ < ... < σₙ} and positive weights w₁,...,wₙ,
    let A∞ be the shuffle algebra with basis a_w, with w running through words in Σ. 
    Define the weight grading by wt(a_σᵢ) := wᵢ. Let A be the sub-algebra generated 
    by a_w with wt(a_w) <= N. This is isomorphic to the polynomial algebra `\QQ[{a_u}]` 
    with `u` a Lyndon word and wt(a_u) <= N. This function constructs the polynomial
    ring A together with a function to convert an arbitrary word w of weight <= N
    to a_w ∈ A.

    INPUT:

    - ``alphabet`` -- a finite list of strings

    - ``weights`` -- a list of positive integers of the same length as the alphabet

    - ``weight_bound`` -- a nonnegative integer N

    OUTPUT: a pair consisting of

    - ``A`` -- the polynomial algebra over `\QQ` in variables ``"a"+w`` where ``w``
        is a Lyndon word in the alphabet of weight <= N
    
    - ``a`` -- a function that takes a word ``w`` of weight <= N and returns the element
        a_w of ``A``
    

    EXAMPLES:

        sage: alphabet = ["1","2","σ"]; weights = [1,1,3]
        sage: A,a = shuffle_algebra_as_polynomial_ring(alphabet,weights,weight_bound = 4); A
        Multivariate Polynomial Ring in a1, a2, a12, a112, a122, aσ, a1112, a1122, a1222, a1σ, a2σ over Rational Field
        sage: a("21")
        a1*a2 - a12
        sage: a("2112")
        -1/2*a12^2 + a2*a112
        sage: a("2") * a("12") == a("212") + 2*a("122")  # shuffle relation
        True
    

    AUTHORS:
        - Martin Lüdtke (2026-11-01)

    """

    Sh = ShuffleAlgebra(QQ, alphabet)
    W  = Words(alphabet)

    # generate all words of weight 0,...,N
    words_of_weight = [[W()]] + [[] for w in range(1,weight_bound+1)]
    for n in range(1,weight_bound+1):
        for x,wx in zip(alphabet,weights):
            if wx <= n:
                words_of_weight[n] += [w + W(x) for w in words_of_weight[n-wx]]
        words_of_weight[n].sort()
    all_words = [ w for n in range(weight_bound+1) for w in words_of_weight[n] ]

    # Construct polynomial ring on Lyndon generators
    lyndon_words = [ w for w in all_words if w.is_lyndon() ]
    lyndon_variable_names = [ "a" + str(w) for w in lyndon_words ]
    A = PolynomialRing(QQ, names=lyndon_variable_names)
    lyndon_variable = dict(zip(lyndon_words, A.gens()))  # Lyndon word -> Lyndon variable

    # construct dict to convert words w to polynomials in Lyndon variables
    word_to_poly = { }
    for w in all_words:
        LF = w.lyndon_factorization()   # w = l_1 ⋯ l_r with l_i Lyndon, l_1 >= ... >= l_r
        # a_l_1 ... a_l_r = sum c_v a_v with v >= w, solve for a_w
        shuffle_prod = Sh(1)
        for l in LF:
            shuffle_prod *= Sh(l)
        coeffs = shuffle_prod.monomial_coefficients()

        # a_w = (a_l_1 ... a_l_r - sum_{v ≠ w} c_v a_v) / c_w
        poly = prod(lyndon_variable[l] for l in LF)
        for v,c in coeffs.items():
            if v != w:
                poly -= c*word_to_poly[v]
        word_to_poly[w] = poly/coeffs[w]

    # converts a word w of weight <= N in the alphabet to the element a_w of A
    a = lambda w: word_to_poly[W(w)]

    return A, a
