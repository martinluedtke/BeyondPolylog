# BeyondPolylog
Computing Kim functions and verifying Kim's Conjecture beyond the polylogarithmic quotient

Sage code for the paper "Mixed Tate Chabauty–Kim beyond the polylog quotient" [[CDCL26]](https://arxiv.org/abs/2610.25365)

The files "Zproots.sage" and "RefinedCK.sage" from [https://github.com/martinluedtke/RefinedCK](https://github.com/martinluedtke/RefinedCK) need to be placed in the working directory to run the code. The code was tested with Sage 10.10.

## 1. Construction of Kim functions for \#S=2

The command `sage abstract_kim_functions.sage` constructs the two motivic Kim functions $F_1^{\mathrm{mot}}$ and $F_2^{\mathrm{mot}}$ from Theorem 5.13 of [CDC26]. These belong to the level 4 Chabauty–Kim ideal $\mathcal{I}\_{Z,4}^{\mathrm{mot}}$ for $Z = \mathrm{Spec}(\mathbb{Z}[1/S])$ with $`\#S=2`$. They are elements of 
```math
A(Z)[\mathrm{log},\mathrm{Li}_1,\mathrm{Li}_2,\mathrm{Li}_3,\mathrm{Li}_{2,1},\mathrm{Li}_4,\mathrm{Li}_{3,1}, \mathrm{Li}_{2,1,1}],
```
of half-weights 10 and 12, respectively, with coefficients in the ring of motivic periods
```math
\mathbb{Q}[a_1, a_2, a_{12}, a_{\sigma}, a_{112}, a_{122}, a_{1112}, a_{1122}, a_{1222}, a_{1\sigma}, a_{2\sigma}] \subseteq A(Z),
```
see §5.5 of the paper for the notation. For any prime $p \not \in S$, the specialization of $F_1^{\mathrm{mot}}$ and $F_2^{\mathrm{mot}}$ along the $p$-adic period homomorphism $\mathrm{per}_p\colon A(Z) \to \mathbb{Q}_p$ yields two Kim functions whose vanishing locus in $`X(\mathbb{Z}_p)`$ contains the level 4 Chabauty–Kim locus $`X(\mathbb{Z}_p)_{Z,4}`$.

The output of the script is saved in `abstract_kim_functions.txt` in human-readable form, and in `abstract_kim_functions.sobj` in Sage-readable form, to be used in `kims_conj_Z_one_sixth.sage` for verifying Kim's Conjecture in the case $`Z = \mathrm{Spec}\, \mathbb{Z}[1/6]`$.

The command `sage weight_bound.sage` is used to prove upper bounds for the dimension of the Chabauty–Kim ideal in a given half-weight. It produces the upper bound
```math
\dim\left(\mathcal{I}_{Z,4}^{\mathrm{mot}} \cap (\mathrm{hwt} \leq 10)\right) \leq 1,
```
which proves Theorem 5.14, stating that $`F_1^{\mathrm{mot}}`$ is the unique motivic Kim function in half-weight $\leq 10$ up to scaling. In particular, there does not exist a Kim function in lower weight.

## 2. Verifying Kim's Conjecture

The scripts `kims_conj_Z_one_half.sage`, `kims_conj_Z_one_third.sage`, `kims_conj_Z_one_sixth.sage` were used to prove Theorem 1.2, i.e., to verify Kim's Conjecture in level 4 for all primes $p < 1000$, $p \in Z$, in the following cases:

  * $`Z = \mathrm{Spec}\,\mathbb{Z}[1/2]`$;
  * $`Z = \mathrm{Spec}\,\mathbb{Z}[1/3]`$;
  * $`Z = \mathrm{Spec}\,\mathbb{Z}[1/6]`$.

More precisely, in the first two cases something slightly stronger is shown: we have inclusions
```math
X(\mathbb{Z}_p) \supseteq X(\mathbb{Z}_p)_{Z,\mathrm{PL},4} \supseteq X(\mathbb{Z}_p)_{Z,\mathrm{PL},4}^{S_3} \supseteq X(\mathbb{Z}_p)_{Z,4}  \supseteq X(Z)
```
with increasingly strong variants of the Chabauty–Kim method: polylogarithmic, S_3-symmetrized polylogarithmic, and full level 4. For $`Z = \mathrm{Spec}\,\mathbb{Z}[1/2]`$, it is verified that Kim's Conjecture holds already for the polylogarithmic level-4 locus:
```math
X(\mathbb{Z}_p)_{Z,\mathrm{PL},4} = \{2,-1,\tfrac12\} = X(\mathbb{Z}[1/2]).
```
For $`Z = \mathrm{Spec}\,\mathbb{Z}[1/3]`$, it is verified that the polylogarithmic Chabauty–Kim locus equals $`\{-1\}`$, which implies that Kim's Conjecture holds for the $`S_3`$-symmetrized polylogarithmic Chabauty–Kim locus:
```math
X(\mathbb{Z}_p)_{Z,\mathrm{PL},4}^{S_3} = \emptyset =  X(\mathbb{Z}[1/3]).
```
For $`Z = \mathrm{Spec}\,\mathbb{Z}[1/6]`$, the full level 4 quotient is needed and it is verified that
```math
X(\mathbb{Z}_p)_{Z,4} = \left\{2, \frac12, -1, 3, \frac13, \frac23, \frac32, -\frac12,-2, 4, \frac14, \frac34, \frac43, -\frac13, -3, 9, \frac19, \frac89, \frac98, -\frac18, -8 \right\} =  X(\mathbb{Z}[1/6]).
```
This uses the two Kim functions produced by `abstract_kim_functions.sage`.

## 3. Multiple polylogarithms and Chabauty–Kim loci in level 4

The file `BeyondPolylog.sage` contains several auxiliary functions for computing $p$-adic multiple polylogarithms and Chabauty–Kim loci in level 4. These can be used to verify instances of Kim's Conjecture as above. 

An example for computing the $p$-adic multiple polylogarithms $`\mathrm{Li}_{2,1}(z)`$, $`\mathrm{Li}_{3,1}(z)`$, $`\mathrm{Li}_{2,1,1}(z)`$:
```sage
sage: p = 7; prec = 15
sage: logp,Lip1,Lip2,Lip3,Lip21,Lip4,Lip31,Lip211 = multiple_polylog_functions(p,prec)
sage: Lip21(10)
7^3 + 6*7^5 + 3*7^6 + 5*7^7 + 5*7^8 + 2*7^9 + 6*7^10 + 7^11 + 2*7^12 + 5*7^13 + 6*7^14 + O(7^15)
sage: Lip31(10)
4*7^5 + 2*7^6 + 7^7 + 7^8 + 4*7^9 + 3*7^10 + 3*7^12 + 5*7^13 + 4*7^14 + O(7^15)
sage: Lip211(10)
3*7^4 + 7^5 + 4*7^6 + 5*7^7 + 3*7^8 + 7^9 + 3*7^10 + 4*7^11 + 7^12 + 2*7^13 + 5*7^14 + O(7^15)

```

Verifying the vanishing of Kim functions on $`X(\mathbb{Z}[1/2]) = \{2,-1,\tfrac12\}`$:
```sage
sage: p = 11; prec = 10
sage: f1,f2,f3,f4 = Z_one_half_kim_functions(p,prec)
sage: f3(2)
O(11^14)
sage: f4(-1)
O(11^14)
sage: f2(1/2)
O(11^14)
```

Verifying that the level 4 Chabauty—Kim locus for $`Z = \mathrm{Spec}\,\mathbb{Z}[1/3]`$ is empty:
```sage
sage: p = 13; prec = 10
sage: li_data =  polylog_data(4,p,prec)
sage: kim_functions = Z_one_third_kim_functions_poly(p,prec,li_data)
sage: compute_depth4_locus(kim_functions,li_data)
[]
```
With only the first two—polylogarithmic—Kim functions, the locus still contains $-1$:
```sage
sage: compute_depth4_locus(kim_functions[:2],li_data)
[12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + 12*13^7 + 12*13^8 + O(13^9)]
```

Verifying Kim's Conjecture for $`Z = \mathrm{Spec}\,\mathbb{Z}[1/6]`$:
```sage
sage: p = 17; prec = 10
sage: abstract_kim_functions = load("abstract_kim_functions.sobj")
sage: li_data =  polylog_data(4,p,prec)
sage: kim_functions = Z_one_sixth_kim_functions_poly(p,prec,abstract_kim_functions,li_data)
sage: locus = compute_depth4_locus(kim_functions,li_data)
sage: len(locus) == 21
True
sage: s_unit_sols = [
    2, -1, 1/2, 
    3, 1/3, 2/3, 3/2, -1/2, -2, 
    4, 1/4, 3/4, 4/3, -1/3, -3, 
    9, 1/9, 8/9, 9/8, -1/8, -8
]
sage: all(z in s_unit_sols for z in locus)
True
```
