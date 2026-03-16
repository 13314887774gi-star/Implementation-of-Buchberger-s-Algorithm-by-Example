# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 21:26:23 2026

@author: wangyunong
"""

import sympy as sp
x, y = sp.symbols('x y')
##Define the set we are testing.
f1 = x**4*y**2
f2 = x**2*y**5+x*y**2
f3 = 3*x**2+x
F = [f1, f2, f3]

##Define the leading terms and leading monomials using graded lex ordering.
def lm_grlex(f):
    poly = sp.Poly(f, x, y)
    lm = max(poly.monoms(), key=lambda mon: (sum(mon), mon))
    return x**lm[0] * y**lm[1] 

def lt_grlex(f):
    poly = sp.Poly(f, x, y)
    lm = max(poly.monoms(), key=lambda mon: (sum(mon), mon))
    coeff = poly.coeff_monomial(x**lm[0] * y**lm[1] )
    return coeff * x**lm[0] * y**lm[1] 

##Find the least common multiple of given monomials
def lcm_mono (f,g):
    lcm = sp.lcm(lm_grlex(f), lm_grlex(g))
    return lcm
##Testing
print("The least common multiple is:", lcm_mono(f1,f2))

##Calculate the Spolynomials given two polynomials
def spoly(f,g):
    s = (lcm_mono(f, g) / lt_grlex(f))* f.as_expr() - (lcm_mono(f, g) / lt_grlex(g))* g.as_expr()
    h = sp.expand(s)
    return h
##Testing
print("The Spolynomial calculated is:", spoly(f1, f2))

#Compute the remainder of f on division of set G.
def compute_remainder(f, G):
    remainder = f
    for g in G:
          remainder = sp.rem(remainder ,g)
    return remainder

##Implement Bucherber's Algorithm
def buchberger(G):
    new_G = list(G)  
    pairs = [(new_G[i], new_G[j]) for i in range(len(new_G)) for j in range(i+1, len(new_G))]
    while pairs:
        f, g = pairs.pop(0)  
        s = spoly(f, g)
        remainder = compute_remainder(s, new_G)  
        if remainder != 0 and remainder not in new_G:  
            new_G.append(remainder)
            for g in new_G[:-1]:  
                pairs.append((remainder, g))

    return new_G


groebner_basis = buchberger(F)
print("Gröbner Basis:", groebner_basis)
                
                
##Test our result
def is_groebner_basis(G):
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            s = spoly(G[i], G[j])
            remainder = compute_remainder(s, G)
            if remainder != 0:
                print(f"remainder on S-poly of {G[i]} and {G[j]} is not zero: {remainder}")
                return False
    return True

if is_groebner_basis(buchberger(F)):
    print("The set is a Groebner basis.")
else:
    print("The set is not a Groebner basis.")
    
G = buchberger(F)    

##Find the corresponding reduced Groebner basis
def normalize_polynomial(f):
    
    poly = sp.Poly(f, x, y)
    lm = lm_grlex(f)
    lc = poly.coeff_monomial(lm)
    return f / lc if lc != 0 else f

def filter_polynomials(poly_list):

    filtered_list = poly_list[:]
    monomials = [lm_grlex(f) for f in poly_list]

    to_remove = set()
    for i in range(len(poly_list)):
        for j in range(len(poly_list)):
            if i != j and sp.rem(monomials[i], monomials[j]) == 0:
                to_remove.add(poly_list[i])
    
    filtered_normalized = [normalize_polynomial(f) for f in filtered_list if f not in to_remove]
    return filtered_normalized


filtered_G = filter_polynomials(G)
print("Reduced Groebner basis is:", filtered_G)


OUTPUT:
The least common multiple is: x**4*y**5
The Spolynomial calculated is: -x**3*y**2
Gröbner Basis: [x**4*y**2, x**2*y**5 + x*y**2, 3*x**2 + x, -x*y**2/9]
The set is a Groebner basis.
Reduced Groebner basis is: [x**2 + x/3, x*y**2]
