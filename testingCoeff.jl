n = 30

using .legendreExpansion.legendrePolynomials
lP = lp_generator(n)

r1 = 1.0
r2 = 0.5

# Let's work with inv_dist function 
using .testFunctions
f1(x) = inv_dist(x,r1,r2)

# Calculating coefficients
using .legendreExpansion.functionExpansion
f1_coefficients = leg_coefficients(f1, lP)

println(f1_coefficients)



