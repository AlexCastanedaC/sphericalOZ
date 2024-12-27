n = 100
x = collect(range(-1.0, 1.0, length = 100001))
include("legendreExpansion.jl")
include("testFunctions.jl")
include("mathTools.jl")

using .legendreExpansion.legendrePolynomials
lP = lp_gen_eval(n, x)

r1 = 1.0
r2 = 0.5
k = 1
# Let's work with inv_dist function 
using .testFunctions
#f1(x) = inv_dist(x,r1,r2)
#f(x) = x 
f2 = yukawa.(x, r1, r2, k)
#f3 = x
#f4 = 0.5.*(3.0 .* x.^2 .-1.0)
# Calculating coefficients
using .legendreExpansion.functionExpansion
#f1_coefficients = leg_coefficients(f1, lP)
f2_coefficients = leg_coeff_arr(x, f2, lP)
#f3_coefficients = leg_coeff_arr(x,f4,lP)
#theoretical_coefficients = [(1/r1)*(r2/r1)^l for l in 0:n]
using .mathTools

theoretical_coeff_2 = [k * (2/pi) * (2 * l + 1) * msbessel_i(l,k*r2) * msbessel_k(l,k*r1) for l in 0:n]

using DataFrames
df = DataFrame(CalculatedCoeff = f2_coefficients, TheoreticalCoeff = theoretical_coeff_2)

println(df)

#using Plots
#default(dpi=300)
#scatter(0:n, f1_coefficients, label= "Calculated Coeff.", xlabel = "Polynomial Order (l)", ylabel="Coeff. Value")
#scatter!(0:n, theoretical_coefficients, label = "Theoretical Coeff.", marker=:cross)
#savefig("Coeff_comparison.png")



