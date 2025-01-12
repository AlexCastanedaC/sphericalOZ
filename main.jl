# Radius of the sphere
R = 100.0
# Number of particles
N = 10  
# x = cos(theta)
x = collect(range(-1.0, 1.0, length = 40001))
# Calculate distance r
r = R .* (2.0 .- 2.0 .* x).^(0.5)

# we will now evaluate the Legendre Polynomials using our x values
include("legendreExpansion.jl")
using .legendreExpansion.legendrePolynomials
n = 100 
lP = lp_gen_eval(n, x)

# now we will define our yukawa potential
include("potentials.jl")
using .potentials
A = 2000.0
kappa = 1.0 / 288.0
alpha = 10.0
u = screenedCoulomb.(r, A, kappa, alpha)

# We will define our bridge
include("bridges.jl")
using .bridges
bridge = bridge_HNC(x)
println("Parameters succesfully defined")
println("Obtaining gamma for OZ equation on a monolayer") 
# solving Ornstein-Zernike equation
# Parameter for gradually turning potential on

function extract_floats_from_file(filename)
    # Open the file for reading
    open(filename, "r") do file
        # Skip the first six rows
        for _ in 1:6
            readline(file)
        end
        # Read the remaining lines and parse them as floats
        floats = Float64[]
        for line in eachline(file)
            push!(floats, parse(Float64, strip(line)))
        end
        return floats
    end
end

include("sphericalOZ.jl")
using .sphericalOZ
gamma_0 = zeros(length(x))
#gamma_0_coeff = extract_floats_from_file("N005.dat")
#print(gamma_0_coeff)
w = collect(range(0.1,1,length=23))

for (i,w_val) in pairs(w)
gamma = ng_solve(x, N, w_val * u, bridge,gamma_0, lP)
println("Factor $i")
global gamma_0 = gamma
end

using .legendreExpansion.functionExpansion
#gamma_0 = series_expansion(gamma_0_coeff, lP) 

#gamma = ng_solve(x, N, u, bridge, gamma_0, lP)

#include("mathTools.jl")
#using .mathTools
#diff = gamma .- gamma_0

#tol = inner_product(diff, diff, x)
#println("Tol = $tol")



using DataFrames

c = exp.(-u .+ gamma_0 .+ bridge) .- gamma_0 .- 1   

h = gamma_0 .+ c

g = h .+ 1.0
gamma_coeff = leg_coeff_arr(x, gamma_0, lP)

h_coeff = leg_coeff_arr(x, h, lP)

g_coeff = leg_coeff_arr(x, g, lP)

df = DataFrame(gamma_l = gamma_coeff, h_l = h_coeff, g_l = g_coeff)
println(df)

