# Radius of the sphere
R = 10.0
# Number of particles
N = 27 
# x = cos(theta)
x = collect(range(-1.0, .999, length = 40001))
# Calculate distance r
r = R .* (2.0 .- 2.0 .* x).^(0.5)

lam = 3.0

sigma = 2 * R / lam

# we will now evaluate the Legendre Polynomials using our x values
include("legendreExpansion.jl")
using .legendreExpansion.legendrePolynomials
n = 50 
lP = lp_gen_eval(n, x)

# now we will define our yukawa potential
include("potentials.jl")
using .potentials
u = softSphere.(r, sigma)

# We will define our bridge
include("bridges.jl")
using .bridges
bridge = bridge_HNC(x)
println("Parameters succesfully defined")
println("Obtaining gamma for OZ equation on a monolayer") 
# solving Ornstein-Zernike equation

include("sphericalOZ.jl")
using .sphericalOZ
gamma_0 = zeros(length(x))

gamma = ng_solve(x, N, u, bridge, gamma_0, lP)

c = exp.(-u .+ gamma .+ bridge) .- gamma .- 1

g_theta = gamma .+ c .+ 1.0

g_theta = reverse(g_theta)  # Reverse the array
theta = reverse(acos.(x)) ./ π .* 180  # Reverse the array after computing acos(x) and converting to degrees

using Plots

plot(x, g_theta, label= "<N> = $N, lambda = $lam", xlabel = "theta", ylabel = "g(theta)", title = "Soft Sphere, HNC", dpi = 300)
savefig("soft_sphere_hnc.png")
#=
using DataFrames

c = exp.(-u .+ gamma_0 .+ bridge) .- gamma_0 .- 1   

h = gamma_0 .+ c

g = h .+ 1.0

h_coeff = leg_coeff_arr(x, h, lP)

g_coeff = leg_coeff_arr(x, g, lP)

df = DataFrame(gamma_l = gamma_coeff, h_l = h_coeff, g_l = g_coeff)
println(df)=#

