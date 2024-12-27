include("testFunctions.jl")
using .testFunctions

r_layers = [1.0, 2.0]
num_per_layer = [5.0, 10.0]

# number of spherical layers
n = length(r_layers)
# x \in [-1.0, 1.0]

x = collect(range(-1.0, 1.0, length = 60001))
r = Array{Float64, 3}(undef, n, n, length(x))

using LinearAlgebra
N = diagm(num_per_layer)

for k in eachindex(x)
	for i in 1:n
		for j in 1:n
			r[i, j, k] = dist(x[k], r_layers[i], r_layers[j])
		end
	end
end

include("potentials.jl")
using .potentials

A = 2000.0
kappa = 1.0 / 288.0
alpha = 10.0
u = similar(r)
u = screenedCoulomb.(r, A, kappa, alpha)

include("bridges.jl")
using .bridges
bridge = bridge_HNC(r)

# we will now evaluate the Legendre Polynomials using our x values
include("legendreExpansion.jl")
using .legendreExpansion.legendrePolynomials
n_lP = 100
lP = lp_gen_eval(n_lP, x)
println("LegPol evaluated")

gamma_0 = zeros(size(r)) .+ 0.01 



using .legendreExpansion.functionExpansion
include("mathTools.jl")
using .mathTools
function simpleIter(n::Number, n_lP::Number,x::Vector{Float64}, N::Matrix{Float64}, gamma_0::AbstractArray, u::AbstractArray, bridge::AbstractArray, lP::AbstractArray)
	c_dir = exp.(-1.0 .* u .+ gamma_0 .+ bridge) .- gamma_0 .- 1.0
	c_l_coeff = Vector{Vector{Float64}}() 
	
	for i in 1:n
		for j in 1:n
			c_l = leg_coeff_arr(x, c_dir[i,j,:], lP)
			push!(c_l_coeff, c_l)
		end
	end
	
	c_l_matrix = reshape(vcat(c_l_coeff...), n, n, n_lP+1)
	gamma = similar(c_l_matrix)
	Iden = diagm(ones(n))
	
	for l in 0:n_lP
		den_l = 1.0 / (2.0 * l + 1)
		c = c_l_matrix[:,:,l+1]
		global gamma = inv(Iden - den_l * c*N)*c - c
	end

	gamma_new = Vector{Vector{Float64}}()

	for i in 1:n
		for j in 1:n
			gamma_ij = series_expansion(gamma[i,j,:],lP)
			push!(gamma_new, gamma_ij)
		end
	end
	
	return reshape(vcat(gamma_new...), n, n, length(gamma_new[1]))
end

#i Solution with mixing parameter

tol = 100.0
count = 0
alpha = 0.1

gamma_iter = gamma_0

while tol > 1e-6
	local gamma = simpleIter(n, n_lP,x, N, gamma_iter,0.1* u, bridge, lP)
	
	local gamma = alpha .* gamma .+ (1.0 - alpha) .* gamma_iter

	println(gamma[:,:,1])
	diff = gamma .- gamma_iter

	tol1 = 0.0

	for i in 1:n
		for j in 1:n
			tol1 += inner_product(diff[i,j,:], diff[i,j,:],x) 
		end
	end

	global tol = sqrt(tol1)
	global gamma_iter = gamma
	global count += 1
	println("count = $count, tol = $tol")
end

println("Convergence reached")
println(gamma_iter[1,1,:])
#println(gamma_iter[1, 2, :])
