include("testFunctions.jl")
using .testFunctions

r_layers = [100.0, 150.0]
num_per_layer = [10.0, 15.0]

# number of spherical layers
n = length(r_layers)
# x \in [-1.0, 1.0]

x = collect(range(-1.0, 1.0, length = 50001))
r = Array{Float64, 3}(undef, n, n, length(x))

using LinearAlgebra
N = diagm(num_per_layer)


for i in 1:n
	for j in 1:n
		r[i, j, :] = dist.(x, r_layers[i], r_layers[j])
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
n_lP = 120 
lP = lp_gen_eval(n_lP, x)
println("LegPol evaluated")

using .legendreExpansion.functionExpansion
include("mathTools.jl")
using .mathTools

function simpleIter(n::Number, n_lP::Number,x::Vector{Float64}, N::Matrix{Float64}, gamma_0::AbstractArray, u::AbstractArray, bridge::AbstractArray, lP::AbstractArray)
	c_dir = exp.(-1.0 .* u .+ gamma_0 .+ bridge) .- gamma_0 .- 1.0
        c_l_coeff = [leg_coeff_arr(x, c_dir[i, j, :], lP) for i in 1:n, j in 1:n]
	c_l_matrix = zeros(n, n, n_lP + 1)
	for i in 1:n, j in 1:n
		c_l_matrix[i, j, :] .= c_l_coeff[i, j]
	end	

	gamma_l = similar(c_l_matrix)
	Iden = I(n) 
	
	for l in 0:n_lP
		den_l = 1.0 / (2.0 * l + 1)
		c = c_l_matrix[:,:,l+1]
		gamma_l[:,:,l+1] =(inv(Iden - den_l.* c*N) * c) - c
	end

	gamma_new_2d = [series_expansion(gamma_l[i, j, :], lP) for i in 1:n, j in 1:n]
	gamma_new = similar(c_dir)

	for i in 1:n, j in 1:n
		gamma_new[i,j,:] .= gamma_new_2d[i, j]
	end
	return gamma_new
end

function inner_product_M(A::Array, B::Array,x::Vector)
	shape_Array = size(A)
	n = shape_Array[1]
	sum = 0.0
	for i in 1:n
		for j in 1:n
                        sum += inner_product(A[i,j,:], B[i,j,:],x) 
                end
        end
	return sum
end

function alpha_solve(n::Int64,n_lP::Int64,x::Vector{Float64}, N::Matrix{Float64}, u::AbstractArray, bridge::AbstractArray, initial_gamma::AbstractArray, legendrePol)
	tol = 100.0
	count = 0
	alpha = 0.1

	gamma_iter = initial_gamma

	while tol > 1e-6
		gamma = simpleIter(n, n_lP,x, N, gamma_iter, u, bridge, legendrePol)
		gamma = alpha .* gamma .+ (1.0 - alpha) .* gamma_iter

		diff = gamma .- gamma_iter

		tol = sqrt(inner_product_M(diff, diff, x))

		gamma_iter = gamma

		count += 1

		println("count = $count, Tol = $tol")
	end

	return gamma_iter

end

# Define the function to compute corrections and Matrix D
function compute_D(d_n::Vector, x::Vector)
	#Calculate the differences d
	n = length(d_n)
	d = [d_n[end] - d_n[end-i] for i in 1:(n-1)]

	# Calculate the matrix D
	D = [inner_product_M(d[i], d[j], x) for i in 1:(n-1), j in 1:(n-1)]

	d_vec = [inner_product_M(d_n[end], d[i], x) for i in 1:(n-1)] 
	#println(D)
	#println(d_vec)

	return D, d_vec
end

function ng_solve(n::Int64,n_lP::Int64,x::Vector{Float64}, N::Matrix{Float64}, u::AbstractArray, bridge::AbstractArray, initial_gamma::AbstractArray, legendrePol)
	gamma = Vector{Array{Float64}}()
	g = Vector{Array{Float64}}()
	push!(gamma, initial_gamma)
	i = 0
	tol_Ng = 1000.0
	l=0
	# Append new arrays to gamma using simpleIter
	while tol_Ng > 1e-6
		for j in 1:5
			push!(gamma,simpleIter(n, n_lP,x, N, gamma[i + j],u, bridge, legendrePol))
		end
	# Extend g with the last three elements of gamma
	append!(g, gamma[end-4:end])
    	# Compute differences d_n
	d_n = [g[i + k] .- gamma[i + k] for k in 1:5]
	# println(d_n[2][1,1,:])
	# Compute D and d_vec
	D, d_vec = compute_D(d_n, x)
	# Solve for c using linear algebra
	c = D \ d_vec
	# Calculate sum_g
	sum_g = zeros(size(g[end]))
	for j in 1:4
		sum_g .+= c[j] .* g[end-j]
	end
	# Compute f_n1
	f_n1 = (1 - sum(c)) .* g[end] .+ sum_g
	# Calculate difference d_n2
	d_n2 = f_n1 .- g[end]
	# Update tolerance
	tol_Ng = sqrt(inner_product_M(d_n2, d_n2, x))
	# Update the last element of gamma
	gamma[end] = f_n1
	# Increment counters
	i += 5
	l += 1
	println("l = $l, Tol = $tol_Ng")
        end

	return gamma[end]

end

gamma_0 = zeros(size(r)) 
 
w = collect(range(0.1,1,length=30))

for (i,w_val) in pairs(w)
gamma = ng_solve(n, n_lP,x, N, w_val .* u, bridge,gamma_0, lP)
println("Factor $i")
global gamma_0 = gamma
end

c_f = exp.(-u .+ gamma_0 .+ bridge) .- gamma_0 .- 1.0
g_theta = gamma_0 .+ c_f .+ 1.0
for i in 1:n, j in 1:n
	reverse!(g_theta[i,j,:])
end

#reverse!(x)
theta = acos.(x) ./ pi .* 180.0

using Plots

plot()

for i in 1:n, j in 1:n
	plot!(theta,g_theta[i,j,:]; ls=:auto, label="g$i$j")
end

plot!(;plot_title = "PCF", legend=:topright,dpi = 300)

savefig("pcf_2layers.png")

#=
using DataFrames

c_0 = exp.(-u[1,2,:] .+ gamma_0[1,2,:] .+ bridge[1,2,:]) .- gamma_0[1,2,:] .- 1

h_0 = gamma_0[1,2,:] .+ c_0

g_0 = h_0 .+ 1.0
gamma_coeff = leg_coeff_arr(x, gamma_0[1,2,:], lP)

h_coeff = leg_coeff_arr(x, h_0, lP)

g_coeff = leg_coeff_arr(x, g_0, lP)

df = DataFrame(gamma_l = gamma_coeff, h_l = h_coeff, g_l = g_coeff)
println(df)
=#
