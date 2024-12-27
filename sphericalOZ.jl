module sphericalOZ
export simplePicard, compute_D, mix_solve, ng_solve
# Defining simple Picard iteration for Ng Algorithm
# simplePicard returns the new gamma
include("legendreExpansion.jl")
using .legendreExpansion.functionExpansion


function simplePicard(x::Vector{Float64}, N::Int64, gamma::Vector{Float64}, u::Vector{Float64}, bridge::Vector{Float64}, legendrePol::Vector{Vector{Float64}})
    # Using closure 
    c_dir = exp.(-1.0 .* u .+ gamma .+ bridge) .- gamma .- 1.0
    
    # Calculating Legendre coefficients for c_dir
    c_l = leg_coeff_arr(x, c_dir, legendrePol)
    l = collect(0:length(c_l)-1)
    
    # Applying OZ equation to obtain gamma_l
    denominator = 2.0 .* l .+ 1.0
    fraction = (N ./ denominator) .* c_l
    adjusted_c_l = c_l ./ (1.0 .- fraction)
    
    # Final expression
    gamma_l = adjusted_c_l .- c_l
    gamma_new = series_expansion(gamma_l, legendrePol)
    
    return gamma_new
end

include("mathTools.jl")
using .mathTools 

# Define the function to compute corrections and Matrix D
function compute_D(d_n::Vector, x_grid::Vector)
	#Calculate the differences d
	n = length(d_n)
	d = [d_n[end] - d_n[end-i] for i in 1:(n-1)]

	# Calculate the matrix D
	D = [inner_product(d[i], d[j], x_grid) for i in 1:(n-1), j in 1:(n-1)]

	d_vec = [inner_product(d_n[end], d[i], x_grid) for i in 1:(n-1)]

	return D, d_vec
end

function mix_solve(x::Vector, N::Number, u::Vector, bridge::Vector, initial_gamma, legendrePol)

	tol = 10.0
	count = 0
	alpha = 0.1
	gamma_iter = initial_gamma
	
	while tol > 1e-5
		gamma = simplePicard(x, N, gamma_iter, u, bridge, legendrePol)
		gamma = alpha .* gamma .+ (1.0 - alpha) .* gamma_iter
		diff = gamma .- gamma_iter
		tol = sqrt(inner_product(diff, diff, x))
		gamma_iter = gamma
		count += 1
		println("count = $count, Tol = $tol")
	end

	return gamma_iter
end

function ng_solve(x::Vector, N::Number, u::Vector, bridge::Vector, initial_gamma, legendrePol)
gamma = Vector{Vector{Float64}}()
g = Vector{Vector{Float64}}()
push!(gamma, initial_gamma)
i = 0
tol_Ng = 10.0
l=0# Append new values to gamma using simple_picard
while tol_Ng > 1e-5
    for j in 1:3
        push!(gamma, simplePicard(x, N, gamma[i + j], u, bridge,legendrePol))	
    end
   
    # Extend g with the last three elements of gamma
    append!(g, gamma[end-2:end])
   # println(g[2])
    # Compute differences d_n
    d_n = [g[i + k] .- gamma[i + k] for k in 1:3] 
   # println(d_n)
    # Compute D and d_vec
    D, d_vec = compute_D(d_n, x)
    #println(D)

    # Solve for c using linear algebra
    c = D \ d_vec

    # Calculate sum_g
    sum_g = zeros(size(g[end]))
    for j in 1:2
        sum_g .+= c[j] .* g[end-j]
    end

    # Compute f_n1
    f_n1 = (1 - sum(c)) .* g[end] .+ sum_g

    # Calculate difference d_n2
    d_n2 = f_n1 .- g[end]

    # Update tolerance
    tol_Ng = sqrt(inner_product(d_n2, d_n2, x))

    # Update the last element of gamma
    gamma[end] = f_n1

    # Increment counters
    i += 3
    l += 1
println("l = $l, Tol = $tol_Ng")
end

return gamma[end]

end

# end of sphericalOZ module
end



