module legendreExpansion

using Polynomials
using FastGaussQuadrature 
export Polynomials,FastGaussQuadrature
include("mathTools.jl")
export mathTools

# submodule to build and evaluate Legendre Polynomials
module legendrePolynomials
using ..Polynomials
export lp_generator, leg_eval, lp_gen_eval

# Function that generates Legendre polynomials up to nth degree
function lp_generator(n::Int)
    legPoly = Vector{Polynomial{BigFloat}}()
    P0 = Polynomial([BigFloat(1.0)])
    push!(legPoly, P0)
    P1 = Polynomial([BigFloat(0.0), BigFloat(1.0)])
    push!(legPoly, P1)
    x = P1

    if n == 0
        return legPoly[1]
    elseif n == 1
        return legPoly
    else
        for k in 2:n
	    c1 = BigFloat((2.0 * k - 1) / k)
	    c2 = BigFloat((k - 1.0) / k)
	    Pk = c1 * x * P1 - c2 * P0  
            push!(legPoly, Pk)
            P0, P1 = P1, Pk
        end
        return legPoly
    end
end

# Function that generates evaluated Legendre Polynomials, 
# assuming the user supplies n > 2 to avoid if evaluation
#=
function lp_gen_eval(n::Int, x::Vector)
    legPoly = Vector{Vector{Float64}}()
    P0 = ones(Float64, size(x))
    push!(legPoly, P0)
    P1 = x
    push!(legPoly, P1)

    for k in 2:n
	c1 = (2.0 * k - 1) / k
	c2 = (k - 1.0) / k
	Pk = c1 .* x .* P1 .- c2 .* P0
	push!(legPoly, Pk)
	P0, P1 = P1, Pk
    end
    return legPoly
end
=#
function lp_gen_eval(n::Int, x::Vector{Float64})
    legPoly = Vector{Vector{Float64}}(undef, n + 1)
    legPoly[1] = ones(Float64, length(x))
    legPoly[2] = x

    for k in 2:n
        c1 = (2.0 * k - 1) / k
        c2 = (k - 1.0) / k
        legPoly[k + 1] = c1 .* x .* legPoly[k] .- c2 .* legPoly[k - 1]
    end

    return legPoly
end
# Function that evaluates Legendre polynomials at x
function leg_eval(x, legPoly)
    return map(p -> p(x), legPoly)
end
# end of submodule
end

# submodule for function expansion
module functionExpansion
#using ..FastGaussQuadrature
using ..mathTools
export leg_coefficients, series_expansion, leg_coeff_arr

function leg_coefficients(f, lP)
    coefficients = Float64[]
   # x, w = gausslegendre(200)
    for i in eachindex(lP)
        integrand(x) = f(x) * lP[i](x)
	integral = simpson(integrand, -1.0, 1.0,100000)
        an = (2*(i-1) + 1) * integral / 2.0
        push!(coefficients, an)
    end
    return coefficients
end

# The following function creates the coefficients from arrays
function leg_coeff_arr(x::Vector, f::Vector, lP::Vector)
    coefficients = Float64[]
    for i in eachindex(lP)
	    integrand = f .* lP[i]
	    integral = simps_arr(integrand, x)
	    an = (2.0 * (i-1) + 1.0) * integral / 2.0
	    push!(coefficients,an)
    end
    return coefficients
end


function series_expansion(coefficients, eval_lp)
    return sum([coefficients[i] * eval_lp[i] for i in eachindex(coefficients)])
end
# end of submodule 
end

# end of legendreExpansion module
end



    
