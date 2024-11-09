module legendreExpansion

using Polynomials
using QuadGK
export Polynomials, QuadGK

# submodule to build and evaluate Legendre Polynomials
module legendrePolynomials
using ..Polynomials
export lp_generator, leg_eval

# Function that generates Legendre polynomials up to nth degree
function lp_generator(n::Int)
    legPoly = Vector{Polynomial{Float64}}()
    P0 = Polynomial([1.0])
    push!(legPoly, P0)
    P1 = Polynomial([0.0, 1.0])
    push!(legPoly, P1)

    if n == 0
        return legPoly[1]
    elseif n == 1
        return legPoly
    else
        for k in 2:n
            Pk = ((2.0*k-1)* Polynomial([0.0, 1.0]) * P1 - (k-1) * P0) / k
            push!(legPoly, Pk)
            P0, P1 = P1, Pk
        end
        return legPoly
    end
end

# Function that evaluates Legendre polynomials at x
function leg_eval(x, legPoly)
    return map(p -> p(x), legPoly)
end
# end of submodule
end

# submodule for function expansion
module functionExpansion
using ..QuadGK
export leg_coefficients, series_expansion

function leg_coefficients(f, lP)
    coefficients = Float64[]
    for i in eachindex(lP)
        integrand(x) = f(x) * lP[i](x)
        integral,_ = quadgk(integrand, -1.0, 1.0)
        an = (2*(i-1) + 1) * integral / 2.0
        push!(coefficients, an)
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



    
