module mathTools

using SpecialFunctions

export msbessel_i, msbessel_k, simpson, simps_arr, inner_product

function simpson(f::Function, a::Number, b::Number, n::Number)
	h = (b-a)/n
	s = f(a) + f(b)
	for i in 1:2:n-1
		s += 4 * f(a + i * h)
	end
	for i in 2:2:n-2
		s += 2 * f(a + i * h)
	end
	return h/3 * s
end

function simps_arr(y::Vector, x::Vector)
	n = length(y) - 1
	n % 2 == 0 || error(" y length (number of intervals) must be odd")
	length(x)-1 == n || error("x and y length must be equal")
	h = (x[end]-x[1])/n
	s = y[1] + 4 * sum(y[2:2:end-1]) + 2 * sum(y[3:2:end-2]) + y[end]
	return h/3.0 *s
end

# define inner product function
function inner_product(f1::Vector, f2::Vector, x::Vector)
	integrand = f1 .* f2
	return simps_arr(integrand, x)
end

# Define modified spherical Bessel functions of the first kind
function msbessel_i(n, x)
	return sqrt(pi/(2*x))*besseli(n+0.5,x)
end

# Define modified spherical Bessel functions of the third kind
function msbessel_k(n, x)
	return sqrt(pi/(2*x))*besselk(n+0.5,x)
end


# end module
end
