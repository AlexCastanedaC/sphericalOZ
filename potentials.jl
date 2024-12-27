module potentials

export screenedCoulomb, softSphere
# screened Coulomb potential

function screenedCoulomb(r, A::Number, kappa::Number, alpha::Number)
	return alpha * A * f_exp(alpha * r) * exp(-1.0 * kappa * r)
end

function f_exp(x)
	if x >= 1e-6
		return (1 - exp(-x))/x
	else
		taylor_f = 1 - 0.5 * x + (1/6) * x^2
		return taylor_f
	end
end

function softSphere(r, sigma)
	return (sigma / r)^6
end
# end module
end
