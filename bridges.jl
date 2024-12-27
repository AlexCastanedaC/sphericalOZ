module bridges

export bridge_HNC, bridge_PY

# Function for the HNC condition
function bridge_HNC(x::AbstractArray)
    return zeros(size(x))
end

# Function for the PY condition
function bridge_PY(gamma_func::AbstractArray)
    return log.(1.0 .+ gamma_func) .- gamma_func
end
# end of module 
end
