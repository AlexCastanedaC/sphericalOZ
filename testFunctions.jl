module testFunctions

export dist, inv_dist, yukawa

# r1 > r2
function dist(x, r1, r2)
    xi = r2 / r1
    return r1 * (1 + xi^2 - 2 * xi * x)^(1/2)
end

function yukawa(x, r1, r2, k)
    return exp(-k * dist(x, r1, r2) ) / dist(x, r1, r2)
end

function inv_dist(x, r1, r2)
    return 1 / dist(x, r1, r2)
end

# ending module
end
