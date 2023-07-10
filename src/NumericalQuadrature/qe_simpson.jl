struct QESimpson <: QuadratureMethod end

function integration_weights!(weights::AbstractVector, x::AbstractVector, ::QESimpson)
    dx = mesh_derivative(x)
    weights .= zero(eltype(weights))
    n = length(x)
    for i in 2:(n - 1)
        fct = abs(mod(i, 2) - 2) * 2
        weights[i] = fct * dx[i]
    end
    if mod(n, 2) == 1
        weights[1] = dx[1]
        weights[n] = dx[n]
    else
        weights[1] = dx[1]
        weights[n - 1] = -dx[n - 1]
    end
    return weights .*= 1 / 3
end
