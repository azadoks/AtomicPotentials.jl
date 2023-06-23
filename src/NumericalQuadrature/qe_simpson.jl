@doc raw"""
QuantumESPRESSO Simpson's (1/3) rule quadrature. The expression is equivalent to the
uniform grid case of `Simpson` for an even number of subintervals. However, the
derivative of the mesh `dx` is used in place of the finite difference between adjacent
mesh points `Δx`.
"""
struct QESimpson <: QuadratureMethod end

function integration_weights!(weights::AbstractVector, x::RadialMesh, ::QESimpson)
    weights .= zero(eltype(weights))
    n = length(mesh)
    for i in 2:(n - 1)
        fct = abs(mod(i, 2) - 2) * 2
        weights[i] = fct * deriv(mesh, i)
    end
    if mod(n, 2) == 1
        weights[1] = deriv(mesh, 1)
        weights[n] = deriv(mesh, n)
    else
        weights[1] = deriv(mesh, 1)
        weights[n - 1] = -deriv(mesh, n - 1)
    end
    return weights .*= 1 / 3
end
