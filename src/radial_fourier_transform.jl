import .NumericalQuadrature: QuadratureMethod, Simpson, integration_weights

@doc raw"""
Compute the radial Fourier transform via numerical quadrature.

Computes the radial Fourier transform
```math
F(q) = 4\pi \int_{0}^{\infty} r f(r) j_l(qr) r dr,
```
of an arbitrarily-spaced periodic sequence using numerical quadrature, where $j_l$ is
the spherical Bessel function of order $l$. The index $l$ may be any non-negative
integer.

See also: [`irft`](@ref)
"""
function rft(
    q::AbstractVector{TQ},
    x::AbstractVector{TX},
    f::AbstractVector{TX},
    l::Int;
    n_x_factors::Int=0,
    quadrature_method::QuadratureMethod=Simpson(),
    kwargs...,
) where {TQ<:Real,TX<:Real}
    weights_ = integration_weights(x, quadrature_method)
    integrand_ = similar(x, TQ)
    jl = fast_sphericalbesselj(l)
    x²f = x .^ (2 - n_x_factors) .* f
    return map(q) do qᵢ
        integrand_ .= x²f .* jl.(qᵢ .* x)
        return 4TQ(π) * dot(weights_, integrand_)
    end
end

@doc raw"""
Compute the inverse radial Fourier transform via numerical quadrature.

Computes the inverse radial Fourier transform
```math
f(r) = \frac{1}{(2\pi)^3} \int_{0}^{\infty} q f(q) j_l(rq) q dr,
```
of an arbitrarily-spaced periodic sequence using numerical quadrature, where $j_l$ is
the spherical Bessel function of order $l$. The index $l$ may be any non-negative
integer.

See also: [`rft`](@ref)
"""
function irft(
    r::AbstractVector{TR}, x::AbstractVector{TX}, F::AbstractVector{TX}, l::Int; kwargs...
) where {TR<:Real,TX<:Real}
    return 1 / (2TR(π))^3 .* rft(r, x, F, l; kwargs...)
end
