abstract type AbstractLocalPotential{S,A} <: AbstractAtomicQuantity{S,A} end
angular_momentum(quantity::AbstractLocalPotential)::Int = 0

## Numerical local potential
struct LocalPotential{S,Numerical} <: AbstractLocalPotential{S,Numerical}
    r::AbstractVector
    f::AbstractVector  # Vloc(r) in real-space; Vloc(q) in Fourier-space
    interpolator  # Vloc(r) in real-space; Vloc(q) in Fourier-Space
    Z
end
function (Vloc::LocalPotential{FourierSpace})(q::T)::T where {T}
    !iszero(q) && return Vloc.interpolator(q)  # Compensating charge background
    return zero(T)
end
function (Vloc::LocalPotential{RealSpace})(r::T)::T where {T}
    !iszero(r) && return Vloc.interpolator(r)  # Divergence at r=0
    return T(-Inf)
end
function fht(
    Vloc::LocalPotential{RealSpace},
    q::AbstractVector,
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    r²f = Vloc.r .* (Vloc.r .* Vloc.f .- -Vloc.Z)  # == r² (Vloc - -Z/r)
    F = fht(Vloc.r, r²f, q, angular_momentum(Vloc), method) .+ 4π .* (-Vloc.Z ./ q .^ 2)
    interpolator = Interpolation.construct_interpolator(
        q[(begin + 1):end], F[(begin + 1):end], interpolation_method
    )
    return _construct_dual_quantity(Vloc; r=q, f=F, interpolator=interpolator)
end
function ifht(
    Vloc::LocalPotential{FourierSpace},
    r::AbstractVector,
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    q²F = q .^ 2 .* Vloc.f .- -Vloc.Z  # == q² (Vloc - -Z/q²)
    f = fht(Vloc.r, q²F, r, angular_momentum(Vloc), method) .+ 4π / (2π)^3 .* (-Vloc.Z ./ r)
    interpolator = Interpolation.construct_interpolator(
        r[(begin + 1):end], f[(begin + 1):end], interpolation_method
    )
    return _construct_dual_quantity(Vloc; r=r, f=f, interpolator=interpolator)
end

## Coulomb local potential
struct CoulombLocalPotential{S,Analytical} <: AbstractLocalPotential{S,Analytical}
    Z
end
function (Vloc::CoulombLocalPotential{RealSpace})(r)
    return -Vloc.Z / r
end
function (Vloc::CoulombLocalPotential{FourierSpace})(q::T) where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    # General atom => Use default Coulomb potential
    # We use int_{R^3} -Z/r e^{-i q⋅x} = 4π / |q|^2
    return -4T(π) * Vloc.Z / q^2
end

## Gaussian local potential
struct GaussianLocalPotential{S,Analytical} <: AbstractLocalPotential{S,Analytical}
    α
    L
end
function (Vloc::GaussianLocalPotential{RealSpace})(r::T)::T where {T}
    return -Vloc.α / (sqrt(2T(π)) * Vloc.L) * exp(-(r / Vloc.L)^2 / 2)
end
function (Vloc::GaussianLocalPotential{FourierSpace})(q::T)::T where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    # = ∫_ℝ³ V(x) exp(-ix⋅q) dx
    return -Vloc.α * exp(-(q * Vloc.L)^2 / 2)
end

## Cohen-Bergstresser local potential
struct CohenBergstresserLocalPotential{FourierSpace,Analytical} <:
       AbstractLocalPotential{FourierSpace,Analytical}
    V_sym  # Map |G|^2 (in units of (2π / lattice_constant)^2) to form factors
    lattice_constant  # Lattice constant (in Bohr) which is assumed
end
function (Vloc::CohenBergstresserLocalPotential{FourierSpace})(q::T)::T where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    # Get |q|^2 in units of (2π / lattice_constant)^2
    qsq_pi = Int(round(q^2 / (2π / el.lattice_constant)^2; digits=2))
    return T(get(el.V_sym, qsq_pi, 0.0))
end
