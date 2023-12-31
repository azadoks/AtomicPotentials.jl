import .NumericalQuadrature: QuadratureMethod, Simpson, integration_weights

abstract type AbstractLocalPotential{S,A} <: AbstractAtomicQuantity{S,A} end
angular_momentum(Vloc::AbstractLocalPotential) = 0
charge_ionic(Vloc::AbstractLocalPotential) = Vloc.Z
energy_correction(T::Type{<:Real}, ::AbstractLocalPotential) = zero(T)

## Numerical local potential
struct LocalPotential{S,Numerical} <: AbstractLocalPotential{S,Numerical}
    r::AbstractVector
    f::AbstractVector  # Vloc(r) in real-space; Vloc(q) in Fourier-space
    interpolator::BSplineKit.SplineWrapper  # Vloc(r) in real-space; Vloc(q) in Fourier-Space
    Z::Real
end

function (Vloc::LocalPotential{FourierSpace})(q::T)::T where {T}
    !iszero(q) && return Interpolation.evaluate(Vloc.interpolator, q)  # Compensating charge background
    return zero(T)
end

function (Vloc::LocalPotential{RealSpace})(r::T)::T where {T}
    !iszero(r) && return Interpolation.evaluate(Vloc.interpolator, r)  # Divergence at r=0
    return T(-Inf)
end

function energy_correction(
    T::Type{<:Real},
    Vloc::LocalPotential{RealSpace},
    quadrature_method::QuadratureMethod=Simpson(),
)
    integrand = Vloc.r .* (Vloc.r .* Vloc.f .- -Vloc.Z)
    weights = integration_weights(Vloc.r, quadrature_method)
    return 4T(π) * dot(weights, integrand)
end

## HGH local potential
struct HghLocalPotential{S,Analytical} <: AbstractLocalPotential{S,Analytical}
    r
    c
    Z
end

# [GTH98] (1)
function (Vloc::HghLocalPotential{RealSpace})(r::T)::T where {T<:Real}
    r += iszero(r) ? eps(T) : zero(T)  # quick hack for the division by zero below
    rr::T = r / Vloc.r
    c = Vloc.c
    return -Vloc.Z / r * erf(rr / sqrt(T(2))) +
           exp(-rr^2 / 2) * (c[1] + c[2] * rr^2 + c[3] * rr^4 + c[4] * rr^6)
end

# [GTH98] (6) except they do it with plane waves normalized by 1/sqrt(Ω).
function (Vloc::HghLocalPotential{FourierSpace})(q::T)::T where {T<:Real}
    iszero(q) && return zero(T)
    x::T = q * Vloc.r
    return _hgh_local_potential_polynomial_fourier(Vloc, x) * exp(-x^2 / 2) / x^2
end

# The local potential of a HGH pseudopotentials in reciprocal space
# can be brought to the form ``Q(t) / (t^2 exp(t^2 / 2))``
# where ``x = r_\text{loc} q`` and `Q`
# is a polynomial of at most degree 8. This function returns `Q`.
@inline function _hgh_local_potential_polynomial_fourier(
    Vloc::HghLocalPotential{FourierSpace}, x::T
) where {T<:Real}
    r::T = Vloc.r
    Z::T = Vloc.Z

    # The polynomial prefactor P(t) (as used inside the { ... } brackets of equation
    # (5) of the HGH98 paper)
    P = (
        Vloc.c[1] +
        Vloc.c[2] * (3 - x^2) +
        Vloc.c[3] * (15 - 10x^2 + x^4) +
        Vloc.c[4] * (105 - 105x^2 + 21x^4 - x^6)
    )

    return 4T(π) * r^2 * (-Z + sqrt(T(π) / 2) * r * x^2 * P)
end

## Coulomb local potential
struct CoulombLocalPotential{S,Analytical} <: AbstractLocalPotential{S,Analytical}
    Z
end

function (Vloc::CoulombLocalPotential{RealSpace})(r::T)::T where {T<:Real}
    return -Vloc.Z / r
end

function (Vloc::CoulombLocalPotential{FourierSpace})(q::T) where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    # General atom => Use default Coulomb potential
    # We use int_{R^3} -Z/r e^{-i q⋅x} = 4π / |q|^2
    return -4T(π) * Vloc.Z / q^2
end

## Coulomb error-function local potential
# Corresponds to a term used by QuantumESPRESSO to perform the Fourier-transform of numeric
# local potentials (i.e. remove and add back the non-local Coulombic part of the potential)
struct CoulombErfLocalPotential{S,Analytical} <: AbstractLocalPotential{S,Analytical}
    Z
end

function (Vloc::CoulombErfLocalPotential{RealSpace})(r)
    return -Vloc.Z / r * erf(r)
end

function (Vloc::CoulombErfLocalPotential{FourierSpace})(q::T) where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    # General atom => Use default Coulomb potential
    # We use int_{R^3} -Z/r erf(r) e^{-i q⋅x} = 4π exp(-|q|^2 / 4) / |q|^2
    return -4T(π) * Vloc.Z / q^2 * exp(-q^2 / 4)
end

## Gaussian local potential
struct GaussianLocalPotential{S,Analytical} <: AbstractLocalPotential{S,Analytical}
    α
    L
end
charge_ionic(Vloc::GaussianLocalPotential) = 0

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
charge_ionic(Vloc::CohenBergstresserLocalPotential) = 4

function (Vloc::CohenBergstresserLocalPotential{FourierSpace})(q::T)::T where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    # Get |q|^2 in units of (2π / lattice_constant)^2
    qsq_pi = Int(round(q^2 / (2π / el.lattice_constant)^2; digits=2))
    return T(get(el.V_sym, qsq_pi, 0.0))
end
