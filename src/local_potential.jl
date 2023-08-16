import .NumericalQuadrature: QuadratureMethod, Simpson, integrate

@doc raw"""
Abstract atomic local potential correction.

This correction is used to reduce ringing in the radial Fourier transform of
numerical atomic local potentials. The local potential almost always exhibits a
Coulomb-like tail as it to decays to zero far from the nucleus, as $r\righarrow\infty$.
Therefore, representing the local potential on a radial grid until it decays to zero is
unfeasible. In order to combat this, a Coulomb-like correction which can be analytically
Fourier-transformed is subtracted from the local potential before applying the radial
Fouier transform (rFT). After the rFT, the Fourier transform of the correction is added
back in:

```math
V_{\mathrm{loc}}(q) = 4π \int_{0}^{\infty} r (V_\mathrm{loc}(r) - c(r)) j_l(qr) r dr
+ c(q)
```

See also: [`CoulombLocalPotentialCorrection`](@ref),
[`ErfCoulombLocalPotentialCorrection`](@ref)
"""
abstract type AbstractLocalPotentialCorrection{S<:EvaluationSpace} end

@doc raw"""
Coulombic atomic local potential correction.

In real space, the correction is defined as
```math
c(r) = \frac{-Z}{r},
```
where $Z$ is the (pseudo-)atomic charge. The Fourier transform is
```math
c(q\gt0) = 4π \frac{-Z}{q^2}
```
with a D.C. component
```
c(q=0) = 0.
```
In the radial Fourier transform, one factor of `r` is multiplied through to avoid
divide-by-zero errors:
```math
V_{\mathrm{loc}}(q) = 4\pi \int_{0}^{\infty} (r V_\mathrm{loc}(r) - r c(r)) j_l(qr) r dr
.
```
So in practice, the correction evaluates to
```math
c(r) = -Z
```
in real space.

See also: [`LocalPotentialCorrection`](@ref), [`ErfCoulombLocalPotentialCorrection`](@ref)
"""
struct CoulombLocalPotentialCorrection{S<:EvaluationSpace} <:
       AbstractLocalPotentialCorrection{S} end

function rft(::CoulombLocalPotentialCorrection{RealSpace}, args...; kwargs...)
    return CoulombLocalPotentialCorrection{FourierSpace}()
end

function irft(::CoulombLocalPotentialCorrection{FourierSpace}, args...; kwargs...)
    return CoulombLocalPotentialCorrection{RealSpace}()
end

function (::CoulombLocalPotentialCorrection{RealSpace})(::T, Z)::T where {T}
    return -Z
end

function (::CoulombLocalPotentialCorrection{FourierSpace})(q::T, Z)::T where {T}
    !iszero(q) && return 4T(π) * -Z / q^2
    return zero(T)
end

@doc raw"""
Error-function-modulated Coulombic local potential correction.

In real space, the correction is defined as
```math
c(r) = -\frac{Z}{r} \erf(r),
```
where $Z$ is the (pseudo-)atomic charge and `\erf` is the error function. The
correction in Fourier space is therefore
```math
c(q) = 4\pi \frac{-Z}{q^2} e^{-q^2 / 4}.
```
In the radial Fourier transform, one factor of `r` is multiplied through to avoid
divide-by-zero errors:
```math
V_{\mathrm{loc}}(q) = 4\pi \int_{0}^{\infty} (r V_\mathrm{loc}(r) - r c(r)) j_l(qr) r dr
.
```
So in practice, the correction evaluates to
```math
c(r) = -Z \erf(r)
```
in real space.

See also: [`LocalPotentialCorrection`](@ref), [`CoulombLocalPotentialCorrection`](@ref)
"""
struct ErfCoulombLocalPotentialCorrection{S<:EvaluationSpace} <:
       AbstractLocalPotentialCorrection{S} end

function rft(::ErfCoulombLocalPotentialCorrection{RealSpace}, args...; kwargs...)
    return ErfCoulombLocalPotentialCorrection{FourierSpace}()
end

function irft(::ErfCoulombLocalPotentialCorrection{FourierSpace}, args...; kwargs...)
    return ErfCoulombLocalPotentialCorrection{RealSpace}()
end

function (::ErfCoulombLocalPotentialCorrection{RealSpace})(r::T, Z)::T where {T}
    return -T(Z) * erf(r)
end

function (::ErfCoulombLocalPotentialCorrection{FourierSpace})(q::T, Z)::T where {T}
    !iszero(q) && return 4T(π) * -T(Z) / q^2 * exp(-q^2 / T(4))
    return zero(T)
end

@doc raw"""
Abstract atomic local potential.
"""
abstract type AbstractLocalPotential{S<:EvaluationSpace,A<:Analyticity,T<:Real} <:
              AbstractQuantity{S,A} end

@doc raw"""
Return the (pseudo-)atomic charge generating the atomic local potential.
"""
ionic_charge(quantity::AbstractLocalPotential) = quantity.Z

@doc raw"""
Compute the correction to the energy of the atomic local potential due to the D.C.
component of the Fourier transform
```math
E = 4\pi \int_0^\infty r (V_\mathrm{local}(r) - \frac{-Z}{r}) r dr.
```
"""
function energy_correction(
    quantity::AbstractLocalPotential{RealSpace,<:Analyticity,T}; kwargs...
)::T where {T}
    return T(0)
end

@doc raw"""
Numerical, i.e. non-analytically-evaluatable, atomic local potential.
"""
struct NumericalLocalPotential{S<:EvaluationSpace,Numerical,T<:Real,V<:AbstractVector{T}} <:
       AbstractLocalPotential{S,Numerical,T}
    "Radial grid."
    x::V
    "Radial function."
    f::V
    "Ionic charge."
    Z::T
end

function NumericalLocalPotential{S}(
    x::V, f, Z
) where {S<:EvaluationSpace,T,V<:AbstractVector{T}}
    return NumericalLocalPotential{S,Numerical,T,V}(x, f, Z)
end

function resample(
    quantity::NumericalLocalPotential{S,Numerical,T,V},
    xp::VXP;
    interpolation_method::InterpolationMethod=LinearInterpolation(),
) where {S,T,V,VXP<:AbstractVector}
    TXP = eltype(VXP)
    fp = similar(xp)
    fp .= interpolate(quantity, interpolation_method).(xp)
    return NumericalLocalPotential{S,Numerical,TXP,VXP}(xp, fp, ionic_charge(quantity))
end

angular_momentum(q::NumericalLocalPotential)::Int = 0

n_x_factors(q::NumericalLocalPotential{RealSpace})::Int = 1

function ionic_charge(q::NumericalLocalPotential{S,A,T,V})::T where {S,A,T,V}
    return q.Z
end

function rft(
    quantity::NumericalLocalPotential{RealSpace},
    q::V;
    correction_method::AbstractLocalPotentialCorrection{RealSpace}=CoulombLocalPotentialCorrection{
        RealSpace
    }(),
    kwargs...,
) where {V<:AbstractVector}
    T = eltype(q)
    f =
        radial_grid(quantity) .* radial_function(quantity) .-
        correction_method.(radial_grid(quantity), ionic_charge(quantity))
    F =
        rft(q, radial_grid(quantity), f, 0; n_x_factors=1, kwargs...) .+
        rft(correction_method).(q, T(ionic_charge(quantity)))
    F[iszero.(q)] .= energy_correction(quantity; kwargs...)
    return NumericalLocalPotential{FourierSpace,Numerical,T,V}(
        q, F, T(ionic_charge(quantity))
    )
end

function irft(
    quantity::NumericalLocalPotential{FourierSpace},
    r::V;
    correction_method::AbstractLocalPotentialCorrection{FourierSpace}=CoulombLocalPotentialCorrection{
        FourierSpace
    }(),
    kwargs...,
) where {V<:AbstractVector}
    T = eltype(r)
    F =
        radial_grid(quantity) .* radial_function(quantity) .-
        correction_method.(radial_grid(quantity), ionic_charge(quantity))
    f =
        1 / (2T(π))^3 .* (
            irft(r, radial_grid(quantity), F, 0; n_x_factors=1, kwargs...) .+
            ift(correction_method).(r, T(ionic_charge(quantity)))
        )
    f[iszero.(r)] .= T(-Inf)
    return NumericalLocalPotential{RealSpace,Numerical,T,V}(r, f, T(ionic_charge(quantity)))
end

function energy_correction(
    quantity::NumericalLocalPotential{RealSpace,Numerical,T};
    quadrature_method::QuadratureMethod=Simpson(),
    kwargs...,
)::T where {T}
    integrand =
        radial_grid(quantity) .*
        (radial_grid(quantity) .* radial_function(quantity) .- -ionic_charge(quantity))
    return 4T(π) * integrate(radial_grid(quantity), integrand, quadrature_method)
end

@doc raw"""
Hartwigsen-Goedecker-Hutter atomic local potential.

The local potential is defined in real space in Eq. 1 of Phys. Rev. B 54 1703:
```math
V_\mathrm{loc}(r) - \frac{-Z}{r} \erf \left( \frac{r}{\sqrt{2} r_\mathrm{loc}} \right)
+ \exp \left[ -\frac{1}{2} \left( \frac{r}{r_\mathrm{loc}} \right)^2 \right] \cross
\left[
c_1 +
c_2 \left( frac{r}{r_\mathrm{loc}} \right)^2 +
c_3 \left( frac{r}{r_\mathrm{loc}} \right)^4 +
c_4 \left( frac{r}{r_\mathrm{loc}} \right)^6
\right].

The local potential is defined in Fourier space in Eq. 6 of Phys. Rev. B 54 1703 in
terms of plane waves normalized by $\frac{1}{\sqrt{\Omega}}$.
It can be brought to the form
```
V_\mathrm{local}(t) = \frac{P(t)}{t^2 \exp \left( t^2 / 2 \right)}
```
where $t = q r_\mathrm{local}$ and $P$ is a polynomial of at most degree 8.

[Phys. Rev. B 54, 1703 (1996)](https://doi.org/10.1103/PhysRevB.54.1703)
[Phys. Rev. B 58, 3651 (1998)](https://doi.org/10.1103/PhysRevB.58.3641)
"""
struct HghLocalPotential{S,Analytical,T<:Real} <: AbstractLocalPotential{S,Analytical,T}
    "Range of the Gaussian ionic charge."
    rₗ::T
    raw"Polynomial coefficients defining $P(\frac{r}{r_{\mathrm{loc}})$."
    cₗ::Union{NTuple{1,T},NTuple{2,T},NTuple{3,T},NTuple{4,T}}
    "Ionic charge."
    Z::T
end
function HghLocalPotential{S}(rₗ::T, cₗ, Z) where {S<:EvaluationSpace,T}
    return HghLocalPotential{S,Analytical,T}(rₗ, cₗ, Z)
end

@doc raw"""
Evaluate the HGH local potential in the space defined by its type parameters.
"""
function (quantity::HghLocalPotential{RealSpace})(r::T)::T where {T<:Real}
    # TODO: just return Inf for r=0?
    r += iszero(r) ? eps(T) : zero(T)  # quick hack for the division by zero below
    rr = r / T(quantity.rₗ)
    cₗ = convert.(T, quantity.cₗ)
    # cₗ contains a variable number of coefficients ranging from 1 to 4.
    # Each entry is the polynomial coefficient of order 2(i-1), where i is the coefficient's
    # index in cₗ.
    # Here, we interleave cₗ with zeros
    #   [c₁, c₂] -> [c₁, 0, c₂, 0]
    # and then drop the last element
    #   [c₁, 0, c₂, 0] -> [c₁, 0, c₂],
    # which yields the correct array of polynomial coefficients for `evalpoly`.
    c = collect(Iterators.flatten(zip(cₗ, zeros(T, length(cₗ)))))[begin:(end - 1)]
    return -ionic_charge(quantity) / r * erf(rr / sqrt(T(2))) +
           exp(-rr^2 / 2) * evalpoly(rr, c)
end

function (quantity::HghLocalPotential{FourierSpace})(q::T)::T where {T<:Real}
    iszero(q) && return energy_correction(quantity)
    x = q * T(quantity.rₗ)
    # The polynomial prefactor P(t) (as used inside the { ... } brackets of equation
    # (5) of the HGH98 paper)
    poly_coeff = ((1,), (3, 0, -1), (15, 0, -10, 0, 1), (105, 0, -105, 0, 21, 0, -1))
    P = sum(zip(quantity.cₗ, poly_coeff)) do (cₗ_i, c_i)
        return T(cₗ_i) * evalpoly(x, c_i)
    end
    return 4T(π) *
           quantity.rₗ^2 *
           (-ionic_charge(quantity) + sqrt(T(π) / 2) * quantity.rₗ * x^2 * P) *
           exp(-x^2 / 2) / x^2
end

function energy_correction(
    quantity::HghLocalPotential{<:EvaluationSpace,Analytical,T}; kwargs...
) where {T}
    coeffs = (1, 3, 15, 105)
    return 4T(π) * (
        ionic_charge(quantity) * quantity.rₗ^2 / 2 +
        sqrt(π / T(2)) * quantity.rₗ^3 * sum(coeffs .* quantity.cₗ)
    )
end

function rft(quantity::HghLocalPotential{RealSpace,A,T}, args...; kwargs...) where {A,T}
    return HghLocalPotential{FourierSpace,A,T}(
        quantity.rₗ, quantity.cₗ, ionic_charge(quantity)
    )
end

function irft(quantity::HghLocalPotential{FourierSpace,A,T}, args...; kwargs...) where {A,T}
    return HghLocalPotential{RealSpace,A,T}(
        quantity.rₗ, quantity.cₗ, ionic_charge(quantity)
    )
end

@doc raw"""
Coulombic atomic local potential.

```math
V_\mathrm{local}(r) = \frac{-Z}{r}.
```

```math
V_\mathrm{local}(q) =
\int_{\mathbb{R^3}} \frac{Z}{\mathbf{r}} e^{-i \mathbf{q} \cdot \mathbf{r}} =
4\pi \frac{-Z}{|q|^2}.
```
"""
struct CoulombLocalPotential{S,Analytical,T} <: AbstractLocalPotential{S,Analytical,T}
    "Ionic charge."
    Z::T
end

@doc raw"""
Evaluate the Coulomb local potential in the space defined by its type parameters.
"""
function (quantity::CoulombLocalPotential{RealSpace})(r::T)::T where {T}
    return -ionic_charge(quantity) / r
end

function (quantity::CoulombLocalPotential{FourierSpace})(q::T)::T where {T}
    iszero(q) && return zero(T)
    return 4T(π) * -T(ionic_charge(quantity)) / q^2
end

function rft(quantity::CoulombLocalPotential{RealSpace,A,T}, args...; kwargs...) where {A,T}
    return CoulombLocalPotential{FourierSpace,A,T}(ionic_charge(quantity))
end

function irft(
    quantity::CoulombLocalPotential{FourierSpace,A,T}, args...; kwargs...
) where {A,T}
    return CoulombLocalPotential{RealSpace,A,T}(ionic_charge(quantity))
end

@doc raw"""
Error function modulated Coulombic atomic local potential.

```math
V_\mathrm{local}(r) = \frac{-Z}{r} \erf(r).
```

```math
V_\mathrm{local}(q) =
\int_{\mathbb{R^3}} \frac{Z}{\mathbf{r}} \erf(r) e^{-i \mathbf{q} \cdot \mathbf{r}} =
4\pi \frac{-Z}{|q|^2} \exp\left( \frac{-|q|^2}{4} \right).
```
"""
struct CoulombErfLocalPotential{S,Analytical,T} <: AbstractLocalPotential{S,Analytical,T}
    "Ionic charge."
    Z::T
end

@doc raw"""
Evaluate the error function modulated Coulomb local potential in the space defined by its
type parameters.
"""
function (quantity::CoulombErfLocalPotential{RealSpace})(r::T)::T where {T}
    return -ionic_charge(quantity) / r * erf(r)
end

function (quantity::CoulombErfLocalPotential{FourierSpace})(q::T)::T where {T}
    iszero(q) && return zero(T)
    return 4T(π) * -T(ionic_charge(quantity)) / q^2 * exp(-q^2 / T(4))
end

function rft(quantity::CoulombErfLocalPotential{RealSpace,A,T}) where {A,T}
    return CoulombErfLocalPotential{FourierSpace,A,T}(ionic_charge(quantity))
end

function irft(
    quantity::CoulombErfLocalPotential{FourierSpace,A,T}, args...; kwargs...
) where {A,T}
    return CoulombErfLocalPotential{RealSpace,A,T}(ionic_charge(quantity))
end

@doc raw"""
Gaussian local potential.

```math
V_\mathrm{local}(r) = \frac{\alpha}{L \sqrt{2\pi}}
\exp \left(-\frac{1}{2} \frac{r^2}{L^2} \right)
```

```math
V_\mathrm{local}(q) = -\alpha \exp \left(-\frac{(q L)^2}{2} \right)
```
"""
struct GaussianLocalPotential{S,Analytical,T<:Real} <:
       AbstractLocalPotential{S,Analytical,T}
    "Scale factor."
    α::T
    "Gaussian width."
    L::T
end

ionic_charge(quantity::GaussianLocalPotential) = 0

@doc raw"""
Evaluate the Gaussian local potential in the space defined by its type parameters.
"""
function (quantity::GaussianLocalPotential{RealSpace})(r::T)::T where {T}
    return -T(quantity.α) / (sqrt(2T(π)) * T(quantity.L)) *
           exp(-(r / T(quantity.L))^2 / T(2))
end

function (quantity::GaussianLocalPotential{FourierSpace})(q::T)::T where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    return -T(quantity.α) * exp(-(q * T(quantity.L))^2 / T(2))
end

function rft(
    quantity::GaussianLocalPotential{RealSpace,A,T}, args...; kwargs...
) where {A,T}
    return GaussianLocalPotential{FourierSpace,A,T}(quantity.α, quantity.L)
end

function irft(
    quantity::GaussianLocalPotential{FourierSpace,A,T}, args...; kwargs...
) where {A,T}
    return GaussianLocalPotential{RealSpace,A,T}(quantity.α, quantity.L)
end

# TODO: document
struct CohenBergstresserLocalPotential{FourierSpace,Analytical,T<:Real} <:
       AbstractLocalPotential{FourierSpace,Analytical,T}
    V_sym::Dict{Int,T}  # Map |G|^2 (in units of (2π / lattice_constant)^2) to form factors
    a::T  # Lattice constant (in Bohr) which is assumed
end

ionic_charge(quantity::CohenBergstresserLocalPotential) = 4

function (quantity::CohenBergstresserLocalPotential{FourierSpace})(q::T)::T where {T}
    iszero(q) && return zero(T)  # Compensating charge background
    # Get |q|^2 in units of (2π / lattice_constant)^2
    qsq_pi = Int(round(q^2 / (2T(π) / quantity.a)^2; digits=2))
    return T(get(el.V_sym, qsq_pi, 0.0))
end
