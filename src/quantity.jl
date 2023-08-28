import .NumericalQuadrature: QuadratureMethod, Simpson, integration_weights
import .Interpolation: InterpolationMethod, Spline, interpolate

@doc raw"""
Abstract evaluation space.
"""
abstract type EvaluationSpace end
@doc raw"""
Singleton representing real/direct space.
"""
struct RealSpace <: EvaluationSpace end
@doc raw"""
Singleton representing Fourier/reciprocal space.
"""
struct FourierSpace <: EvaluationSpace end

@doc raw"""
Abstract analyticity.
"""
abstract type Analyticity end
@doc raw"""
Singleton representing an analytically-evaluatable atomic quantity.
"""
struct Analytical <: Analyticity end
@doc raw"""
Singleton representing a numerical, i.e. non-analytically-evaluatable, atomic quantity.
"""
struct Numerical <: Analyticity end

@doc raw"""
Abstract atomic potential quantity.
"""
abstract type AbstractQuantity{S<:EvaluationSpace,A<:Analyticity} end

Base.broadcastable(qty::AbstractQuantity) = Ref(qty)
Base.isempty(qty::AbstractQuantity) = false

@doc raw"""
Return the angular momentum of an atomic potential quantity.
"""
angular_momentum(quantity::AbstractQuantity)::Int = quantity.l

@doc raw"""
Return the radial grid of a numerical atomic potential quantity.
"""
radial_grid(quantity::AbstractQuantity{<:EvaluationSpace,Numerical}) = quantity.x

@doc raw"""
Return the radial function of a numerical atomic potential quantity.
"""
radial_function(quantity::AbstractQuantity{<:EvaluationSpace,Numerical}) = quantity.f

@doc raw"""
Return the number of times the radial grid is pre-multiplied into the radial function
of an atomic potential quantity.
"""
n_x_factors(quantity::AbstractQuantity{<:EvaluationSpace,Numerical})::Int =
    quantity.n_x_factors
n_x_factors(quantity::AbstractQuantity{<:EvaluationSpace,Analytical})::Int = 0

@doc raw"""
Construct an interpolation of the radial function of an atomic potential quantity.
"""
function Interpolation.interpolate(
    quantity::AbstractQuantity{<:EvaluationSpace,Numerical},
    interpolation_method::InterpolationMethod=Spline(4),
    args...;
    kwargs...,
)
    return interpolate(
        radial_grid(quantity),
        radial_function(quantity),
        interpolation_method,
        args...;
        kwargs...,
    )
end

@doc raw"""
Return a callable which evaluates the radial function of an atomic potential quantity.
"""
function evaluate(
    quantity::AbstractQuantity{<:EvaluationSpace,Numerical}, args...; kwargs...
)
    return interpolate(quantity, args..., kwargs...)
end
function evaluate(
    quantity::AbstractQuantity{<:EvaluationSpace,Analytical}, args...; kwargs...
)
    return quantity
end

@doc raw"""
Resample a numeric atomic potential quantity onto a uniform radial grid with spacing
at most dx.
"""
function resample(
    quantity::AbstractQuantity{<:EvaluationSpace,Numerical}, dxp::Real; kwargs...
)
    x = radial_grid(q)
    x_min, x_max = extrema(x)
    n_xp = ceil(Int, (x_max - x_min) / dx)
    xp = range(x_min, x_max, n_xp)
    return resample(quantity, xp; kwargs...)
end

function resample(
    quantity::AbstractQuantity{<:EvaluationSpace,Numerical}, n_xp::Integer; kwargs...
)
    x = radial_grid(q)
    x_min, x_max = extrema(x)
    xp = range(x_min, x_max, n_xp)
    return resample(quantity, xp; kwargs...)
end

## Symmetric array support

LinearAlgebra.symmetric_type(::Type{T}) where {T<:AbstractQuantity} = T
LinearAlgebra.symmetric(quantity::AbstractQuantity, ::Symbol) = quantity
LinearAlgebra.transpose(quantity::AbstractQuantity) = quantity

@doc raw"""
Numerical, i.e. non-analytically-evaluatable, atomic quantity.

See also: [`AbstractQuantity`](@ref), [`AtomicLocalPotential`](@ref)
"""
struct NumericalQuantity{S<:EvaluationSpace,Numerical,T<:Real,V<:AbstractVector{T}} <:
       AbstractQuantity{S,Numerical}
    "Radial grid."
    x::V
    "Radial function."
    f::V
    "Angular momentum."
    l::Int
    "Number of times the radial grid is pre-multiplied into the radial function."
    n_x_factors::Int
end

function NumericalQuantity{S}(
    x::V, f::V, l::Integer, n_x_factors::Integer=0
) where {S<:EvaluationSpace,T<:Real,V<:AbstractVector{T}}
    return NumericalQuantity{S,Numerical,T,V}(x, f, l, n_x_factors)
end

function Base.zero(::Type{<:NumericalQuantity{S,A,T,V}}) where {S,A,T,V}
    return NumericalQuantity{S,A,T,V}(V(), V(), -1, 0)
end

function Base.iszero(quantity::NumericalQuantity)
    return isempty(quantity.x) || isempty(quantity.f) || quantity.l == -1
end

function Base.show(io::IO, quantity::NumericalQuantity{RealSpace,Numerical,T}) where {T}
    return print(io, "f{$(T)}[$(length(quantity.x)),l=$(quantity.l)](r)")
end

function Base.show(io::IO, quantity::NumericalQuantity{FourierSpace,Numerical,T}) where {T}
    return print(io, "F{$(T)}[$(length(quantity.x)),l=$(quantity.l)](q)")
end

function Base.convert(::Type{T}, x::NumericalQuantity{S}) where {T<:Real,S<:EvaluationSpace}
    return NumericalQuantity{S}(convert.(T, x.x), convert.(T, x.f), x.l, x.n_x_factors)
end

function Adapt.adapt_structure(to, x::NumericalQuantity{S}) where {S<:EvaluationSpace}
    return NumericalQuantity{S}(adapt(to, x.x), adapt(to, x.f), x.l, x.n_x_factors)
end

@doc raw"""
Re-sample a numerical atomic quantity's radial function on a new set of radial points
`xp` by interpolating using `interpolation_method`.

See also: [`InterpolationMethod`](@ref), [`Linear`](@ref), [`Spline`](@ref),
"""
function resample(
    quantity::NumericalQuantity{S,Numerical,T,V},
    xp::VXP;
    interpolation_method::InterpolationMethod=Spline(4),
) where {S,T,V,TXP<:Real,VXP<:AbstractVector{TXP}}
    itp = interpolate(quantity, interpolation_method)
    fp = map(xp) do xpi
        itp(xpi)
    end
    return NumericalQuantity{S,Numerical,TXP,VXP}(
        xp, fp, angular_momentum(quantity), n_x_factors(quantity)
    )
end

function rft(
    quantity::NumericalQuantity{RealSpace}, q::V; kwargs...
) where {V<:AbstractVector}
    T = eltype(q)
    iszero(quantity) && return zero(NumericalQuantity{FourierSpace,Numerical,T,V})
    F::V = rft(
        q,
        radial_grid(quantity),
        radial_function(quantity),
        angular_momentum(quantity);
        n_x_factors=n_x_factors(quantity),
        kwargs...,
    )
    return NumericalQuantity{FourierSpace,Numerical,T,V}(
        q, F, angular_momentum(quantity), 0
    )
end

function irft(
    quantity::NumericalQuantity{FourierSpace}, r::V; kwargs...
) where {V<:AbstractVector}
    T = eltype(r)
    iszero(quantity) && return zero(NumericalQuantity{RealSpace,Numerical,T,V})
    f::V = irft(
        r,
        radial_grid(quantity),
        radial_function(quantity),
        angular_momentum(quantity);
        n_x_factors=n_x_factors(quantity),
        kwargs...,
    )
    return NumericalQuantity{FourierSpace,Numerical,T,V}(
        r, f, angular_momentum(quantity), 0
    )
end

@doc raw"""
Hartwigsen-Goedecker-Hutter non-local potential "β"-projectors.

The non-local projectors are defined in Real space in Eq. 3 of Phys. Rev. B 54 1703 as:
```math
\beta_n^l(r) = \frac{
    \sqrt{2} r^{l + 2(n-1)} \exp \left( -\frac{r^2}{2r_l^2} \right)
}{
    r_l^{l + (4n-1)/2} \sqrt{\Gamma\left(l + \frac{4n-1}{2} \right)}
},
```
where $\Gamma$ is the gamma function.

The non-local projectors are defined in Fourier space in Eq.s 5-15 of Phys. Rev. B 54
1703, although there is a minor error: the first 8 in Eq. 8 should not be under the
square-root. This can be confirmed by cross-reference with Pys. Rev. B 58 3651.

[Phys. Rev. B 54, 1703 (1996)](https://doi.org/10.1103/PhysRevB.54.1703)
[Phys. Rev. B 58, 3651 (1998)](https://doi.org/10.1103/PhysRevB.58.3641)
"""
struct HghProjector{S<:EvaluationSpace,Analytical,T<:Real} <: AbstractQuantity{S,Analytical}
    ""
    "Gaussian width."
    rₚ::T
    "Index."
    n::Int
    "Angular momentum."
    l::Int
end
function HghProjector{S}(rₚ::T, n, l) where {S<:EvaluationSpace,T}
    return HghProjector{S,Analytical,T}(rₚ, n, l)
end

function Base.convert(::Type{T}, x::HghProjector{S}) where {T<:Real,S<:EvaluationSpace}
    return HghProjector{S}(convert.(T, x.rₚ), x.n, x.l)
end

@doc raw"""
Evaluate an HGH non-local projector in the space defined by its type parameters.
"""
function (quantity::HghProjector{RealSpace})(r::T)::T where {T<:Real}
    n = quantity.n
    l = quantity.l
    rₚ = quantity.rₚ
    num = sqrt(T(2)) * r^(l + 2(n - 1)) * exp(-r^2 / 2rₚ^2)
    denom = rₚ^(l + (4n - 1) / T(2)) * sqrt(gamma(l + (4n - 1) / T(2)))
    return num / denom
end

function (quantity::HghProjector{FourierSpace})(q::T)::T where {T<:Real}
    n = quantity.n
    l = quantity.l
    x = q * T(quantity.rₚ)
    c = 4T(π)^(5 / T(4)) * sqrt(T(2^(l + 1)) * quantity.rₚ^3)
    e = exp(-x^2 / T(2))
    ec = T(e) * T(c)
    if n == 1
        (l == 0) && return ec
        (l == 1) && return ec / sqrt(T(3)) * evalpoly(x, (0, 1))
        (l == 2) && return ec / sqrt(T(15)) * evalpoly(x, (0, 0, 1))
        (l == 3) && return ec / sqrt(T(105)) * evalpoly(x, (0, 0, 0, 1))
    elseif n == 2
        (l == 0) && return ec * 2 / sqrt(T(15)) * evalpoly(x, (3, 0, -1))
        (l == 1) && return ec * 2 / sqrt(T(105)) * evalpoly(x, (0, 5, 0, -1))
        (l == 2) && return ec * 2 / (3 * sqrt(T(105))) * evalpoly(x, (0, 0, 7, 0, -1))
    elseif n == 3
        (l == 0) && return ec * 4 / (3 * sqrt(T(105))) * evalpoly(x, (15, 0, -10, 0, 1))
        (l == 1) && return ec * 4 / (3 * sqrt(T(1155))) * evalpoly(x, (0, 35, 0, -14, 0, 1))
    else
        throw(ArgumentError("Not implemented for l=$l and n=$n"))
    end
end

function rft(quantity::HghProjector{RealSpace,A,T}, args...; kwargs...) where {A,T}
    return HghProjector{FourierSpace,A,T}(
        quantity.rₚ, quantity.n, angular_momentum(quantity)
    )
end

function irft(quantity::HghProjector{FourierSpace,A,T}, args...; kwargs...) where {A,T}
    return HghProjector{RealSpace,A,T}(quantity.rₚ, quantity.n, angular_momentum(quantity))
end

@doc raw"""
Hydrogenic orbital.

!!! note
    Hydrogenic orbitals do not have an analytical Fourier transform. In order to adhere
    to the interface, they must be constructed with a radial mesh in real space.
    This mesh is only used to compute the radial Fourier transform; evaluation in
    real space is done analytically. Likewise, the inverse radial Fourier transform
    simply reconstructs a `RealSpace`-parameterized structure and does not perform
    the numerical inverse radial Fourier transform.
"""
struct HydrogenicProjector{
    S<:EvaluationSpace,A<:Analyticity,T<:Real,V<:AbstractVector{T}
} <: AbstractQuantity{S,A}
    "Radial grid."
    x::V
    "Radial function."
    f::V
    "Angular momentum."
    l::Int
    "Principal quantum number."
    n::Int
    "Orbital radius."
    α::T
    function HydrogenicProjector{RealSpace}(
        r::V, l::Integer, n::Integer, α::T
    ) where {T<:Real,V<:AbstractVector{T}}
        f = similar(r)
        f .= hydrogenic_projector_radial.(r, n, α)
        return new{RealSpace,Analytical,T,V}(r, f, l, n, α)
    end
    function HydrogenicProjector{FourierSpace}(
        r::V, f::V, l::Integer, n::Integer, α::T
    ) where {T<:Real,V<:AbstractVector{T}}
        return new{FourierSpace,Numerical,T,V}(r, f, l, n, α)
    end
end

function Base.convert(::Type{T}, x::HydrogenicProjector{RealSpace}) where {T<:Real}
    return HydrogenicProjector{RealSpace}(convert.(T, x.x), x.l, x.n, convert(T, x.α))
end

function Base.convert(::Type{T}, x::HydrogenicProjector{FourierSpace}) where {T<:Real}
    return HydrogenicProjector{FourierSpace}(
        convert.(T, x.x), convert.(T, x.f), x.l, x.n, convert(T, x.α)
    )
end

function Adapt.adapt_structure(to, x::HydrogenicProjector{RealSpace})
    return HydrogenicProjector{RealSpace}(adapt(to, x.x), x.l, x.n, x.α)
end

function Adapt.adapt_structure(to, x::HydrogenicProjector{FourierSpace})
    return HydrogenicProjector{FourierSpace}(adapt(to, x.x), adapt(to, x.f), x.l, x.n, x.α)
end

# TODO: cite the source of the expressions!
@doc raw"""
Evaluate a hydrogenic orbital in real space analytically.
"""
function (quantity::HydrogenicProjector{RealSpace,Analytical})(r::T)::T where {T<:Real}
    return hydrogenic_projector_radial(r, quantity.n, quantity.α)
end

@doc raw"""
Evaluate a hydrogenic orbital with principal quantum number `n=1` and radius `α`.

```math
R_1(r) = 2 α^(3/2) e^(-\alpha r)
```
"""
function hydrogenic_projector_radial_1(r::T, α::Real)::T where {T}
    return 2 * α^(3 / 2) * exp(-α * r)
end

@doc raw"""
Evaluate a hydrogenic orbital with principal quantum number `n=2` and radius `α`.

```math
R_2(r) = 2^(-3/2) α^(3/2) (2 - \alpha r) e^(-\alpha r / 2)
```
"""
function hydrogenic_projector_radial_2(r::T, α::Real)::T where {T}
    return 2^(-3 / 2) * α^(3 / 2) * (2 - α * r) * exp(-α * r / 2)
end

@doc raw"""
Evaluate a hydrogenic orbital with principal quantum number `n=3` and radius `α`.

```math
R_3(r) =
\sqrt{\frac{4}{27}} α^(3/2)
(1 - \frac{2}{3} \alpha r + \frac{2}{27} (\alpha r)^2)
e^(-\alpha r / 3)
```
"""
function hydrogenic_projector_radial_3(r::T, α::Real)::T where {T}
    return sqrt(4 / 27) *
           α^(3 / 2) *
           (1 - 2 / 3 * α * r + 2 / 27 * α^2 * r^2) *
           exp(-α * r / 3)
end

@doc raw"""
Return a callable which evaluates a hydrogenic orbital in real space.
"""
function hydrogenic_projector_radial_n(n::Integer, α::Real)
    (n == 1) && return Base.Fix2(hydrogenic_projector_radial_1, α)
    (n == 2) && return Base.Fix2(hydrogenic_projector_radial_2, α)
    (n == 3) && return Base.Fix2(hydrogenic_projector_radial_3, α)
    throw(ArgumentError("n=$(n) is not supported"))
end

@doc raw"""
Evaluate a hydrogenic orbital in real space.
"""
function hydrogenic_projector_radial(r::T, n::Integer, α::Real)::T where {T<:Real}
    (n == 1) && return hydrogenic_projector_radial_1(r, α)
    (n == 2) && return hydrogenic_projector_radial_2(r, α)
    (n == 3) && return hydrogenic_projector_radial_3(r, α)
    throw(ArgumentError("n=$(n) is not supported"))
end

function resample(
    quantity::HydrogenicProjector{RealSpace,Analytical,T,V}, xp::VXP; kwargs...
) where {T,V,VXP<:AbstractVector}
    TXP = eltype(VXP)
    return HydrogenicProjector{FourierSpace,Analytical,TXP,VXP}(
        xp, angular_momentum(quantity), quantiy.n, quantity.α
    )
end

function resample(
    quantity::HydrogenicProjector{FourierSpace,Numerical,T,V},
    xp::VXP;
    interpolation_method::InterpolationMethod=Spline(4),
) where {T,V,VXP<:AbstractVector}
    TXP = eltype(VXP)
    fp = similar(xp)
    fp .= interpolation_method(quantity).(xp)
    return HydrogenicProjector{FourierSpace,Numerical,TXP,VXP}(
        xp, fp, angular_momentum(quantity), quantiy.n, quantity.α
    )
end

n_x_factors(::HydrogenicProjector)::Int = 0
function rft(quantity::HydrogenicProjector{RealSpace}, q::V; kwargs...) where {V}
    T = eltype(q)
    F::V = rft(
        q,
        radial_grid(quantity),
        radial_function(quantity),
        angular_momentum(quantity);
        kwargs...,
    )
    return HydrogenicProjector{FourierSpace,Numerical,T,V}(
        q, F, angular_momentum(quantity), quantity.n, quantity.α
    )
end

function irft(quantity::HydrogenicProjector{FourierSpace}, r::V; kwargs...) where {V}
    T = eltype(r)
    return HydrogenicProjector{RealSpace,Analytical,T,V}(
        r, angular_momentum(quantity), quantity.n, quantity.α
    )
end

struct GaussianChargeDensity{S<:EvaluationSpace,Analytical,T<:Real} <:
       AbstractQuantity{S,Analytical}
    Z::T
    L::T
end
function GaussianChargeDensity{S}(Z::Real, L::T) where {S<:EvaluationSpace,T<:Real}
    return GaussianChargeDensity{S,Analytical,T}(Z, L)
end
function GaussianChargeDensity{S}(
    n_elec_core::Integer, n_elec_valence::Integer
) where {S<:EvaluationSpace}
    L = gaussian_density_decay_length(n_elec_core, n_elec_valence)
    return GaussianChargeDensity{S}(n_elec_valence, L)
end

function (qty::GaussianChargeDensity{RealSpace})(r::T)::T where {T<:Real}
    return T(qty.Z) * exp(-(r / (T(2) * T(qty.L)))^2) / (sqrt(T(2)) * T(qty.L))
end

function (qty::GaussianChargeDensity{FourierSpace})(q::T)::T where {T<:Real}
    return T(qty.Z) * exp(-(q * T(qty.L))^2)
end

function rft(qty::GaussianChargeDensity{RealSpace}, args...; kwargs...)
    return GaussianChargeDensity{FourierSpace}(qty.Z, qty.L)
end

function irft(qty::GaussianChargeDensity{FourierSpace}, args...; kwargs...)
    return GaussianChargeDensity{RealSpace}(qty.Z, qty.L)
end

# Get the lengthscale of the valence density for an atom with `n_elec_core` core
# and `n_elec_valence` valence electrons.
function gaussian_density_decay_length(n_elec_core, n_elec_valence)
    # Adapted from ABINIT/src/32_util/m_atomdata.F90,
    # from which also the data has been taken.

    n_elec_valence = round(Int, n_elec_valence)
    if n_elec_valence == 0
        return 0.0
    end

    data = if n_elec_core < 0.5
        # Bare ions: Adjusted on 1H and 2He only
        [0.6, 0.4, 0.3, 0.25, 0.2]
    elseif n_elec_core < 2.5
        # 1s2 core: Adjusted on 3Li, 6C, 7N, and 8O
        [1.8, 1.4, 1.0, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3]
    elseif n_elec_core < 10.5
        # Ne core (1s2 2s2 2p6): Adjusted on 11na, 13al, 14si and 17cl
        [2.0, 1.6, 1.25, 1.1, 1.0, 0.9, 0.8, 0.7, 0.7, 0.7, 0.6]
    elseif n_elec_core < 12.5
        # Mg core (1s2 2s2 2p6 3s2): Adjusted on 19k, and on n_elec_core==10
        [1.9, 1.5, 1.15, 1.0, 0.9, 0.8, 0.7, 0.6, 0.6, 0.6, 0.5]
    elseif n_elec_core < 18.5
        # Ar core (Ne + 3s2 3p6): Adjusted on 20ca, 25mn and 30zn
        [2.0, 1.8, 1.5, 1.2, 1.0, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.65, 0.6]
    elseif n_elec_core < 28.5
        # Full 3rd shell core (Ar + 3d10): Adjusted on 31ga, 34se and 38sr
        [1.5, 1.25, 1.15, 1.05, 1.00, 0.95, 0.95, 0.9, 0.9, 0.85, 0.85, 0.80, 0.8, 0.75, 0.7]
    elseif n_elec_core < 36.5
        # Krypton core (Ar + 3d10 4s2 4p6): Adjusted on 39y, 42mo and 48cd
        [2.0, 2.00, 1.60, 1.40, 1.25, 1.10, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.7]
    else
        # For the remaining elements, consider a function of n_elec_valence only
        [2.0, 2.00, 1.55, 1.25, 1.15, 1.10, 1.05, 1.0, 0.95, 0.9, 0.85, 0.85, 0.8]
    end
    return data[min(n_elec_valence, length(data))]
end
