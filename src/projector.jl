abstract type AbstractProjector{S,A} <: AbstractAtomicQuantity{S,A} end

### KB projectors
abstract type AbstractKleinmanBylanderProjector{S,A} <: AbstractProjector{S,A} end

# Numerical KB projector
struct KleinmanBylanderProjector{S,Numerical} <:
       AbstractKleinmanBylanderProjector{S,Numerical}
    r::AbstractVector
    f::AbstractVector  # r²β(r) in real-space; β(q) in Fourier-space
    interpolator::BSplineKit.SplineWrapper  # r²β(r) in real-space; β(q) in Fourier-space
    n::Int
    l::Int
    j::Real
end
function KleinmanBylanderProjector{S}(
    r::AbstractVector, f::AbstractVector, interpolator, n::Integer, l::Integer, j::Real
) where {S<:EvaluationSpace}
    return KleinmanBylanderProjector{S,Numerical}(r, f, interpolator, n, l, j)
end

## HGH KB projector
struct HghKleinmanBylanderProjector{S,Analytical} <:
       AbstractKleinmanBylanderProjector{S,Analytical}
    r::Real
    n::Int
    l::Int
end
function hgh_projector_polynomial(
    P::HghKleinmanBylanderProjector{FourierSpace}, t::T
) where {T<:Real}
    c::T = 4T(π)^(5 / T(4)) * sqrt(T(2^(P.l + 1)) * P.r^3)

    # Note: In the (l == 0 && i == 2) case the HGH paper has an error.
    #       The first 8 in equation (8) should not be under the sqrt-sign
    #       This is the right version (as shown in the GTH paper)
    (P.l == 0 && P.n == 1) && return convert(T, c)
    (P.l == 0 && P.n == 2) && return c * 2 / sqrt(T(15)) * (3 - t^2)
    (P.l == 0 && P.n == 3) && return c * 4 / 3sqrt(T(105)) * (15 - 10t^2 + t^4)
    #
    (P.l == 1 && P.n == 1) && return c * 1 / sqrt(T(3)) * t
    (P.l == 1 && P.n == 2) && return c * 2 / sqrt(T(105)) * t * (5 - t^2)
    (P.l == 1 && P.n == 3) && return c * 4 / 3sqrt(T(1155)) * t * (35 - 14t^2 + t^4)
    #
    (P.l == 2 && P.n == 1) && return c * 1 / sqrt(T(15)) * t^2
    (P.l == 2 && P.n == 2) && return c * 2 / 3sqrt(T(105)) * t^2 * (7 - t^2)
    #
    (P.l == 3 && P.n == 1) && return c * 1 / sqrt(T(105)) * t^3

    throw(ArgumentError("Not implemented for l=$l and n=$n"))
end
function (P::HghKleinmanBylanderProjector{RealSpace})(r::T)::T where {T<:Real}
    ired::T = (4 * P.n - 1) / T(2)
    return sqrt(T(2)) * r^(P.l + 2(P.n - 1)) * exp(-r^2 / 2P.r^2) / P.r^(P.l + ired) /
           sqrt(gamma(P.l + ired))
end
function (P::HghKleinmanBylanderProjector{FourierSpace})(q::T) where {T<:Real}
    x::T = q * P.r
    return hgh_projector_polynomial(P, x) * exp(-x^2 / T(2))
end

### State projectors
abstract type AbstractStateProjector{S,A} <: AbstractProjector{S,A} end

## Numerical state projector
struct StateProjector{S,Numerical} <: AbstractStateProjector{S,Numerical}
    r::AbstractVector
    f::AbstractVector  # r²χ(r) in real-space; χ(q) in Fourier-space
    interpolator::BSplineKit.SplineWrapper  # r²χ(r) in real-space; χ(q) in Fourier-space
    n::Int
    l::Int
    j::Real
end

## Hydrogenic state projector
function hydrogenic_projector_radial_1(r::T, α::Real)::T where {T}
    return 2 * α^(3 / 2) * exp(-α * r) * r^2
end
function hydrogenic_projector_radial_2(r::T, α::Real)::T where {T}
    return 2^(-3 / 2) * α^(3 / 2) * (2 - α * r) * exp(-α * r / 2) * r^2
end
function hydrogenic_projector_radial_3(r::T, α::Real)::T where {T}
    return sqrt(4 / 27) *
           α^(3 / 2) *
           (1 - 2 / 3 * α * r + 2 / 27 * α^2 * r^2) *
           exp(-α * r / 3) *
           r^2
end
function hydrogenic_projector_radial_n(n::Integer, α::Real)
    if n == 1
        return Base.Fix2(hydrogenic_projector_radial_1, α)
    elseif n == 2
        return Base.Fix2(hydrogenic_projector_radial_2, α)
    elseif n == 3
        return Base.Fix2(hydrogenic_projector_radial_3, α)
    else
        throw(ArgumentError("n=$(n) is not supported"))
    end
end
struct HydrogenicProjector{S,A} <: AbstractStateProjector{S,A}
    r::Union{Nothing,AbstractVector}
    f::Union{Nothing,AbstractVector}  # r²P(r) in real-space; P(q) in Fourier-space
    interpolator
    n::Int
    l::Int
    α::Real
    function HydrogenicProjector{RealSpace}(r::AbstractVector, n::Int, l::Int, α::Real=1.0)
        interpolator = hydrogenic_projector_radial_n(n, α)
        f = interpolator.(r)
        return new{RealSpace,Numerical}(r, f, interpolator, n, l, α)
    end
    function HydrogenicProjector{FourierSpace,Numerical}(
        q::AbstractVector, F::AbstractVector, interpolator, n::Int, l::Int, α::Real=1.0
    )
        return new{FourierSpace,Numerical}(q, F, interpolator, n, l, α)
    end
end
function iht(proj::HydrogenicProjector, r::AbstractVector, args...; kwargs...)
    return HydrogenicProjector{RealSpace}(r, proj.n, proj.l, proj.α)
end
