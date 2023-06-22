abstract type AbstractChargeDensity{S,A} <: AbstractAtomicQuantity{S,A} end
angular_momentum(quantity::AbstractChargeDensity)::Int = 0

abstract type AbstractValenceChargeDensity{S,A} <: AbstractChargeDensity{S,A} end

struct ValenceChargeDensity{S,Numeric} <: AbstractValenceChargeDensity{S,Numeric}
    r::AbstractVector
    f::AbstractVector
    interpolator
    l::Int
end

struct GaussianChargeDensity{S,Analytical} <: AbstractChargeDensity{S,Analytical}
    Z
    L
end
function (ρ::GaussianChargeDensity{RealSpace})(r::T)::T where {T<:Real}
    ρ.Z * exp(-(r / (2 * ρ.L))^2) / (sqrt(2) * ρ.L)
end
function (ρ::GaussianChargeDensity{FourierSpace})(q::T)::T where {T<:Real}
    ρ.Z * exp(-(q * ρ.L)^2)
end

abstract type AbstractCoreChargeDensity{S,A} <: AbstractChargeDensity{S,A} end

struct CoreChargeDensity{S,Numeric} <: AbstractCoreChargeDensity{S,Numeric}
    r::AbstractVector
    f::AbstractVector
    interpolator
end
