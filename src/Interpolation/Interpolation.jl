module Interpolation
using BSplineKit: BSplineKit
import BSplineKit: Natural

export InterpolationMethod
export Spline
export Natural  # Re-export from BSplineKit
export evaluate

abstract type InterpolationMethod end
function construct_interpolator end

struct Spline <: InterpolationMethod
    order::Int
    bc::Union{Nothing,BSplineKit.Natural}
end
function Spline(order::Integer; bc=BSplineKit.Natural())
    return Spline(order, bc)
end

function construct_interpolator(x::AbstractVector, y::AbstractVector, method::Spline)
    return BSplineKit.interpolate(x, y, BSplineKit.BSplineOrder(method.order), method.bc)
end
function evaluate(i::BSplineKit.SplineWrapper, x::T)::T where {T<:Real}
    return i.spline(x)
end
end
