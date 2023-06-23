module Interpolation
using BSplineKit: BSplineKit

export InterpolationMethod
export Spline

abstract type InterpolationMethod end
function construct_interpolator end

struct Spline <: InterpolationMethod
    order::Int
end
function construct_interpolator(x::AbstractVector, y::AbstractVector, method::Spline)
    return BSplineKit.interpolate(
        x, y, BSplineKit.BSplineOrder(method.order), BSplineKit.Natural()
    )
end
end
