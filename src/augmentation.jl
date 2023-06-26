struct AugmentationFunction{S,A} <: AbstractAtomicQuantity{S,A}
    r::AbstractVector
    f::AbstractVector
    interpolator
    n::Int
    m::Int
    l::Int
end
angular_momentum(aug::AugmentationFunction) = aug.l

function AugmentationFunction{S}(
    r::AbstractVector, f::AbstractVector, interpolator, n::Integer, m::Integer, l::Integer
) where {S<:EvaluationSpace}
    return AugmentationFunction{S,Numerical}(r, f, interpolator, n, m, l)
end
function fht(
    ::AugmentationFunction{RealSpace,Numerical}, q::AbstractVector, args...; kwargs...
)
    return error("`fht` not implemented for `AugmentationFunction`")
end
function ifht(
    ::AugmentationFunction{RealSpace,Numerical}, r::AbstractVector, args...; kwargs...
)
    return error("`ifht` not implemented for `AugmentationFunction`")
end
