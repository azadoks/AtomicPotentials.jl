struct AugmentationFunction{S,Numerical} <: AbstractAtomicQuantity{S,Numerical}
    r::AbstractVector
    f::AbstractVector
    interpolator::BSplineKit.SplineWrapper
    n::Int
    m::Int
    l::Int
end

struct Augmentation{S<:EvaluationSpace,QT<:AugmentationFunction{S}}
    Q::OffsetVector{Matrix{QT},Vector{Matrix{QT}}}
    q::OffsetVector{Matrix,Vector{Matrix}}
end
function _apply(aug::Augmentation, f::Function, args...; kwargs...)
    Q = map(aug.Q) do Ql
        map(Ql) do Qlnm
            f(Qlnm, args...; kwargs...)
        end
    end
    return Augmentation(Q, aug.q)
end
_apply(::Nothing, ::Function, args...; kwargs...) = nothing
