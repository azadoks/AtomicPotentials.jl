function _truncindex(y::AbstractVector{T}; rtol=√eps(T)) where {T<:Real}
    y_max = maximum(y)
    for i in lastindex(y):-1:firstindex(y)
        if abs(y[i] / y_max) > rtol
            return i
        end
    end
    return lastindex(y)
end

function Base.truncate(y::AbstractVector; kwargs...)
    return y[firstindex(y):_truncindex(y, rtol)]
end

function Base.truncate(y::AbstractVector, args...; kwargs...)
    itrunc = firstindex(y):_truncindex(y; kwargs...)
    return y[itrunc], map(Base.Fix2(getindex, itrunc), args)...
end

function Base.truncate(quantity::AbstractAtomicQuantity{S,Numerical}; kwargs...) where {S}
    f, r = Base.truncate(quantity.f, quantity.r; kwargs...)
    return _construct_similar_quantity(quantity; f=f, r=r)
end

Base.truncate(::Nothing, args...; kwargs...) = nothing

function Base.truncate(
    x::Union{NonlocalPotential{S,Numerical},Augmentation{S,Numerical},AtomicPotential},
    args...;
    kwargs...,
) where {S}
    return _apply(x, truncate, args...; kwargs...)
end
