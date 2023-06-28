function _truncindex(y::AbstractVector, rtol)
    y_max = maximum(y)
    for i in lastindex(y):-1:firstindex(y)
        if abs(y[i] / y_max) > rtol
            return i
        end
    end
    return lastindex(y)
end

function Base.truncate(y::AbstractVector{T}; rtol=√eps(T)) where {T<:Real}
    return y[firstindex(y):_truncindex(y, rtol)]
end

function Base.truncate(y::AbstractVector{T}, args...; rtol=√eps(T)) where {T<:Real}
    itrunc = firstindex(y):_truncindex(y, rtol)
    return y[itrunc], map(Base.Fix2(getindex, itrunc), args)...
end

function Base.truncate(quantity::AbstractAtomicQuantity{S,Numerical}; kwargs...) where {S}
    f, r = Base.truncate(quantity.f, quantity.r; kwargs...)
    return _construct_similar_quantity(quantity; f=f, r=r)
end

function Base.truncate(
    x::Union{NonLocalPotential{S,Numerical},Augmentation{S,Numerical},AtomicPotential},
    args...;
    kwargs...,
) where {S}
    return _apply(x, truncate, args...; kwargs...)
end
