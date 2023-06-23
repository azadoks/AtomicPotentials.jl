function interpolate_onto(
    quantity::AbstractAtomicQuantity{S,Numerical}, r::AbstractVector
) where {S}
    f = quantity.interpolator.(r)
    constructor = typeof(quantity)
    return _construct_similar_quantity(quantity, constructor; r=r, f=f)
end
