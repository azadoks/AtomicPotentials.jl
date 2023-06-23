## Interpolation
function interpolate_onto(
    quantity::AbstractAtomicQuantity{S,A}, r::AbstractVector
) where {S,A}
    kwargs = pairs((r=r, f=quantity.(r)))
    fields = OrderedDict(
        name => getproperty(quantity, name) for name in propertynames(quantity)
    )
    for (key, value) in kwargs
        @assert key in keys(fields)
        fields[key] = value
    end
    return typeof(quantity){S,Numeric}(values(fields)...)
end
