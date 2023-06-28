interpolate_onto(::Nothing, args...; kwargs...) = nothing

function interpolate_onto(
    quantity::AbstractAtomicQuantity{S,Numerical}, r::AbstractVector
) where {S}
    f = quantity.interpolator.(r)
    constructor = typeof(quantity)
    return _construct_similar_quantity(quantity, constructor; r=r, f=f)
end

function interpolate_onto(quantity::AbstractAtomicQuantity{S,Numerical}, Δr::Real) where {S}
    r_min, r_max = extrema(quantity.r)
    # Ensure that the endpoints of the radial mesh stay in the new radial mesh
    # by massaging Δr to be _at most_ the value provided, i.e. reduce the spacing
    # Δr until a range with an integer number of intervals containing the endpoints
    # can be made.
    r = range(; start=r_min, stop=r_max, length=ceil(Int, (r_max - r_min) / Δr))
    return interpolate_onto(quantity, r)
end

function interpolate_onto(
    x::Union{NonLocalPotential{S,Numerical},Augmentation{S,Numerical},AtomicPotential},
    args...;
    kwargs...,
) where {S}
    return _apply(x, interpolate_onto, args...; kwargs...)
end
