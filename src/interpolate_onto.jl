function interpolate_onto(
    quantity::AbstractAtomicQuantity{S,Numerical}, r::AbstractVector
)::typeof(quantity) where {S}
    f = quantity.interpolator.(r)
    constructor = typeof(quantity)
    return _construct_similar_quantity(quantity, constructor; r=r, f=f)
end
function interpolate_onto(
    quantity::AbstractAtomicQuantity{S,Numerical}, Δr::Real
)::typeof(quantity) where {S}
    r_min, r_max = extrema(quantity.r)
    # Ensure that the endpoints of the radial mesh stay in the new radial mesh
    # by massaging Δr to be _at most_ the value provided, i.e. reduce the spacing
    # Δr until a range with an integer number of intervals containing the endpoints
    # can be made.
    r = range(; start=r_min, stop=r_max, length=ceil(Int, (r_max - r_min) / Δr))
    return interpolate_onto(quantity, r)
end

#* This is mainly useful for Hydrogenic orbitals which are analytical in real-space but
#* have no analytical Fourier transform. In this case, `fht(::HydrogenicProjector{Real,Analytical})`
#* should
#* 1. interpolate onto some reasonable radial mesh (using this function)
#*     a. use the requested q-grid to back-calculate a 'good' default r-mesh spacing
#*     b. set rmax to a value where the projector ≈ 0
#* 2. perform the numeric Fourier transform
#* 3. return a {FourierSpace,Numerical} quantity
# function interpolate_onto(
#     quantity::AbstractAtomicQuantity{S,Analytical}, r::AbstractVector
# ) where {S}
#     f = quantity.interpolator.(r)
#     f .= requires_r²(numerical_type) ? r .^ 2 .* f : f  # TODO: implement requires_r²
#     return _construct_numerical_quantity(quantity, r, f)  # TODO: implement _construct_numerical_quantity
# end
