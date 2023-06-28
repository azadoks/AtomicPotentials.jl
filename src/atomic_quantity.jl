abstract type EvaluationSpace end
abstract type RealSpace <: EvaluationSpace end  # Is in real / direct space
abstract type FourierSpace <: EvaluationSpace end  # Is in Fourier / reciprocal space
dual_space_of(::Type{RealSpace})::Type = FourierSpace
dual_space_of(::Type{FourierSpace})::Type = RealSpace

abstract type Analyticity end
abstract type Analytical <: Analyticity end  # Has an analytical (inverse) Fourier transform
abstract type Numerical <: Analyticity end  # Must be numerically (inversely) Fourier transformed

abstract type AbstractAtomicQuantity{S<:EvaluationSpace,A<:Analyticity} end
angular_momentum(quantity::AbstractAtomicQuantity)::Int = quantity.l

function Base.show(io::IO, quantity::AbstractAtomicQuantity{S,Numerical}) where {S}
    r_min, r_max = extrema(quantity.r)
    r_length = length(quantity.r)
    r_eltype = eltype(quantity.r)
    return print(
        io, "$(typeof(quantity))(r=$r_length-element $(r_eltype)[$r_min...$r_max])"
    )
end

## Quantity evaluators
function (quantity::AbstractAtomicQuantity{S,Numerical})(x::T)::T where {S,T}
    return Interpolation.evaluate(quantity.interpolator, x)
end

## Methods for generically constructing new quantities with different type parameters
dual_space_of(::AbstractAtomicQuantity{S}) where {S<:EvaluationSpace} = dual_space_of(S)

# TODO: replace with `similar` + `copyto!` dispatches if possible (requires identical
# TODO: memory layout, I think, so maybe not?)
@generated function _bare_type_of(::Type{T}) where {T<:AbstractAtomicQuantity}
    return getfield(parentmodule(T), nameof(T))
end
function _bare_type_of(::T) where {T<:AbstractAtomicQuantity}
    return _bare_type_of(T)
end

_dual_type_of(::Type{Nothing}) = Nothing
function _dual_type_of(::Type{T}) where {T}
    S, A = T.parameters
    return _bare_type_of(T){dual_space_of(S),A}
end
_dual_type_of(::Nothing) = Nothing
function _dual_type_of(quantity::AbstractAtomicQuantity{S,A}) where {S,A}
    return _bare_type_of(quantity){dual_space_of(S),A}
end

function _construct_similar_quantity(
    quantity::AbstractAtomicQuantity,
    constructor::Type{<:AbstractAtomicQuantity}=typeof(quantity);
    kwargs...,
)
    # Take all of the property names and property values from the quantity instance and put
    # them into an ordered dictionary
    properties = OrderedDict(
        name => getproperty(quantity, name) for name in propertynames(quantity)
    )
    # Replace key-value pairs in the properties with the new values provided in the keyword
    # arguments
    for (key, value) in kwargs
        @assert key in keys(properties)
        properties[key] = value
    end
    return constructor(values(properties)...)
end

function _construct_dual_quantity(quantity::AbstractAtomicQuantity; kwargs...)
    return _construct_similar_quantity(quantity, _dual_type_of(quantity); kwargs...)
end
# TODO \
