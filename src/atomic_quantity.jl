using OrderedCollections

abstract type EvaluationSpace end
abstract type RealSpace <: EvaluationSpace end  # Is in real / direct space
abstract type FourierSpace <: EvaluationSpace end  # Is in Fourier / reciprocal space
dual_space_of(::Type{RealSpace})::Type = FourierSpace
dual_space_of(::Type{FourierSpace})::Type = RealSpace

abstract type Analyticity end
abstract type Analytical <: Analyticity end  # Has an analytical (inverse) Fourier transform
abstract type Numerical <: Analyticity end  # Must be numerically (inversely) Fourier transformed

abstract type AbstractAtomicQuantity{S<:EvaluationSpace, A<:Analyticity} end
# angular_momentum(quantity::AbstractAtomicQuantity)::Int

## Quantity evaluators
function (quantity::AbstractAtomicQuantity{S,Numerical})(x) where {S}
    return quantity.interpolator(x)
end

## Methods for generically constructing new quantities with different type parameters
dual_space_of(::AbstractAtomicQuantity{S}) where {S<:EvaluationSpace} = dual_space_of(S)

# TODO: replace with `similar` + `copyto!` dispatches if possible (requires identical
# TODO: memory layout, I think)
@generated function bare_constructor_of(::Type{T}) where {T<:AbstractAtomicQuantity}
    return getfield(parentmodule(T), nameof(T))
end
function bare_constructor_of(::T) where {T<:AbstractAtomicQuantity}
    return bare_constructor_of(T)
end

# TODO: function dual_constructor_of(::Type{T}) end
function dual_constructor_of(quantity::AbstractAtomicQuantity{S,A}) where {S,A}
    return bare_constructor_of(quantity){dual_space_of(S),A}
end

function construct_dual_quantity(quantity::AbstractAtomicQuantity{S,A}; kwargs...) where {S<:EvaluationSpace,A<:Analyticity}
    # Take all of the property names and property values from the quantity instance and put
    # them into an ordered dictionary
    properties = OrderedDict(name => getproperty(quantity, name) for name in propertynames(quantity))
    # Replace key-value pairs in the properties with the new values provided in the keyword
    # arguments
    for (key, value) in kwargs
        @assert key in keys(properties)
        properties[key] = value
    end
    # Get the constructor for the type of quantity provided and use it to construct a new
    # instance with the dual evaluation space type parameter
    return dual_constructor_of(quantity)(values(properties)...)
end
# TODO \
