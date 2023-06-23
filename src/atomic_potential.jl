struct AtomicPotential{S<:EvaluationSpace}
    identifier::Any
    symbol::Union{Nothing,Symbol}
    quantities::Vector{AbstractAtomicQuantity{S}}
end

function get_quantities(a::AtomicPotential, T::Type{<:AbstractAtomicQuantity})
    return filter(Base.Fix2(isa, T), a.quantities)
end
get_quantities(f::Function, a::AtomicPotential) = filter(f, a.quantities)
