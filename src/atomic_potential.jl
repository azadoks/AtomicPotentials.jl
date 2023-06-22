struct AtomicPotential{S<:EvaluationSpace}
    identifier::Any
    symbol::Union{Nothing,Symbol}
    quantities::Vector{AbstractAtomicQuantity{S}}
end

get_quantities(a::AtomicPotential, T::Type{<:AbstractAtomicQuantity}) = filter(Base.Fix2(isa, T), a.quantities)
get_quantities(f::Function, a::AtomicPotential) = filter(f, a.quantities)
