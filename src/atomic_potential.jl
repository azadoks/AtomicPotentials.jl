struct AtomicPotential{
    S<:EvaluationSpace,
    L<:AbstractLocalPotential{S},
    NL<:Union{Nothing,AbstractNonLocalPotential{S}},
    VD<:Union{Nothing,AbstractValenceChargeDensity{S}},
    CD<:Union{Nothing,AbstractCoreChargeDensity{S}},
    SP<:Union{Nothing,AbstractStateProjector{S}},
}
    identifier::Any
    symbol::Union{Nothing,Symbol}
    local_potential::L
    nonlocal_potential::NL
    valence_density::VD
    core_density::CD
    states::OffsetVector{Vector{SP},Vector{Vector{SP}}}
end
function AtomicPotential(
    local_potential;
    identifier="",
    symbol=:X,
    nonlocal_potential=nothing,
    valence_density=nothing,
    core_density=nothing,
    states=OffsetVector(Vector{Nothing}[], 0:-1),
)
    return AtomicPotential(
        identifier,
        symbol,
        local_potential,
        nonlocal_potential,
        valence_density,
        core_density,
        states,
    )
end
Base.Broadcast.broadcastable(potential::AtomicPotential) = Ref(potential)
function Base.show(io::IO, ::MIME"text/plain", pot::AtomicPotential)
    @printf io "%032s: %s\n" "identifier" pot.identifier
    @printf io "%032s: %s\n" "element" pot.symbol
    @printf io "%032s: %s\n" "local potential" pot.local_potential
    @printf io "%032s: %s\n" "non-local potential" pot.nonlocal_potential
    @printf io "%032s: %s\n" "valence density" pot.valence_density
    @printf io "%032s: %s\n" "core density" pot.core_density
    @printf io "%032s: %s\n" "states" pot.states
end

function _apply(potential::AtomicPotential, f::Function, args...; kwargs...)
    states = map(potential.states) do states_l
        map(states_l) do state_ln
            return f(state_ln, args...; kwargs...)
        end
    end
    return AtomicPotential(
        potential.identifier,
        potential.symbol,
        f(potential.local_potential, args...; kwargs...),
        f(potential.nonlocal_potential, args...; kwargs...),
        f(potential.valence_density, args...; kwargs...),
        f(potential.core_density, args...; kwargs...),
        states,
    )
end
# TODO: metaprogramming?
function fht(potential::AtomicPotential{RealSpace}, args...; kwargs...)
    return _apply(potential, fht, args...; kwargs...)
end
function ifht(potential::AtomicPotential{FourierSpace}, args...; kwargs...)
    return _apply(potential, ifht, args...; kwargs...)
end
function interpolate_onto(potential::AP, args...; kwargs...)::AP where {AP<:AtomicPotential}
    return _apply(potential, interpolate_onto, args...; kwargs...)
end
function truncate(potential::AP, args...; kwargs...)::AP where {AP<:AtomicPotential}
    return _apply(potential, truncate, args...; kwargs...)
end
