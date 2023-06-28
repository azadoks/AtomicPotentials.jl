struct AtomicPotential{
    S<:EvaluationSpace,
    L<:AbstractLocalPotential{S},
    NL<:Union{Nothing,NonLocalPotential{S}},
    VD<:Union{Nothing,AbstractChargeDensity{S}},
    CD<:Union{Nothing,AbstractChargeDensity{S}},
    SP<:Union{Nothing,AbstractStateProjector{S}},
    AUG<:Union{Nothing,Augmentation{S}},
}
    identifier::Any
    symbol::Union{Nothing,Symbol}
    local_potential::L
    nonlocal_potential::NL
    valence_density::VD
    core_density::CD
    states::OffsetVector{Vector{SP},Vector{Vector{SP}}}
    augmentation::AUG
end

function AtomicPotential(
    local_potential;
    identifier="",
    symbol=:X,
    nonlocal_potential=nothing,
    valence_density=nothing,
    core_density=nothing,
    states=OffsetVector(Vector{Nothing}[], 0:-1),
    augmentation=nothing,
)
    return AtomicPotential(
        identifier,
        symbol,
        local_potential,
        nonlocal_potential,
        valence_density,
        core_density,
        states,
        augmentation,
    )
end

for (op, Tin, Tout) in ((:ht, :RealSpace, :FourierSpace), (:iht, :RealSpace, :FourierSpace))
    eval(
        quote
            function $(op)(
                pot::AtomicPotential{$(Tin)}, args...; kwargs...
            )::AtomicPotential{$(Tout)}
                states = map(pot.states) do states_l
                    map(states_l) do state_ln
                        return $(op)(state_ln, args...; kwargs...)
                    end
                end
                return AtomicPotential(
                    pot.identifier,
                    pot.symbol,
                    $(op)(pot.local_potential, args...; kwargs...),
                    $(op)(pot.nonlocal_potential, args...; kwargs...),
                    $(op)(pot.valence_density, args...; kwargs...),
                    $(op)(pot.core_density, args...; kwargs...),
                    states,
                    $(op)(pot.augmentation, args...; kwargs...),
                )
            end
        end,
    )
end

Base.Broadcast.broadcastable(pot::AtomicPotential) = Ref(pot)

function Base.show(io::IO, ::MIME"text/plain", pot::AtomicPotential)
    @printf io "%032s: %s\n" "identifier" pot.identifier
    @printf io "%032s: %s\n" "element" pot.symbol
    @printf io "%032s: %s\n" "local potential" pot.local_potential
    @printf io "%032s: %s\n" "non-local potential" pot.nonlocal_potential
    @printf io "%032s: %s\n" "valence density" pot.valence_density
    @printf io "%032s: %s\n" "core density" pot.core_density
    @printf io "%032s: %s\n" "states" pot.states
    @printf io "%032s: %s\n" "augmentation" pot.augmentation
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
        f(potential.augmentation, args...; kwargs...),
    )
end

charge_ionic(pot::AtomicPotential) = charge_ionic(pot.local_potential)
charge_nuclear(pot::AtomicPotential) = PeriodicTable.elements[pot.symbol].number
n_elec_valence(pot::AtomicPotential) = round(Int, charge_ionic(pot))
n_elec_core(pot::AtomicPotential) = round(Int, charge_nuclear(pot) - charge_ionic(pot))
function energy_correction(T::Type{<:Real}, pot::AtomicPotential, args...; kwargs...)
    return energy_correction(T, pot.local_potential, args...; kwargs...)
end
function count_n_proj_radial(pot::AtomicPotential, args...)
    return count_n_proj_radial(pot.nonlocal_potential, args...)
end
count_n_proj(pot::AtomicPotential, args...) = count_n_proj(pot.nonlocal_potential, args...)
function count_n_proj(pots::AbstractVector{<:AtomicPotential}, positions::AbstractVector)
    return dot(count_n_proj.(pots), length.(positions))
end
lmax(pot::AtomicPotential) = lmax(pot.nonlocal_potential)
angular_momenta(pot::AtomicPotential) = angular_momenta(pot.nonlocal_potential)
