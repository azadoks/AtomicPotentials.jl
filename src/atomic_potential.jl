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

function _apply(pot::AtomicPotential, f::Function, args...; kwargs...)
    states = map(pot.states) do states_l
        map(states_l) do state_ln
            return f(state_ln, args...; kwargs...)
        end
    end
    return AtomicPotential(
        pot.identifier,
        pot.symbol,
        f(pot.local_potential, args...; kwargs...),
        f(pot.nonlocal_potential, args...; kwargs...),
        f(pot.valence_density, args...; kwargs...),
        f(pot.core_density, args...; kwargs...),
        states,
        f(pot.augmentation, args...; kwargs...),
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
function max_r_length(pot::AtomicPotential)
    return maximum(
        max_r_length,
        [
            pot.local_potential,
            pot.nonlocal_potential,
            pot.valence_density,
            pot.core_density,
            pot.states,
            pot.augmentation,
        ];
        init=0,
    )
end
max_r_length(::Nothing) = 0
max_r_length(::AbstractAtomicQuantity{S,Analytical}) where {S} = 0
max_r_length(x::AbstractAtomicQuantity{S,Numerical}) where {S} = length(x.r)
function max_r_length(Vnl::NonLocalPotential{S,Numerical}) where {S}
    return maximum(βl -> maximum(max_r_length, βl; init=0), Vnl.β; init=0)
end
function max_r_length(aug::Augmentation)
    return maximum(Ql -> maximum(max_r_length, Ql; init=0), aug.Q; init=0)
end
max_r_length(x::AbstractVector) = maximum(max_r_length, x; init=0)
