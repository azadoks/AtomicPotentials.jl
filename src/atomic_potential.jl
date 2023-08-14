struct AtomicPotential{
    S<:EvaluationSpace,
    L<:AbstractLocalPotential{S},
    NL<:Union{Nothing,NonlocalPotential{S}},
    VD<:Union{Nothing,AbstractChargeDensity{S}},
    CD<:Union{Nothing,AbstractChargeDensity{S}},
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

# Constructor to fetch types
function AtomicPotential(
    identifier,
    symbol,
    local_potential,
    nonlocal_potential,
    valence_density,
    core_density,
    states,
)
    S = typeof(local_potential).parameters[1]
    return AtomicPotential{
        S,
        typeof(local_potential),
        typeof(nonlocal_potential),
        typeof(valence_density),
        typeof(core_density),
        eltype(eltype((states))),
    }(
        identifier,
        symbol,
        local_potential,
        nonlocal_potential,
        valence_density,
        core_density,
        states,
    )
end

# Minimal constructor
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

Base.Broadcast.broadcastable(pot::AtomicPotential) = Ref(pot)

identifier(pot::AtomicPotential) = pot.identifier
symbol(pot::AtomicPotential) = pot.symbol
local_potential(pot::AtomicPotential) = pot.local_potential
nonlocal_potential(pot::AtomicPotential) = pot.nonlocal_potential
valence_density(pot::AtomicPotential) = pot.valence_density
core_density(pot::AtomicPotential) = pot.core_density
states(pot::AtomicPotential) = pot.states

charge_ionic(pot::AtomicPotential) = charge_ionic(local_potential(pot))
charge_nuclear(pot::AtomicPotential) = PeriodicTable.elements[pot.symbol].number

n_elec_valence(pot::AtomicPotential) = round(Int, charge_ionic(pot))
n_elec_core(pot::AtomicPotential) = round(Int, charge_nuclear(pot) - charge_ionic(pot))

function n_kbprojectors_radial(pot::AtomicPotential, args...)
    return n_kbprojectors_radia(nonlocal_potential(pot), args...)
end
count_n_proj_radial = n_kbprojectors_radial  # For DFTK transition

function n_kbprojectors_angular(pot::AtomicPotential, args...)
    return n_kbprojectors_angula(nonlocal_potential(pot), args...)
end
function n_kbprojectors_angular(
    pots::AbstractVector{<:AtomicPotential}, positions::AbstractVector
)
    return dot(n_kbprojectors_angular.(pots), length.(positions))
end
count_n_proj = n_kbprojectors_angular  # For DFTK transition

lmin(pot::AtomicPotential) = lmin(nonlocal_potential(pot))
lmin(::Nothing) = -1

lmax(pot::AtomicPotential) = lmax(nonlocal_potential(pot))
lmax(::Nothing) = -2

# Get a unit range from the minimum to maximum angular momentum in an atomic
# potential (effectively, in its non-local potential, if present)
angular_momenta(pot::AtomicPotential) = angular_momenta(nonlocal_potential(pot))

function energy_correction(T::Type{<:Real}, pot::AtomicPotential, args...; kwargs...)
    return energy_correction(T, local_potential(pot), args...; kwargs...)
end

# Find the maximum length of the radial mesh vector `r` over all the quantities
# in an atomic potential.
# Used for allocating correctly-sized working vectors.
function max_r_length(pot::AtomicPotential)
    return maximum(
        max_r_length,
        [
            pot.local_potential,
            pot.nonlocal_potential,
            pot.valence_density,
            pot.core_density,
            pot.states,
        ];
        init=0,
    )
end
max_r_length(::Nothing) = 0
max_r_length(x::AbstractVector) = maximum(max_r_length, x; init=0)
max_r_length(x::AbstractAtomicQuantity{S,Numerical}) where {S} = length(x.r)
max_r_length(::AbstractAtomicQuantity{S,Analytical}) where {S} = 0
function max_r_length(aug::Augmentation)
    return maximum(Ql -> maximum(max_r_length, Ql; init=0), aug.Q; init=0)
end
function max_r_length(Vnl::NonlocalPotential{S,Numerical}) where {S}
    max_n_β = maximum(βl -> maximum(max_r_length, βl; init=0), Vnl.β; init=0)
    max_n_aug = max_r_length(Vnl.augmentation)
    return max(max_n_β, max_n_aug)
end

function Base.show(io::IO, ::MIME"text/plain", pot::AtomicPotential)
    @printf io "%032s: %s\n" "identifier" pot.identifier
    @printf io "%032s: %s\n" "element" pot.symbol
    @printf io "%032s: %s\n" "local potential" pot.local_potential
    @printf io "%032s: %s\n" "non-local potential" pot.nonlocal_potential
    @printf io "%032s: %s\n" "valence density" pot.valence_density
    @printf io "%032s: %s\n" "core density" pot.core_density
    @printf io "%032s: %s\n" "states" pot.states
end

function Base.map(f, pot::AtomicPotential, args...; kwargs...)
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
    )
end

_apply(pot::AtomicPotential, f, args...; kwargs...) = map(f, pot, args...; kwargs...)
