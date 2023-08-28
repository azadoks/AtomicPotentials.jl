@doc raw"""
Atomic potential.
"""
struct AtomicPotential{
    S<:EvaluationSpace,
    TLP<:AbstractLocalPotential{S},
    TNL<:Union{Nothing,NonLocalPotential{S}},
    TCD<:Union{Nothing,AbstractQuantity{S}},
    TID<:Union{Nothing,AbstractQuantity{S}},
    VVO<:AbstractVector{<:Union{Nothing,AbstractVector{<:AbstractQuantity{S}}}},
}
    "Local potential."
    local_potential::TLP
    "Non-local potential in Kleinman-Bylander separable form."
    non_local_potential::TNL
    "Core charge density for non-linear core correction."
    core_charge_density::TCD
    "Ionic / pseudo-valence charge density for charge superposition / projection."
    ionic_charge_density::TID
    "Atomic / pseudo-atomic orbitals for linear combination / projection."
    orbitals::VVO
end

## Constructors

function AtomicPotential(
    local_potential,
    non_local_potential=nothing,
    core_charge_density=nothing,
    ionic_charge_density=nothing,
    orbitals=Nothing[],
)
    return AtomicPotential(
        local_potential,
        non_local_potential,
        core_charge_density,
        ionic_charge_density,
        orbitals,
    )
end

## Base dispatches

Base.broadcastable(potential::AtomicPotential) = Ref(potential)
Base.isempty(::AtomicPotential) = false

function Base.convert(::Type{T}, x::AtomicPotential) where {T<:Real}
    local_ = convert(T, x.local_potential)
    non_local = _convert_optional(T, x.non_local_potential)
    core = _convert_optional(T, x.core_charge_density)
    ionic = _convert_optional(T, x.ionic_charge_density)
    orbitals = map(o_l -> map(o -> convert(T, o), o_l), x.orbitals)
    return AtomicPotential(local_, non_local, core, ionic, orbitals)
end

# TODO: Base.show

## Adapt

function Adapt.adapt_structure(to, x::AtomicPotential)
    local_ = adapt(to, x.local_potential)
    non_local = adapt(to, x.non_local_potential)
    core = adapt(to, x.core_charge_density)
    ionic = adapt(to, x.ionic_charge_density)
    orbitals = adapt(to, map(orbital -> adapt(to, orbital), x.orbitals))
    return AtomicPotential(local_, non_local, core, ionic, orbitals)
end

## Existence checks

function hasquantity(potential::AtomicPotential, name::Symbol)
    !in(name, propertynames(potential)) && return false
    quantity = getproperty(potential, name)
    return !isnothing(quantity) && !isempty(quantity)
end

## Local potential passthroughs

function energy_correction(potential::AtomicPotential{RealSpace}; kwargs...)
    return energy_correction(potential.local_potential; kwargs...)
end

## Non-local potential passthroughs

function n_projector_radials(potential::AtomicPotential, args...)
    vnl = potential.non_local_potential
    return isnothing(vnl) ? 0 : n_projector_radials(vnl, args...)
end

function n_projector_angulars(potential::AtomicPotential, args...)
    vnl = potential.non_local_potential
    return isnothing(vnl) ? 0 : n_projector_angulars(vnl, args...)
end

function maximum_angular_momentum(potential::AtomicPotential)
    vnl = potential.non_local_potential
    return isnothing(vnl) ? 0 : maximum_angular_momentum(vnl)
end

## Orbital functions

function n_orbital_radials(potential::AtomicPotential, l::Integer)
    l + 1 > length(potential.orbitals) && return 0
    return length(potential.orbitals[l + 1])
end

function n_orbital_radials(potential::AtomicPotential)
    return sum(length, potential.orbitals; init=0)
end

function n_orbital_angulars(potential::AtomicPotential, l::Integer)
    return n_orbital_radials(potential, l) * (2l + 1)
end

function n_orbital_angulars(potential::AtomicPotential)
    return sum(l -> n_orbital_angulars(potential, l), 0:maximum_angular_momentum(potential))
end

## Fourier transforms

function rft(potential::AtomicPotential{RealSpace}, q::AbstractVector; kwargs...)
    local_ = rft(potential.local_potential, q; kwargs...)
    non_local = rft(potential.non_local_potential, q; kwargs...)
    core = rft(potential.core_charge_density, q; kwargs...)
    ionic = rft(potential.ionic_charge_density, q; kwargs...)
    orbitals = map(f_l -> map(f -> rft(f, q; kwargs...), f_l), potential.orbitals)
    return AtomicPotential(local_, non_local, core, ionic, orbitals)
end

function irft(potential::AtomicPotential{FourierSpace}, r::AbstractVector; kwargs...)
    local_ = irft(potential.local_potential, r; kwargs...)
    non_local = irft(potential.non_local_potential, r; kwargs...)
    core = irft(potential.core_charge_density, r; kwargs...)
    ionic = irft(potential.ionic_charge_density, r; kwargs...)
    orbitals = map(F_l -> map(F -> irft(F, r; kwargs...), F_l), potential.orbitals)
    return AtomicPotential(local_, non_local, core, ionic, orbitals)
end
