# (inverse) radial fourier transform fallbacks so that quantities / containers which are
# optional (like the augmentation in the non-local potential don't require loads of
# if-else statements
rft(::Nothing, args...; kwargs...) = nothing
irft(::Nothing, args...; kwargs...) = nothing

# dispatches for vectors of quantities (e.g. orbitals, projectors)
# function rft(
#     quantities::AbstractVector{<:Union{Nothing,AbstractQuantity{RealSpace}}},
#     q::AbstractVector,;
#     kwargs...,
# )
#     return map(quantity -> rft(quantity, q; kwargs...), quantities)
# end

# function irft(
#     quantities::AbstractVector{<:Union{Nothing,AbstractQuantity{RealSpace}}},
#     r::AbstractVector;
#     kwargs...,
# )
#     return map(quantity -> irft(quantity, r; kwargs...), quantities)
# end

@doc raw"""
Ultrasoft augmentation contribution to the charge density.
"""
struct Augmentation{
    S<:EvaluationSpace,
    TF<:AbstractQuantity{S},
    VF<:AbstractVector{TF},
    TC<:Real,
    AC<:AbstractArray{Real,3},
}
    "Augmentation functions."
    functions::VF
    "Augmentation charges / coupling constants"
    coupling::AC
end

function rft(augmentation::Augmentation{RealSpace}, q::AbstractVector; kwargs...)
    functions = map(f -> rft(f, q; kwargs...), augmentation.functions)
    return Augmentation(functions, augmentation.coupling)
end

function irft(augmentation::Augmentation{FourierSpace}, r::AbstractVector; kwargs...)
    functions = map(f -> irft(f, r; kwargs...), augmentation.functions)
    return Augmentation(functions, augmentation.coupling)
end

@doc raw"""
Non-local part of a Kleinman-Bylander separable form pseudo-potential.
"""
struct NonLocalPotential{
    S<:EvaluationSpace,
    TP<:AbstractQuantity{S},
    VP<:AbstractVector{TP},
    TC<:Real,
    MC<:AbstractMatrix{TC},
    TA<:Union{Nothing,Augmentation{S}},
}
    "Non-local projectors."
    projectors::VP
    "Kleinman-Bylander energies / coupling constants."
    coupling::MC
    "Augmentation contribution to the charge density."
    augmentation::TA
end
function NonLocalPotential(
    projectors::VP, coupling, augmentation=nothing
) where {S<:EvaluationSpace,TP<:AbstractQuantity{S},VP<:AbstractVector{TP}}
    return NonLocalPotential(projectors, coupling, augmentation)
end

function rft(non_local::NonLocalPotential{RealSpace}, q::AbstractVector; kwargs...)
    projectors = map(f -> rft(f, q; kwargs...), non_local.projectors)
    augmentation = rft(non_local.augmentation, q; kwargs...)
    return NonLocalPotential(projectors, non_local.coupling, augmentation)
end

function irft(non_local::NonLocalPotential{FourierSpace}, r::AbstractVector; kwargs...)
    projectors = map(f -> irft(f, r; kwargs...), non_local.projectors)
    augmentation = irft(non_local.augmentation, r; kwargs...)
    return NonLocalPotential(projectors, non_local.coupling, augmentation)
end

@doc raw"""
Atomic potential.
"""
struct AtomicPotential{
    S<:EvaluationSpace,
    TLP<:AbstractLocalPotential{S},
    TNL<:Union{Nothing,NonLocalPotential{S}},
    TCD<:Union{Nothing,AbstractQuantity{S}},
    TID<:Union{Nothing,AbstractQuantity{S}},
    TO<:Union{Nothing,AbstractQuantity{S}},
    VO<:AbstractVector{TO},
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
    orbitals::VO
end
function AtomicPotential(
    local_potential::AbstractLocalPotential{S},
    non_local_potential=nothing,
    core_charge_density=nothing,
    ionic_charge_density=nothing,
    orbitals=Nothing[],
) where {S<:EvaluationSpace}
    return AtomicPotential{S}(
        local_potential,
        non_local_potential,
        core_charge_density,
        ionic_charge_density,
        orbitals,
    )
end

function rft(potential::AtomicPotential{RealSpace}, q::AbstractVector; kwargs...)
    local_potential = rft(potential.local_potential, q; kwargs...)
    non_local_potential = rft(potential.non_local_potential, q; kwargs...)
    core_charge_density = rft(potential.core_charge_density, q; kwargs...)
    ionic_charge = rft(potential.ionic_charge_density, q; kwargs...)
    orbitals = map(f -> rft(f, q; kwargs...), potential.orbitals)
    return AtomicPotential(
        local_potential, non_local_potential, core_charge_density, ionic_charge, orbitals
    )
end

function irft(potential::AtomicPotential{FourierSpace}, r::AbstractVector; kwargs...)
    local_potential = irft(potential.local_potential, r; kwargs...)
    non_local_potential = irft(potential.non_local_potential, r; kwargs...)
    core_charge_density = irft(potential.core_charge_density, r; kwargs...)
    ionic_charge = irft(potential.ionic_charge_density, r; kwargs...)
    orbitals = map(f -> irft(f, r; kwargs...), potential.orbitals)
    return AtomicPotential(
        local_potential, non_local_potential, core_charge_density, ionic_charge, orbitals
    )
end
