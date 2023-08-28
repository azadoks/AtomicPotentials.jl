@doc raw"""
Non-local part of a Kleinman-Bylander separable form pseudo-potential.
"""
struct NonLocalPotential{
    S<:EvaluationSpace,
    VP<:AbstractVector{<:AbstractVector{<:AbstractQuantity{S}}},
    VC<:AbstractVector{<:AbstractMatrix{<:Real}},
    TA<:Union{Nothing,Augmentation{S}},
}
    "Non-local projectors."
    projectors::VP
    "Kleinman-Bylander energies / coupling constants."
    coupling::VC
    "Augmentation contribution to the charge density."
    augmentation::TA
end
function NonLocalPotential(projectors, coupling, augmentation=nothing)
    return NonLocalPotential(projectors, coupling, augmentation)
end

Base.broadcastable(qty::NonLocalPotential) = Ref(qty)
Base.isempty(qty::NonLocalPotential) = false

function Base.show(io::IO, ::MIME"text/plain", quantity::NonLocalPotential)
    println(io, "NonLocalPotential")
    for (l_index, (D_l, P_l)) in enumerate(zip(quantity.coupling, quantity.projectors))
        println(io, "\tl=$(l_index-1)")
        println(io, "\t\tD: ")
        println(io, "\t\t\t$(D_l)")
        println(io, "\t\tP:")
        print("\t\t\t")
        for P_li in P_l
            print(io, "$(P_li) ")
        end
        println()
    end
end

function Base.convert(::Type{T}, x::NonLocalPotential) where {T<:Real}
    return NonLocalPotential(
        map(P_l -> map(P -> convert(T, P), P_l), x.projectors),
        map(D_l -> map(D -> convert(T, D), D_l), x.coupling),
        _convert_optional(T, x.augmentation),
    )
end

function Adapt.adapt(to, x::NonLocalPotential)
    return NonLocalPotential(
        adapt(to, x.projectors), adapt(to, x.coupling), adapt(to, augmentation)
    )
end

function n_projector_radials(quantity::NonLocalPotential, l::Integer)
    l > maximum_angular_momentum(quantity) && return 0
    return length(quantity.projectors[l + 1])
end
n_projector_radials(quantity::NonLocalPotential) = length(quantity.projectors)

function n_projector_angulars(quantity::NonLocalPotential, l::Integer)
    return n_projector_radials(quantity, l) * (2l + 1)
end
function n_projector_angulars(quantity::NonLocalPotential)
    return sum(projector -> 2 * angular_momentum(projector) + 1, quantity.projectors)
end

function maximum_angular_momentum(quantity::NonLocalPotential)
    return maximum(
        projectors_l -> maximum(angular_momentum, projectors_l), quantity.projectors
    )
end

function rft(non_local::NonLocalPotential{RealSpace}, q::AbstractVector; kwargs...)
    projectors = map(f_l -> map(f -> rft(f, q; kwargs...), f_l), non_local.projectors)
    augmentation = rft(non_local.augmentation, q; kwargs...)
    return NonLocalPotential(projectors, non_local.coupling, augmentation)
end

function irft(non_local::NonLocalPotential{FourierSpace}, r::AbstractVector; kwargs...)
    projectors = map(F_l -> map(F -> irft(F, r; kwargs...), F_l), non_local.projectors)
    augmentation = irft(non_local.augmentation, r; kwargs...)
    return NonLocalPotential(projectors, non_local.coupling, augmentation)
end
