struct NonlocalPotential{
    S<:EvaluationSpace,
    A<:Analyticity,
    P<:AbstractKleinmanBylanderProjector{S,A},
    AUG<:Union{Nothing,Augmentation{S}},
}
    β::OffsetVector{Vector{P},Vector{Vector{P}}}
    D::OffsetVector{Matrix,Vector{Matrix}}
    augmentation::AUG
end

Base.Broadcast.broadcastable(Vnl::NonlocalPotential) = Ref(Vnl)

function Base.show(io::IO, Vnl::NonlocalPotential)
    n_projectors = sum(length, Vnl.β)
    return print(io, "$(typeof(Vnl))($(n_projectors)x$(eltype(eltype(Vnl.β))))")
end

function _apply(Vnl::NonlocalPotential, f::Function, args...; kwargs...)
    β = map(Vnl.β) do βl
        map(βl) do βln
            f(βln, args...; kwargs...)
        end
    end
    augmentation = _apply(Vnl.augmentation, f, args...; kwargs...)
    return NonlocalPotential(β, Vnl.D, augmentation)
end

kbprojectors(Vnl::NonlocalPotential) = Vnl.β
kbprojectors(Vnl::NonlocalPotential, l) = Vnl.β[l]
kbprojectors(Vnl::NonlocalPotential, l, n) = kbprojectors(Vnl, l)[n]

kbenergies(Vnl::NonlocalPotential) = Vnl.D
kbenergies(Vnl::NonlocalPotential, l) = Vnl.D[l]
kbenergies(Vnl::NonlocalPotential, l, n) = kbenergies(Vnl, l)[n]

augmentation(Vnl::NonlocalPotential) = Vnl.augmentation

n_kbprojectors_radial(Vnl::NonlocalPotential, l::Integer) = length(kbprojectors(Vnl, l))
n_kbprojectors_radial(Vnl::NonlocalPotential) = map(length, kbprojectors(Vnl))
n_kbprojectors_radial(::Nothing) = 0

function n_kbprojectors_angular(Vnl::NonlocalPotential, l::Integer)
    return n_kbprojectors_radial(Vnl, l) * (2l + 1)
end
function n_kbprojectors_angular(Vnl::NonlocalPotential)
    return sum(Base.Fix1(n_kbprojectors_angular, Vnl), angular_momenta(Vnl))
end
n_kbprojectors_angular(::Nothing) = 0

lmin(Vnl::NonlocalPotential) = firstindex(Vnl.β)
lmax(Vnl::NonlocalPotential) = lastindex(Vnl.β)

angular_momenta(Vnl::NonlocalPotential) = lmin(Vnl):lmax(Vnl)
