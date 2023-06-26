abstract type AbstractNonLocalPotential{S<:EvaluationSpace,A<:Analyticity} end
Base.Broadcast.broadcastable(Vnl::AbstractNonLocalPotential) = Ref(Vnl)
function fht(Vnl::AbstractNonLocalPotential{RealSpace}, args...; kwargs...)
    return _apply(Vnl, fht, args...; kwargs...)
end
function ifht(Vnl::AbstractNonLocalPotential{FourierSpace}, args...; kwargs...)
    return _apply(Vnl, ifht, args...; kwargs...)
end
function interpolate_onto(
    Vnl::AbstractNonLocalPotential{S,Numerical}, args...; kwargs...
) where {S}
    return _apply(Vnl, interpolate_onto, args...; kwargs...)
end
function truncate(Vnl::AbstractNonLocalPotential{S,Numerical}, args...; kwargs...) where {S}
    return _apply(Vnl, truncate, args...; kwargs...)
end

struct NormConservingNonLocalPotential{
    S<:EvaluationSpace,A<:Analyticity,P<:AbstractKleinmanBylanderProjector{S,A}
} <: AbstractNonLocalPotential{S,A}
    β::OffsetVector{Vector{P},Vector{Vector{P}}}
    D::OffsetVector{Matrix,Vector{Matrix}}
end
function Base.show(io::IO, Vnl::NormConservingNonLocalPotential)
    n_projectors = sum(length, Vnl.β)
    return print(io, "$(typeof(Vnl))($(n_projectors)x$(eltype(eltype(Vnl.β))))")
end
function _apply(Vnl::NormConservingNonLocalPotential, f::Function, args...; kwargs...)
    β = map(Vnl.β) do βl
        map(βl) do βln
            f(βln, args...; kwargs...)
        end
    end
    return NormConservingNonLocalPotential(β, Vnl.D)
end

struct UltrasoftNonLocalPotential{
    S<:EvaluationSpace,
    A<:Analyticity,
    P<:AbstractKleinmanBylanderProjector{S,A},
    QT<:AugmentationFunction{S,A},
} <: AbstractNonLocalPotential{S,A}
    β::OffsetVector{Vector{P},Vector{Vector{P}}}
    D::OffsetVector{Matrix,Vector{Matrix}}
    Q::OffsetVector{Matrix{QT},Vector{Matrix{QT}}}
    q::OffsetVector{Matrix,Vector{Matrix}}
end
function Base.show(io::IO, Vnl::UltrasoftNonLocalPotential)
    n_projectors = sum(length, Vnl.β)
    n_augmentation = sum(length, Vnl.Q)
    return print(
        io,
        "$(typeof(Vnl))($(n_projectors)x$(eltype(eltype(Vnl.β))),",
        " $(n_augmentation)x$(eltype(eltype(Vnl.Q))))",
    )
end
function _apply(Vnl::UltrasoftNonLocalPotential, f::Function, args...; kwargs...)
    β = map(Vnl.β) do βl
        map(βl) do βln
            f(βln, args...; kwargs...)
        end
    end
    Q = map(Vnl.Q) do Ql
        map(Ql) do Qlnm
            f(Qlnm, args...; kwargs...)
        end
    end
    return UltrasoftNonLocalPotential(β, Vnl.D, Q, Vnl.q)
end
