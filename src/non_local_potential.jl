struct NonLocalPotential{
    S<:EvaluationSpace,A<:Analyticity,P<:AbstractKleinmanBylanderProjector{S,A}
}
    β::OffsetVector{Vector{P},Vector{Vector{P}}}
    D::OffsetVector{Matrix,Vector{Matrix}}
    # Q::Vector{AugmentationFunction{S,A}}  # TODO: ultrasoft potentials not implemented
    # q::AbstractVector{AbstractMatrix}  # TODO: ultrasoft potentials not implemented
end
Base.Broadcast.broadcastable(Vnl::NonLocalPotential) = Ref(Vnl)
function _apply(Vnl::NonLocalPotential, f::Function, args...; kwargs...)
    β = map(Vnl.β) do βl
        map(βl) do βln
            f(βln, args...; kwargs...)
        end
    end
    return NonLocalPotential(β, Vnl.D)
end
function fht(Vnl::NonLocalPotential{RealSpace}, args...; kwargs...)
    return _apply(Vnl, fht, args...; kwargs...)
end
function ifht(Vnl::NonLocalPotential{FourierSpace}, args...; kwargs...)
    return _apply(Vnl, ifht, args...; kwargs...)
end
function interpolate_onto(
    Vnl::NonLocalPotential{S,Numerical,P}, args...; kwargs...
)::NonLocalPotential{S,Numerical,P} where {S,P}
    return _apply(Vnl, interpolate_onto, args...; kwargs...)
end
function truncate(
    Vnl::NonLocalPotential{S,Numerical,P}, args...; kwargs...
)::NonLocalPotential{S,Numerical,P} where {S,P}
    return _apply(Vnl, truncate, args...; kwargs...)
end
