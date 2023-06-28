struct NonLocalPotential{
    S<:EvaluationSpace,A<:Analyticity,P<:AbstractKleinmanBylanderProjector{S,A}
}
    β::OffsetVector{Vector{P},Vector{Vector{P}}}
    D::OffsetVector{Matrix,Vector{Matrix}}
end

Base.Broadcast.broadcastable(Vnl::NonLocalPotential) = Ref(Vnl)

function Base.show(io::IO, Vnl::NonLocalPotential)
    n_projectors = sum(length, Vnl.β)
    return print(io, "$(typeof(Vnl))($(n_projectors)x$(eltype(eltype(Vnl.β))))")
end

function _apply(Vnl::NonLocalPotential, f::Function, args...; kwargs...)
    β = map(Vnl.β) do βl
        map(βl) do βln
            f(βln, args...; kwargs...)
        end
    end
    return NonLocalPotential(β, Vnl.D)
end
