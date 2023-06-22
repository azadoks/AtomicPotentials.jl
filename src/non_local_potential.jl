struct NonLocalPotential{S,A} <: AbstractAtomicQuantity{S,A}
    β::Vector{AbstractKleinmanBylanderProjector{S,A}}
    D::AbstractMatrix
    Q::Vector{AugmentationFunction{S,A}}
    q::AbstractVector{AbstractMatrix}
end
