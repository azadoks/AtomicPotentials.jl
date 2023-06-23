struct NonLocalPotential{S}
    β::Vector{AbstractKleinmanBylanderProjector{S}}
    D::AbstractMatrix
    # Q::Vector{AugmentationFunction{S,A}}  # TODO: ultrasoft potentials not implemented
    # q::AbstractVector{AbstractMatrix}  # TODO: ultrasoft potentials not implemented
end

function fht(
    quantity::NonLocalPotential{RealSpace}, q::AbstractVector, method::QuadratureMethod
)
    β = [fht(βi, q, method) for βi in quantity.β]
    return NonLocalPotential{FourierSpace}(β, quantity.D)
end
function fht(quantity::NonLocalPotential{RealSpace})
    β = [fht(βi) for βi in quantity.β]
    return NonLocalPotential{FourierSpace}(β, quantity.D)
end

function ifht(
    quantity::NonLocalPotential{FourierSpace}, r::AbstractVector, method::QuadratureMethod
)
    β = [ifht(βi, r, method) for βi in quantity.β]
    return NonLocalPotential{RealSpace}(β, quantity.D)
end
function ifht(quantity::NonLocalPotential{FourierSpace})
    β = [ifht(βi) for βi in quantity.β]
    return NonLocalPotential{RealSpace}(β, quantity.D)
end
