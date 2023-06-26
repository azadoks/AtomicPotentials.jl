struct NonLocalPotential{
    S<:EvaluationSpace,A<:Analyticity,P<:AbstractKleinmanBylanderProjector{S,A}
} <: AbstractAtomicQuantity{S,A}
    β::OffsetVector{Vector{P},Vector{Vector{P}}}
    D::OffsetVector{Matrix,Vector{Matrix}}
    # Q::Vector{AugmentationFunction{S,A}}  # TODO: ultrasoft potentials not implemented
    # q::AbstractVector{AbstractMatrix}  # TODO: ultrasoft potentials not implemented
end
function Base.show(io::IO, Vnl::NonLocalPotential)
    return print(io, "$(typeof(Vnl))(β=$(Vnl.β), D=$(Vnl.D))")
end
function _dual_type_of(::Type{T}) where {T<:NonLocalPotential}
    S, A, P = T.parameters
    return _bare_type_of(T){dual_space_of(S),A,_dual_type_of(P)}
end
function _dual_type_of(quantity::NonLocalPotential{S,A,P}) where {S,A,P}
    return _bare_type_of(quantity){dual_space_of(S),A,_dual_type_of(P)}
end
function _apply(Vnl::NonLocalPotential, f::Function, args...; kwargs...)
    β = map(Vnl.β) do βl
        map(βl) do βln
            f(βln, args...; kwargs...)
        end
    end
    return NonLocalPotential(β, Vnl.D)
end
function fht(
    Vnl::NonLocalPotential{RealSpace,A,P}, q::AbstractVector, args...; kwargs...
)::NonLocalPotential{FourierSpace,A,_dual_type_of(P)} where {A,P}
    return _apply(Vnl, fht, q::AbstractVector, args...; kwargs...)
end
function ifht(
    Vnl::NonLocalPotential{FourierSpace,A,P}, r::AbstractVector, args...; kwargs...
)::NonLocalPotential{RealSpace,A,_dual_type_of(P)} where {A,P}
    return _apply(Vnl, ifht, r::AbstractVector, args...; kwargs...)
end
function interpolate_onto(
    Vnl::NonLocalPotential{S,Numerical,P}, Δr::Real, args...; kwargs...
)::NonLocalPotential{S,Numerical,P} where {S,P}
    return _apply(Vnl, interpolate_onto, Δr, args...; kwargs...)
end
function interpolate_onto(
    Vnl::NonLocalPotential{S,Numerical,P}, r::AbstractVector, args...; kwargs...
)::NonLocalPotential{S,Numerical,P} where {S,P}
    return _apply(Vnl, interpolate_onto, r, args...; kwargs...)
end
function truncate(
    Vnl::NonLocalPotential{S,Numerical,P}, args...; kwargs...
)::NonLocalPotential{S,Numerical,P} where {S,P}
    return _apply(Vnl, truncate, args...; kwargs...)
end

# function fht(
#     quantity::NonLocalPotential{RealSpace,P},
#     q::AbstractVector,
#     quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
#     interpolation_method::Interpolation.InterpolationMethod,
# ) where {P}
#     β = map(quantity.β) do βl
#         map(βl) do βln
#             fht(βln, q, quadrature_method, interpolation_method)
#         end
#     end
#     return NonLocalPotential{FourierSpace,_dual_constructor_of(P)}(β, quantity.D)
# end
# function fht(quantity::NonLocalPotential{RealSpace,P}) where {P}
#     β = map(quantity.β) do βl
#         map(fht, βl)
#     end
#     return NonLocalPotential{FourierSpace,_dual_constructor_of(P)}(β, quantity.D)
# end

# function ifht(
#     quantity::NonLocalPotential{FourierSpace,P},
#     r::AbstractVector,
#     quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
#     interpolation_method::Interpolation.InterpolationMethod,
# ) where {P}
#     β = map(quantity.β) do βl
#         map(βl) do βln
#             ifht(βln, r, quadrature_method, interpolation_method)
#         end
#     end
#     return NonLocalPotential{RealSpace,_dual_constructor_of(P)}(β, quantity.D)
# end
# function ifht(quantity::NonLocalPotential{FourierSpace,P}) where {P}
#     β = map(quantity.β) do βl
#         map(ifht, βl)
#     end
#     return NonLocalPotential{RealSpace,_dual_constructor_of(P)}(β, quantity.D)
# end
