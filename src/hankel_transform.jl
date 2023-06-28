ht(::Nothing, args...; kwargs...) = nothing

# args... are included for interface consistency with numeric `ht`, which requires
# a q-point mesh and an integration method
function ht(quantity::AbstractAtomicQuantity{RealSpace,Analytical}, args...)
    return _construct_dual_quantity(quantity)
end

function ht(
    quantity::AbstractAtomicQuantity{RealSpace,Numerical},
    q::AbstractVector,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType=NumericalQuadrature.Simpson,
    interpolation_method::Interpolation.InterpolationMethod=Interpolation.Spline(4),
)
    F = ht(quantity.r, quantity.f, q, angular_momentum(quantity), quadrature_method)
    interpolator = Interpolation.construct_interpolator(q, F, interpolation_method)
    return _construct_dual_quantity(quantity; r=q, f=F, interpolator=interpolator)
end

function ht(
    r::AbstractVector,
    r²f::AbstractVector,
    q::AbstractVector,
    l::Int,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
)
    # The phase factor (-i)^l is not included
    @assert length(r) == length(r²f) "`r` and `r²f` must be the same length"

    weights_ = zero(r)
    integrand_ = zero(r)

    NumericalQuadrature.integration_weights!(weights_, r, quadrature_method)
    jₗ = fast_sphericalbesselj(l)

    F = map(q) do qi
        integrand_ .= r²f .* jₗ.(qi .* r)
        return 4π * dot(weights_, integrand_)
    end
    return F
end

function ht(
    r::AbstractVector,
    r²f,
    q::AbstractVector,
    l::Int,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
)
    return ht(r, r²f.(r), q, l, quadrature_method)
end

function ht(
    x::Union{
        NonLocalPotential{RealSpace},Augmentation{RealSpace},AtomicPotential{RealSpace}
    },
    args...;
    kwargs...,
)
    return _apply(x, ht, args...; kwargs...)
end

function ht(
    ::AugmentationFunction{RealSpace,Numerical}, q::AbstractVector, args...; kwargs...
)
    return error("`ht` not implemented for `AugmentationFunction`")
end

iht(::Nothing, args...; kwargs...) = nothing

# args... are included for interface consistency with numeric `iht`, which requires
# an r-point mesh and an integration method
function iht(quantity::AbstractAtomicQuantity{FourierSpace,Analytical}, args...)
    return _construct_dual_quantity(quantity)
end

function iht(
    quantity::AbstractAtomicQuantity{FourierSpace,Numerical},
    r::AbstractVector,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType=NumericalQuadrature.Simpson,
    interpolation_method::Interpolation.InterpolationMethod=Interpolation.Spline(4),
)
    f = iht(quantity.r, quantity.f, r, angular_momentum(quantity), quadrature_method)
    r²f = r .^ 2 .* f
    interpolator = Interpolation.construct_interpolator(r, r²f, interpolation_method)
    return _construct_dual_quantity(quantity; r=r, f=r²f, interpolator=interpolator)
end

function iht(
    q::AbstractVector,
    F::AbstractVector,
    r::AbstractVector,
    l::Int,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
)
    # The phase factor (i)^l is not included.
    # The inverse Fourier-Hankel transform is the same as the direct transform
    # except for a normalization factor (1/(2π)³ and the sign of the phase factor
    # (i)ˡ vs. (-i)ˡ
    f = ht(q, q .^ 2 .* F, r, l, quadrature_method) ./ (2π)^3
    return f
end

function iht(
    q::AbstractVector,
    F,
    r::AbstractVector,
    l::Int,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
)
    return iht(q, F.(q), r, l, quadrature_method)
end

function iht(
    x::Union{
        NonLocalPotential{FourierSpace},
        Augmentation{FourierSpace},
        AtomicPotential{FourierSpace},
    },
    args...;
    kwargs...,
)
    return _apply(x, iht, args...; kwargs...)
end

function iht(
    ::AugmentationFunction{FourierSpace,Numerical}, r::AbstractVector, args...; kwargs...
)
    return error("`iht` not implemented for `AugmentationFunction`")
end