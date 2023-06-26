using LinearAlgebra

fht(::Nothing, args...; kwargs...) = nothing
# args... are included for interface consistency with numeric `fht`, which requires
# a q-point mesh and an integration method
function fht(quantity::AbstractAtomicQuantity{RealSpace,Analytical}, args...)
    return _construct_dual_quantity(quantity)
end
function fht(
    quantity::AbstractAtomicQuantity{RealSpace,Numerical},
    q::AbstractVector,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType=NumericalQuadrature.Simpson,
    interpolation_method::Interpolation.InterpolationMethod=Interpolation.Spline(4),
)
    F = fht(quantity.r, quantity.f, q, angular_momentum(quantity), quadrature_method)
    interpolator = Interpolation.construct_interpolator(q, F, interpolation_method)
    return _construct_dual_quantity(quantity; r=q, f=F, interpolator=interpolator)
end
function fht(
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
function fht(
    r::AbstractVector,
    r²f,
    q::AbstractVector,
    l::Int,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
)
    # r²f_ = zero(r)
    # r²f_ .= r²f.(r)
    return fht(r, r²f.(r), q, l, quadrature_method)
end

ifht(::Nothing, args...; kwargs...) = nothing
# args... are included for interface consistency with numeric `ifht`, which requires
# an r-point mesh and an integration method
function ifht(quantity::AbstractAtomicQuantity{FourierSpace,Analytical}, args...)
    return _construct_dual_quantity(quantity)
end
function ifht(
    quantity::AbstractAtomicQuantity{FourierSpace,Numerical},
    r::AbstractVector,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType=NumericalQuadrature.Simpson,
    interpolation_method::Interpolation.InterpolationMethod=Interpolation.Spline(4),
)
    f = ifht(quantity.r, quantity.f, r, angular_momentum(quantity), quadrature_method)
    r²f = r .^ 2 .* f
    interpolator = Interpolation.construct_interpolator(r, r²f, interpolation_method)
    return _construct_dual_quantity(quantity; r=r, f=r²f, interpolator=interpolator)
end
function ifht(
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
    # q²F = q .^ 2 .* F
    f = fht(q, q .^ 2 .* F, r, l, quadrature_method) ./ (2π)^3
    return f
end
function ifht(
    q::AbstractVector,
    F,
    r::AbstractVector,
    l::Int,
    quadrature_method::NumericalQuadrature.QuadratureMethodOrType,
)
    # F_ = zero(q)
    # F_ .= F.(q)
    return ifht(q, F.(q), r, l, quadrature_method)
end
