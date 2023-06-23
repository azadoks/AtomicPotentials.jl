# args... are included for interface consistency with numeric `fht`, which requires
# a q-point mesh and an integration method
function fht(
    quantity::AbstractAtomicQuantity{RealSpace,Analytical}, args...
)::AbstractAtomicQuantity{FourierSpace,Analytical}
    return construct_dual_quantity(quantity)
end
function fht(
    quantity::AbstractAtomicQuantity{RealSpace,Numerical},
    q::AbstractVector,
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    F = fht(quantity.r, quantity.f, q, angular_momentum(quantity), method)
    return construct_dual_quantity(quantity; r=q, f=F)
end
function fht(
    r::AbstractVector,
    r²f::AbstractVector,
    q::AbstractVector,
    l::Int,
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    # The phase factor (-i)^l is not included
    @assert length(r) == length(r²f) "`r` and `r²f` must be the same length"

    weights_ = zeros(r)
    integrand_ = zeros(r)

    integration_weights!(weights_, r, method)
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
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    r²f_ = zeros(r)
    r²f_ .= r²f.(r)
    return fht(r, r²f_, q, l, method)
end

# args... are included for interface consistency with numeric `ifht`, which requires
# an r-point mesh and an integration method
function ifht(
    ::AbstractAtomicQuantity{FourierSpace,Analytical}, args...
)::AbstractAtomicQuantity{RealSpace,Analytical}
    return construct_dual_quantity(quantity)
end
function ifht(
    quantity::AbstractAtomicQuantity{FourierSpace,Numerical},
    r::AbstractVector,
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    f = ifht(
        quantity.r, quantity.r .^ 2 .* quantity.f, r, angular_momentum(quantity), method
    )
    r²f = r .^ 2 .* f
    return construct_dual_quantity(quantity; r=r, f=r²f)
end
function ifht(
    q::AbstractVector,
    F::AbstractVector,
    r::AbstractVector,
    l::Int,
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    # The phase factor (i)^l is not included.
    # The inverse Fourier-Hankel transform is the same as the direct transform
    # except for a normalization factor (1/(2π)³ and the sign of the phase factor
    # (i)ˡ vs. (-i)ˡ
    q²F = q .^ 2 .* F
    f = fht(q, q²F, r, l, method) ./ (2π)^3
    return f
end
function ifht(
    q::AbstractVector,
    F,
    r::AbstractVector,
    l::Int,
    method::NumericalQuadrature.QuadratureMethodOrType,
)
    F_ = zeros(q)
    F_ .= F.(q)
    return ifht(q, F_, r, l, method)
end
