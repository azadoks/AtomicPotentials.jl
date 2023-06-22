function fht(quantity::AbstractAtomicQuantity{RealSpace,Analytical})::AbstractAtomicQuantity{FourierSpace,Analytical}
    return construct_dual_quantity(quantity)
end
function fht(quantity::AbstractAtomicQuantity{RealSpace,Numerical}, q::AbstractVector)
    F = fht(quantity.r, quantity.f, q, angular_momentum(quantity))
    return construct_dual_quantity(quantity; r=q, f=F)
end
function fht(r::AbstractVector, r²f::AbstractVector, q::AbstractVector, l::Int)
    # The phase factor (-i)^l is not included
    @assert length(r) == length(r²f) "`r` and `r²f` must be the same length"

    weights_ = zeros(r)
    q_dot_r_ = zeros(r)
    integrand_ = zeros(r)

    integration_weights!(weights_, r, integration_method)
    jₗ = fast_sphericalbesselj(l)

    F = map(q) do qi
        q_dot_r_ .= qi .* r
        integrand_ .= r²f .* jₗ.(q_dot_r_)
        return 4π * dot(weights_, integrand_)
    end
    return F
end
function fht(r::AbstractVector, r²f, q::AbstractVector, l::Int)
    r²f_ = zeros(r)
    r²f_ .= r²f.(r)
    return fht(r, r²f_, q, l)
end

function ifht(::AbstractAtomicQuantity{FourierSpace,Analytical})::AbstractAtomicQuantity{RealSpace,Analytical}
    return construct_dual_quantity(quantity)
end
function ifht(quantity::AbstractAtomicQuantity{FourierSpace,Numerical}, r::AbstractVector)
    f = ifht(quantity.r, quantity.f, r, angular_momentum(quantity))
    r²f = r.^2 .* f
    return construct_dual_quantity(quantity; r=r, f=r²f)
end
function ifht(q::AbstractVector, q²F::AbstractVector, r::AbstractVector, l::Int)
    # The phase factor (i)^l is not included.
    # The inverse Fourier-Hankel transform is the same as the direct transform
    # except for a normalization factor (1/(2π)³ and the sign of the phase factor
    # (i)ˡ vs. (-i)ˡ
    f = fht(q, q²F, r, l) ./ (2π)^3
    return f
end
function ifht(q::AbstractVector, q²F, r::AbstractVector, l::Int)
    q²F_ = zeros(r)
    q²F_ .= q²F.(q)
    return ifht(q, q²F_, r, l)
end
