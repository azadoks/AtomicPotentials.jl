import .NumericalQuadrature: QuadratureMethod, Trapezoid, integration_weights!
import .Interpolation: InterpolationMethod, Spline, construct_interpolator

function ht!(
    weights_::Vector{T},
    integrand_::Vector{TT},
    r::AbstractVector{T},
    r²f::AbstractVector{T},
    q::AbstractVector{TT},
    l::Int,
    quadrature_method::QuadratureMethod,
)::Vector{TT} where {T<:Real,TT<:Real}
    weights_ = @view weights_[eachindex(r)]
    integrand_ = @view integrand_[eachindex(r)]
    integration_weights!(weights_, r, quadrature_method)
    jl = fast_sphericalbesselj(l)
    F = map(q) do qi
        integrand_ .= r²f .* jl.(qi .* r)
        return 4π * dot(weights_, integrand_)
    end
    return F
end

ht!(::AbstractVector, ::AbstractVector, ::Nothing, args...; kwargs...) = nothing

function ht(
    r::AbstractVector{T},
    r²f::AbstractVector{T},
    q::AbstractVector{TT},
    l::Int,
    quadrature_method::QuadratureMethod,
)::AbstractVector{TT} where {T<:Real,TT<:Real}
    # The phase factor (-i)^l is not included
    @assert length(r) == length(r²f) "`r` and `r²f` must be the same length"

    weights_ = similar(r, T)
    integrand_ = similar(r, TT)
    F = similar(q, TT)

    integration_weights!(weights_, r, quadrature_method)

    map!(F, q) do qi
        integrand_ .= r²f .* fast_sphericalbesselj.(l, qi .* r)
        return 4π * dot(weights_, integrand_)
    end
    return F
end

function ht(
    r::AbstractVector, r²f, q::AbstractVector, l::Int, quadrature_method::QuadratureMethod
)
    return ht(r, r²f.(r), q, l, quadrature_method)
end

ht(::Nothing, args...; kwargs...) = nothing

# args... are included for interface consistency with numeric `ht`, which requires
# a q-point mesh and an integration method
function ht(quantity::AbstractAtomicQuantity{RealSpace,Analytical}, args...)
    return _construct_dual_quantity(quantity)
end

function ht(
    quantity::AbstractAtomicQuantity{RealSpace,Numerical},
    q::AbstractVector,
    quadrature_method::QuadratureMethod=Trapezoid(),
    interpolation_method::InterpolationMethod=Spline(4),
)::AbstractAtomicQuantity{FourierSpace,Numerical}
    F = ht(quantity.r, quantity.f, q, angular_momentum(quantity), quadrature_method)
    interpolator = construct_interpolator(q, F, interpolation_method)
    return _construct_dual_quantity(quantity; r=q, f=F, interpolator=interpolator)
end

function ht(
    Vloc::LocalPotential{RealSpace,Numerical},
    q::AbstractVector{T},
    quadrature_method::QuadratureMethod=Trapezoid(),
    interpolation_method::InterpolationMethod=Spline(4),
)::LocalPotential{FourierSpace,Numerical} where {T<:Real}
    weights_ = Vector{T}(undef, length(Vloc.r))
    integrand_ = Vector{T}(undef, length(Vloc.r))
    return ht!(weights_, integrand_, Vloc, q, quadrature_method, interpolation_method)
end

function ht!(
    weights_::AbstractVector,
    integrand_::AbstractVector{T},
    Vloc::LocalPotential{RealSpace,Numerical},
    q::AbstractVector{T},
    quadrature_method::QuadratureMethod=Trapezoid(),
    interpolation_method::InterpolationMethod=Spline(4),
)::LocalPotential{FourierSpace,Numerical} where {T<:Real}
    f = Vloc.r .* Vloc.f .- -Vloc.Z
    integration_weights!(weights_, Vloc.r, quadrature_method)
    F = map(q) do qi
        integrand_ .= f .* sin.(qi .* Vloc.r)
        integral = dot(weights_, integrand_)
        return 4T(π) * (integral / qi + -Vloc.Z / qi .^ 2)
    end
    interpolator = construct_interpolator(
        q[(begin + 1):end], F[(begin + 1):end], interpolation_method
    )
    return LocalPotential{FourierSpace,Numerical}(q, F, interpolator, Vloc.Z)

    # r²f = Vloc.r .* (Vloc.r .* Vloc.f .- -Vloc.Z)  # == r² (Vloc - -Z/r)
    # F =
    #     ht!(
    #         weights_, integrand_, Vloc.r, r²f, q, angular_momentum(Vloc), quadrature_method
    #     ) .+ 4π .* (-Vloc.Z ./ q .^ 2)
    # interpolator = construct_interpolator(
    #     q[(begin + 1):end], F[(begin + 1):end], interpolation_method
    # )
    # return LocalPotential{FourierSpace,Numerical}(q, F, interpolator, Vloc.Z)
end

function iht(
    quantity::AbstractAtomicQuantity{FourierSpace,Numerical},
    r::AbstractVector,
    quadrature_method::QuadratureMethod=Trapezoid(),
    interpolation_method::InterpolationMethod=Spline(4),
)::AbstractAtomicQuantity{RealSpace,Numerical}
    f = iht(quantity.r, quantity.f, r, angular_momentum(quantity), quadrature_method)
    r²f = r .^ 2 .* f
    interpolator = construct_interpolator(r, r²f, interpolation_method)
    return _construct_dual_quantity(quantity; r=r, f=r²f, interpolator=interpolator)
end

function iht(
    q::AbstractVector,
    F::AbstractVector,
    r::AbstractVector,
    l::Int,
    quadrature_method::QuadratureMethod,
)
    # The phase factor (i)^l is not included.
    # The inverse Hankel transform is the same as the direct transform
    # except for a normalization factor (1/(2π)³ and the sign of the phase factor
    # (i)ˡ vs. (-i)ˡ
    f = ht(q, q .^ 2 .* F, r, l, quadrature_method) ./ (2π)^3
    return f
end

function iht!(
    weights_::AbstractVector,
    integrand_::AbstractVector,
    q::AbstractVector,
    F::AbstractVector,
    r::AbstractVector,
    l::Int,
    quadrature_method::QuadratureMethod,
)
    # The phase factor (i)^l is not included.
    # The inverse Hankel transform is the same as the direct transform
    # except for a normalization factor (1/(2π)³ and the sign of the phase factor
    # (i)ˡ vs. (-i)ˡ
    f = ht!(weights_, integrand_, q, q .^ 2 .* F, r, l, quadrature_method) ./ (2π)^3
    return f
end

iht!(::AbstractVector, ::AbstractVector, ::Nothing, args...; kwargs...) = nothing

iht(::Nothing, args...; kwargs...) = nothing

# args... are included for interface consistency with numeric `iht`, which requires
# an r-point mesh and an integration method
function iht(
    quantity::AbstractAtomicQuantity{FourierSpace,Analytical}, args...
)::AbstractAtomicQuantity{RealSpace,Analytical}
    return _construct_dual_quantity(quantity)
end

function iht(
    q::AbstractVector, F, r::AbstractVector, l::Int, quadrature_method::QuadratureMethod
)
    return iht(q, F.(q), r, l, quadrature_method)
end

function iht(
    Vloc::LocalPotential{FourierSpace},
    r::AbstractVector,
    quadrature_method::QuadratureMethod=Trapezoid(),
    interpolation_method::InterpolationMethod=Spline(4),
)::LocalPotential{RealSpace}
    q²F = Vloc.r .^ 2 .* Vloc.f .- -Vloc.Z  # == q² (Vloc - -Z/q²)
    f =
        ht(Vloc.r, q²F, r, angular_momentum(Vloc), quadrature_method) .+
        4π / (2π)^3 .* (-Vloc.Z ./ r)
    interpolator = construct_interpolator(
        r[(begin + 1):end], f[(begin + 1):end], interpolation_method
    )
    return LocalPotential{RealSpace,Numerical}(r, f, interpolator, Vloc.Z)
end

function iht!(
    weights_::AbstractVector,
    integrand_::AbstractVector,
    Vloc::LocalPotential{FourierSpace},
    r::AbstractVector,
    quadrature_method::QuadratureMethod=Trapezoid(),
    interpolation_method::InterpolationMethod=Spline(4),
)::LocalPotential{RealSpace}
    q²F = Vloc.r .^ 2 .* Vloc.f .- -Vloc.Z  # == q² (Vloc - -Z/q²)
    f =
        ht!(
            weights_, integrand_, Vloc.r, q²F, r, angular_momentum(Vloc), quadrature_method
        ) .+ 4π / (2π)^3 .* (-Vloc.Z ./ r)
    interpolator = construct_interpolator(
        r[(begin + 1):end], f[(begin + 1):end], interpolation_method
    )
    return LocalPotential{RealSpace,Numerical}(r, f, interpolator, Vloc.Z)
end

for (op, S) in ((:ht, :RealSpace), (:iht, :FourierSpace))
    eval(
        quote
            function $(op)(
                x::Union{NonlocalPotential{$(S)},Augmentation{$(S)},AtomicPotential{$(S)}},
                args...;
                kwargs...,
            )
                return _apply(x, $(op), args...; kwargs...)
            end
        end,
    )
end
