import .NumericalQuadrature: QuadratureMethod, Trapezoid
import .Interpolation: InterpolationMethod, Spline, construct_interpolator

for (AQ, args) in (
    (:KleinmanBylanderProjector, (:n, :l, :j)),
    (:ChargeDensity, ()),
    (:AugmentationFunction, (:n, :m, :l)),
    (:StateProjector, (:n, :l, :j)),
)
    eval(
        quote
            function ht!(
                weights_::AbstractVector,
                integrand_::AbstractVector{T},
                quantity::$(AQ){RealSpace,Numerical},
                q::AbstractVector{T},
                quadrature_method::QuadratureMethod=Trapezoid(),
                interpolation_method::InterpolationMethod=Spline(4),
            ) where {T<:Real}
                r = quantity.r
                f = quantity.f
                l = angular_momentum(quantity)
                F = ht!(weights_, integrand_, r, f, q, l, quadrature_method)
                interpolator = construct_interpolator(q, F, interpolation_method)
                other_args = [getfield(quantity, arg) for arg in $(args)]
                return $(AQ){FourierSpace,Numerical}(q, F, interpolator, other_args...)
            end
            function ht(
                quantity::$(AQ){RealSpace,Numerical},
                q::AbstractVector{T},
                quadrature_method::QuadratureMethod=Trapezoid(),
                interpolation_method::InterpolationMethod=Spline(4),
            ) where {T<:Real}
                n = max_r_length(quantity)
                #* this assumes that all Numerical quantities have an element type
                #* of Float64
                weights_ = Vector{Float64}(undef, n)
                integrand_ = Vector{T}(undef, n)
                return ht!(
                    weights_,
                    integrand_,
                    quantity,
                    q,
                    quadrature_method,
                    interpolation_method,
                )
            end
            function iht!(
                weights_::AbstractVector,
                integrand_::AbstractVector{T},
                quantity::$(AQ){FourierSpace,Numerical},
                r::AbstractVector{T},
                quadrature_method::QuadratureMethod=Trapezoid(),
                interpolation_method::InterpolationMethod=Spline(4),
            ) where {T<:Real}
                q = quantity.r
                F = quantity.f
                l = angular_momentum(quantity)
                f = iht!(weights_, integrand_, q, F, r, l, quadrature_method)
                f .*= r .^ 2
                interpolator = construct_interpolator(r, f, interpolation_method)
                other_args = [getfield(quantity, arg) for arg in $(args)]
                return $(AQ){RealSpace,Numerical}(r, f, interpolator, other_args...)
            end
            function iht(
                quantity::$(AQ){RealSpace,Numerical},
                q::AbstractVector{T},
                quadrature_method::QuadratureMethod=Trapezoid(),
                interpolation_method::InterpolationMethod=Spline(4),
            ) where {T<:Real}
                n = max_r_length(quantity)
                #* this assumes that all Numerical quantities have an element type
                #* of Float64
                weights_ = Vector{Float64}(undef, n)
                integrand_ = Vector{T}(undef, n)
                return iht!(
                    weights_,
                    integrand_,
                    quantity,
                    q,
                    quadrature_method,
                    interpolation_method,
                )
            end
        end,
    )
end

for (op!, Sin, Sout) in
    ((:ht!, :RealSpace, :FourierSpace), (:iht!, :FourierSpace, :RealSpace))
    #! format: off
    AFin = eval(quote AugmentationFunction{$(Sin),Numerical} end)
    AFout = eval(quote AugmentationFunction{$(Sout),Numerical} end)
    AUGin = eval(quote Augmentation{$(Sin),$(AFin)} end)
    AUGout = eval(quote Augmentation{$(Sout),$(AFout)} end)
    eval(quote
        function $(op!)(
            weights_::AbstractVector,
            integrand_::AbstractVector{T},
            aug::$(AUGin),
            q::AbstractVector{T},
            args...;
            kwargs...,
        ) where {T<:Real}
            Q = map(aug.Q) do Ql
                map(Ql) do Qlnm
                    return $(op!)(weights_, integrand_, Qlnm, q, args...; kwargs...)
                end
            end
            return $(AUGout)(Q, aug.q)
        end
    end)
    #! format: on
end

for (op!, Sin, Sout) in
    ((:ht!, :RealSpace, :FourierSpace), (:iht!, :FourierSpace, :RealSpace))
    #! format: off
    KBin = eval(quote KleinmanBylanderProjector{$(Sin),Numerical} end)
    KBout = eval(quote KleinmanBylanderProjector{$(Sout),Numerical} end)
    AFin = eval(quote AugmentationFunction{$(Sin),Numerical} end)
    AFout = eval(quote AugmentationFunction{$(Sout),Numerical} end)

    nothing_io = (:Nothing, :Nothing)

    aug_io = (
        eval(quote Augmentation{$(Sin),$(AFin)} end),
        eval(quote Augmentation{$(Sout),$(AFout)} end),
    )

    for (AUGin, AUGout) in (nothing_io, aug_io)
        NLin = eval(quote NonlocalPotential{$(Sin),Numerical,$(KBin),$(AUGin)} end)
        NLout = eval(quote NonlocalPotential{$(Sout),Numerical,$(KBout),$(AUGout)} end)
        eval(quote
            function $(op!)(
                weights_::AbstractVector,
                integrand_::AbstractVector{T},
                Vnl::$(NLin),
                q::AbstractVector{T},
                args...;
                kwargs...,
            ) where {T<:Real}
                β = map(Vnl.β) do βl
                    map(βl) do βln
                        return $(op!)(weights_, integrand_, βln, q, args...; kwargs...)
                    end
                end
                aug = $(op!)(weights_, integrand_, Vnl.augmentation, q, args...; kwargs...)
                return $(NLout)(β, Vnl.D, aug)
            end
        end)
    end
    #! format: on
end

for (op!, op, Sin, Sout) in
    ((:ht!, :ht, :RealSpace, :FourierSpace), (:iht!, :iht, :FourierSpace, :RealSpace))
#! format: off
    Lin = eval(quote LocalPotential{$(Sin),Numerical} end)
    Lout = eval(quote LocalPotential{$(Sout),Numerical} end)
    KBin = eval(quote KleinmanBylanderProjector{$(Sin),Numerical} end)
    KBout = eval(quote KleinmanBylanderProjector{$(Sout),Numerical} end)
    AFin = eval(quote AugmentationFunction{$(Sin),Numerical} end)
    AFout = eval(quote AugmentationFunction{$(Sout),Numerical} end)

    nothing_io = (:Nothing, :Nothing)

    charge_io = (
        eval(quote ChargeDensity{$(Sin),Numerical} end),
        eval(quote ChargeDensity{$(Sout),Numerical} end),
    )
    state_io = (
        eval(quote StateProjector{$(Sin),Numerical} end),
        eval(quote StateProjector{$(Sout),Numerical} end),
    )

    aug_io = (
        eval(quote Augmentation{$(Sin),$(AFin)} end),
        eval(quote Augmentation{$(Sout),$(AFout)} end),
    )

    for (VDin, VDout) in (nothing_io, charge_io), (CDin, CDout) in (nothing_io, charge_io)
        for (SPin, SPout) in (nothing_io, state_io), (AUGin, AUGout) in (nothing_io, aug_io)
        #! format: off
        NLin = eval(quote NonlocalPotential{$(Sin),Numerical,$(KBin),$(AUGin)} end)
        NLout = eval(quote NonlocalPotential{$(Sout),Numerical,$(KBout),$(AUGout)} end)
        eval(quote
            function $(op)(
                pot::AtomicPotential{$(Sin),$(Lin),$(NLin),$(VDin),$(CDin),$(SPin)},
                q::AbstractVector{T},
                args...;
                kwargs...,
            ) where {T<:Real}
                n = max_r_length(pot)
                #* this assumes that all Numerical quantities have an element type
                #* of Float64
                weights_ = Vector{Float64}(undef, n)
                integrand_ = Vector{T}(undef, n)

                states = map(pot.states) do states_l
                    map(states_l) do state_ln
                        return $(op!)(weights_, integrand_, state_ln, q, args...; kwargs...)
                    end
                end
                return AtomicPotential{
                    $(Sout),$(Lout),$(NLout),$(VDout),$(CDout),$(SPout)
                }(
                    pot.identifier,
                    pot.symbol,
                    $(op!)(weights_, integrand_, pot.local_potential, q, args...; kwargs...),
                    $(op!)(weights_, integrand_, pot.nonlocal_potential, q, args...; kwargs...),
                    $(op!)(weights_, integrand_, pot.valence_density, q, args...; kwargs...),
                    $(op!)(weights_, integrand_, pot.core_density, q, args...; kwargs...),
                    states
                )
            end
        end)
        #! format: on
        end
    end
end
