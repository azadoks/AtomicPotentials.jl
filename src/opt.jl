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
                quadrature_method::NumericalQuadrature.QuadratureMethodOrType=NumericalQuadrature.Simpson,
                interpolation_method::Interpolation.InterpolationMethod=Interpolation.Spline(
                    4
                ),
            ) where {T<:Real}
                F = ht!(
                    weights_,
                    integrand_,
                    quantity.r,
                    quantity.f,
                    q,
                    angular_momentum(quantity),
                    quadrature_method,
                )
                interpolator = Interpolation.construct_interpolator(
                    q, F, interpolation_method
                )
                return $(AQ){FourierSpace,Numerical}(
                    q, F, interpolator, [getfield(quantity, arg) for arg in $(args)]...
                )
            end
            function iht!(
                weights_::AbstractVector,
                integrand_::AbstractVector{T},
                quantity::$(AQ){FourierSpace,Numerical},
                r::AbstractVector{T},
                quadrature_method::NumericalQuadrature.QuadratureMethodOrType=NumericalQuadrature.Simpson,
                interpolation_method::Interpolation.InterpolationMethod=Interpolation.Spline(
                    4
                ),
            ) where {T<:Real}
                f = iht!(
                    weights_,
                    integrand_,
                    quantity.r,
                    quantity.f,
                    r,
                    angular_momentum(quantity),
                    quadrature_method,
                )
                f .*= r .^ 2
                interpolator = Interpolation.construct_interpolator(
                    r, f, interpolation_method
                )
                return $(AQ){RealSpace,Numerical}(
                    r, f, interpolator, [getfield(quantity, arg) for arg in $(args)]...
                )
            end
        end,
    )
end

for (op!, op, Sin, Sout) in
    ((:ht!, :ht, :RealSpace, :FourierSpace), (:iht!, :iht, :FourierSpace, :RealSpace))
    eval(
        quote
            function $(op!)(
                weights_::AbstractVector,
                integrand_::AbstractVector{T},
                Vnl::NonLocalPotential{
                    $(Sin),Numerical,KleinmanBylanderProjector{$(Sin),Numerical}
                },
                q::AbstractVector{T},
                args...;
                kwargs...,
            )::NonLocalPotential{
                $(Sout),Numerical,KleinmanBylanderProjector{$(Sout),Numerical}
            } where {T<:Real}
                β = map(Vnl.β) do βl
                    map(βl) do βln
                        return $(op!)(weights_, integrand_, βln, q, args...; kwargs...)
                    end
                end
                return NonLocalPotential(β, Vnl.D)
            end
        end,
    )
    for (VDin, VDout) in (
        (:Nothing, :Nothing),
        (
            eval(
                quote
                    ChargeDensity{$(Sin),Numerical}
                end,
            ),
            eval(
                quote
                    ChargeDensity{$(Sout),Numerical}
                end,
            ),
        ),
    )
        for (CDin, CDout) in (
            (:Nothing, :Nothing),
            (
                eval(
                    quote
                        ChargeDensity{$(Sin),Numerical}
                    end,
                ),
                eval(
                    quote
                        ChargeDensity{$(Sout),Numerical}
                    end,
                ),
            ),
        )
            for (SPin, SPout) in (
                (:Nothing, :Nothing),
                (
                    eval(
                        quote
                            StateProjector{$(Sin),Numerical}
                        end,
                    ),
                    eval(
                        quote
                            StateProjector{$(Sout),Numerical}
                        end,
                    ),
                ),
            )
                for (AUGin, AUGout) in (
                    (:Nothing, :Nothing),
                    (
                        eval(
                            quote
                                Augmentation{$(Sin),AugmentationFunction{$(Sin),Numerical}}
                            end,
                        ),
                        eval(
                            quote
                                Augmentation{
                                    $(Sout),AugmentationFunction{$(Sout),Numerical}
                                }
                            end,
                        ),
                    ),
                )
                    eval(
                        quote
                            function $(op)(
                                potential::AtomicPotential{
                                    $(Sin),
                                    LocalPotential{$(Sin),Numerical},
                                    NonLocalPotential{
                                        $(Sin),
                                        Numerical,
                                        KleinmanBylanderProjector{$(Sin),Numerical},
                                    },
                                    $(VDin),
                                    $(CDin),
                                    $(SPin),
                                    $(AUGin),
                                },
                                q::AbstractVector{T},
                                args...;
                                kwargs...,
                            ) where {T<:Real}
                                # ::AtomicPotential{
                                #     $(Sout),
                                #     LocalPotential{$(Sout),Numerical},
                                #     NonLocalPotential{
                                #         $(Sout),
                                #         Numerical,
                                #         KleinmanBylanderProjector{$(Sout),Numerical},
                                #     },
                                #     $(VDout),
                                #     $(CDout),
                                #     $(SPout),
                                #     $(AUGout),
                                # }
                                n = max_r_length(potential)
                                weights_ = Vector{Float64}(undef, n)
                                integrand_ = Vector{T}(undef, n)

                                states = map(potential.states) do states_l
                                    map(states_l) do state_ln
                                        return $(op!)(
                                            weights_,
                                            integrand_,
                                            state_ln,
                                            q,
                                            args...;
                                            kwargs...,
                                        )
                                    end
                                end
                                return AtomicPotential{
                                    $(Sout),
                                    LocalPotential{$(Sout),Numerical},
                                    NonLocalPotential{
                                        $(Sout),
                                        Numerical,
                                        KleinmanBylanderProjector{$(Sout),Numerical},
                                    },
                                    $(VDout),
                                    $(CDout),
                                    $(SPout),
                                    $(AUGout),
                                }(
                                    potential.identifier,
                                    potential.symbol,
                                    $(op!)(
                                        weights_,
                                        integrand_,
                                        potential.local_potential,
                                        q,
                                        args...;
                                        kwargs...,
                                    ),
                                    $(op!)(
                                        weights_,
                                        integrand_,
                                        potential.nonlocal_potential,
                                        q,
                                        args...;
                                        kwargs...,
                                    ),
                                    $(op!)(
                                        weights_,
                                        integrand_,
                                        potential.valence_density,
                                        q,
                                        args...;
                                        kwargs...,
                                    ),
                                    $(op!)(
                                        weights_,
                                        integrand_,
                                        potential.core_density,
                                        q,
                                        args...;
                                        kwargs...,
                                    ),
                                    states,
                                    $(op!)(
                                        weights_,
                                        integrand_,
                                        potential.augmentation,
                                        q,
                                        args...;
                                        kwargs...,
                                    ),
                                )
                            end
                        end,
                    )
                end
            end
        end
    end
end
