for (AQ, args) in (
    (:KleinmanBylanderProjector, (:n, :l, :j)),
    (:ChargeDensity, ()),
    (:AugmentationFunction, (:n, :m, :l)),
)
    eval(
        quote
            function ht(
                quantity::$(AQ){RealSpace,Numerical},
                q::AbstractVector,
                quadrature_method::NumericalQuadrature.QuadratureMethodOrType=NumericalQuadrature.Simpson,
                interpolation_method::Interpolation.InterpolationMethod=Interpolation.Spline(
                    4
                ),
            )
                F = ht(
                    quantity.r, quantity.f, q, angular_momentum(quantity), quadrature_method
                )
                interpolator = Interpolation.construct_interpolator(
                    q, F, interpolation_method
                )
                return $(AQ){FourierSpace,Numerical}(
                    q, F, interpolator, [getfield(quantity, arg) for arg in $(args)]...
                )
            end
        end,
    )
end

for (op, Sin, Sout) in ((:ht, :RealSpace, :FourierSpace), (:iht, :FourierSpace, :RealSpace))
    eval(
        quote
            function $(op)(
                Vnl::NonLocalPotential{
                    $(Sin),Numerical,KleinmanBylanderProjector{$(Sin),Numerical}
                },
                args...;
                kwargs...,
            )::NonLocalPotential{
                $(Sout),Numerical,KleinmanBylanderProjector{$(Sout),Numerical}
            }
                β = map(Vnl.β) do βl
                    map(βl) do βln
                        return $(op)(βln, args...; kwargs...)
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
                                args...;
                                kwargs...,
                            )::AtomicPotential{
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
                            }
                                states = map(potential.states) do states_l
                                    map(states_l) do state_ln
                                        return $(op)(state_ln, args...; kwargs...)
                                    end
                                end
                                return AtomicPotential(
                                    potential.identifier,
                                    potential.symbol,
                                    $(op)(potential.local_potential, args...; kwargs...),
                                    $(op)(potential.nonlocal_potential, args...; kwargs...),
                                    $(op)(potential.valence_density, args...; kwargs...),
                                    $(op)(potential.core_density, args...; kwargs...),
                                    states,
                                    $(op)(potential.augmentation, args...; kwargs...),
                                )
                            end
                        end,
                    )
                end
            end
        end
    end
end
