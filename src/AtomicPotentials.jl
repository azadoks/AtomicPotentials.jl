module AtomicPotentials

using BSplineKit
using LinearAlgebra
using OffsetArrays
using OrderedCollections
using Polynomials
using PrecompileTools
using Printf

using PeriodicTable: PeriodicTable
import PseudoPotentialIOExperimental as PseudoPotentialIO
import Bessels: gamma
import SpecialFunctions: erf

## Submodules
# Numerical integration / quadrature
export NumericalQuadrature
include("NumericalQuadrature/NumericalQuadrature.jl")
# Interpolation interface
export Interpolation
include("Interpolation/Interpolation.jl")

# Fast spherical Bessel functions of the first kind
include("fast_sphericalbesselj.jl")

## Flag structs, abstract types, fundamental functions
# Evaluation space flags
export EvaluationSpace
export RealSpace
export FourierSpace
# Analyticity flags
export Analyticity
export Analytical
export Numerical
# Abstract type and important associated functions
export AbstractAtomicQuantity
include("atomic_quantity.jl")

## Local potentials
export AbstractLocalPotential
export LocalPotential
export HghLocalPotential
export CoulombLocalPotential
export GaussianLocalPotential
export CohenBergstresserLocalPotential
export energy_correction
include("local_potential.jl")

## Projectors and projector-like quantitities
export AbstractProjector
# Kleinman-Bylander projectors
export AbstractKleinmanBylanderProjector
export KleinmanBylanderProjector
export HghKleinmanBylanderProjector
# Atomic state / pseudo-orbital / pseudo-wavefunction
export AbstractStateProjector
export StateProjector
export HydrogenicProjector
include("projector.jl")

## Augmentation functions
export AugmentationFunction
include("augmentation.jl")

## Non-local potentials
export NonLocalPotential
include("non_local_potential.jl")

## Charge densities
export AbstractChargeDensity
export ChargeDensity
export GaussianChargeDensity
include("charge_density.jl")

## Atomic potential -- collection of atomic quantities
export AtomicPotential
export charge_ionic
export charge_nuclear
export n_elec_valence
export n_elec_core
export count_n_proj_radial
export count_n_proj
export lmax
export angular_momenta
include("atomic_potential.jl")
include("pseudopotentialio.jl")

## Hankel transforms
export ht  # Hankel transform
export iht  # Inverse Hankel transform
include("hankel_transform.jl")

## Interpolation
export interpolate_onto
include("interpolate_onto.jl")

## Truncation
include("truncate.jl")

include("opt.jl")

@setup_workload begin
    psp_files = []
    for psp_file_tuple in [
        ("pd_nc_sr_pbe_standard_0.4.1_psp8", "B.psp8"),
        ("pd_nc_sr_pbe_standard_0.4.1_upf", "B.upf"),
        ("hgh_lda_hgh", "si-q4.hgh"),
    ]
        push!(psp_files, PseudoPotentialIO.load_psp_file(psp_file_tuple...))
    end
    @compile_workload begin
        for psp_file in psp_files
            pot = AtomicPotential(psp_file)
            pot_q = ht(pot, 0.0:1.0:10.0)
            pot′ = iht(pot_q, 0.0:0.1:1.0)
        end
    end
end

end # module AtomicPotentials
