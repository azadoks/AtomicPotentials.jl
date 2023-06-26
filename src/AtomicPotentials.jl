module AtomicPotentials

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
export interpolate_onto
include("atomic_quantity.jl")

## Fourier-Hankel transforms
export fht  # Fourier-Hankel transform
export ifht  # Inverse Fourier-Hankel transform
include("fourier_hankel_transform.jl")

## Interpolation
export interpolate_onto
include("interpolate_onto.jl")

## Truncation
include("truncate.jl")

## Local potentials
export AbstractLocalPotential
export LocalPotential
export HghLocalPotential
export CoulombLocalPotential
export GaussianLocalPotential
export CohenBergstresserLocalPotential
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
export AbstractNonLocalPotential
export NormConservingNonLocalPotential
export UltrasoftNonLocalPotential
include("non_local_potential.jl")

## Charge densities
export AbstractChargeDensity
# Valence charge densities
export AbstractValenceChargeDensity
export ValenceChargeDensity
export GaussianChargeDensity
# Core charge density
export AbstractCoreChargeDensity
export CoreChargeDensity
include("charge_density.jl")

## Atomic potential -- collection of atomic quantities
export AtomicPotential
export get_quantities
include("atomic_potential.jl")
include("pseudopotentialio.jl")

@setup_workload begin
    psp_files = []
    for psp_file_tuple in [
        ("pd_nc_sr_pbe_standard_0.4.1_psp8", "Si.psp8"),
        ("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf"),
        ("hgh_lda_hgh", "si-q4.hgh"),
    ]
        push!(psp_files, PseudoPotentialIO.load_psp_file(psp_file_tuple...))
    end
    @compile_workload begin
        for psp_file in psp_files
            pot = AtomicPotential(psp_file)
            pot_q = fht(pot, 0.0:1.0:10.0)
            pot′ = ifht(pot_q, 0.0:0.1:1.0)
        end
    end
end

end # module AtomicPotentials
