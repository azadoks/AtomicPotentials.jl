module AtomicPotentials
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

## Atomic potential -- collection of atomic quantities
export AtomicPotential
export get_quantities
include("atomic_potential.jl")

## Fourier-Hankel transforms
export fht  # Fourier-Hankel transform
export ifht  # Inverse Fourier-Hankel transform
include("fourier_hankel_transform.jl")

## Interpolation
export interpolate_onto
include("interpolation.jl")

## Local potentials
export AbstractLocalPotential
export LocalPotential
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
export NonLocalPotential
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
end # module AtomicPotentials
