module AtomicPotentials

using LinearAlgebra
using Polynomials
using PrecompileTools
using Printf

using PeriodicTable: PeriodicTable
import PseudoPotentialIOExperimental as PseudoPotentialIO
import Bessels: gamma
import SpecialFunctions: erf

## Submodules
# Classification of radial meshes
export MeshClassification
include("MeshClassification/MeshClassification.jl")
# Numerical integration / quadrature
export NumericalQuadrature
include("NumericalQuadrature/NumericalQuadrature.jl")
# Interpolation interface
export Interpolation
include("Interpolation/Interpolation.jl")

## Common functions
# Fast spherical Bessel functions of the first kind
include("fast_sphericalbesselj.jl")
export fast_sphericalbesselj
# Radial Fourier transform
include("radial_fourier_transform.jl")
export rft
export irft

## Generic atomic quantities
# Evaluation space singletons
export EvaluationSpace
export RealSpace
export FourierSpace
# Analyticity singletons
export Analyticity
export Analytical
export Numerical
# Abstract type and interface
export AbstractQuantity
export angular_momentum
export radial_grid
export radial_function
export n_x_factors
export evaluate
# TODO: maybe these two should be provided by Interpolation?
export interpolate
export resample
# Numerical atomic quantity
export NumericalQuantity
# Analytical atomic quantities
export HghProjector
export HydrogenicProjector
include("quantity.jl")

## Local potentials
# Abstract type and interface
export AbstractLocalPotential
export ionic_charge
export energy_correction
# Concrete types
export NumericalLocalPotential
export HghLocalPotential
export CoulombLocalPotential
export GaussianLocalPotential
export CohenBergstresserLocalPotential
# Local potential correction
export AbstractLocalPotentialCorrection
export CoulombLocalPotentialCorrection
export CoulombErfLocalPotentialCorrection
include("local_potential.jl")

## Container types
export Augmentation
export NonLocalPotential
export AtomicPotential
include("container.jl")

include("pseudopotentialio.jl")
end
