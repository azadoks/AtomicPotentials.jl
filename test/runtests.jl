using Test
using LinearAlgebra
using Random
using Aqua
using AtomicPotentials

Random.seed!(0)

TAGS = ARGS
isempty(TAGS) && (TAGS = ["all"])
TAGSREF = Ref(TAGS)

println("\tRunning tests (TAGS = $(join(TAGS, ", "))).")

# include("fixtures.jl")
@testset "AtomicPotentials.jl" begin
    # if any(in.(("all", "aqua"), TAGSREF))
    #     include("aqua.jl")
    # end

    if any(in.(("all", "numerical_quadrature"), TAGSREF))
        include("NumericalQuadrature/runtests.jl")
    end
    # include("atomic_potential.jl")
    # include("atomic_quantity.jl")
    # include("augmentation.jl")
    # include("charge_density.jl")
    # include("fast_sphericalbesselj.jl")
    # include("hankel_transform.jl")
    # include("interpolate_onto.jl")
    # include("local_potential.jl")
    # include("non_local_potential.jl")
    # include("projector.jl")
    # include("pseudopotentialio.jl")
    # include("truncate.jl")
end
