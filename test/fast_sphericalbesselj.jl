import AtomicPotentials: fast_sphericalbesselj
import Bessels: sphericalbesselj

@testset "fast_sphericalbesselj" begin
    for l in 0:5
        for x in (rand(100) .* 100)
            @test sphericalbesselj(l, x) ≈ fast_sphericalbesselj(l)(x) atol = 1e-8
        end
    end
    @test_throws "6 not supported" fast_sphericalbesselj(6)
end
