@testset "KleinmanBylanderProjector" begin end

@testset "HghKleinmanBylanderProjector" begin
    r_test = 200:-0.1:0.0
    i_test = randperm(length(HGH_POTENTIALS))[1:30]
    qs = [0.01, 0.5, 2.5, 5.0, 10.0, 22.0]
    @testset "$(pot.identifier)" for pot_r in HGH_POTENTIALS[i_test]
        pot_q = ht(pot_r)
        non_local_r = pot_r.nonlocal_potential
        @testset "l=$(l)" for l in eachindex(nonlocal_r.β)
            @testset "n=$(n)" for n in eachindex(nonlocal_r.β[l])
                integrand(r, q) = 4π * r^2 * β_r(r) * sphericalbesselj(l, q * r)

                β_r = nonlocal_r.β[l][n]
                β_q = ht(β_r)
                pot_β_q = pot_q.nonlocal_potential.β[l][n]

                r_max = r_test[findfirst(r -> !isapprox(0, β_r(r)), r_test)]
                β_q_ref = map(qs) do q
                    quadgk(Base.Fix2(integrand, q), 0.0, rmax)[1]
                end

                @testset "q=$q" for (i, q) in enumerate(qs)
                    @test β_q(q) ≈ β_q_ref[i] rtol = 1e-12 atol = 1e-12
                    @test pot_β_q(q) ≈ β_q_ref[i] rtol = 1e-12 atol = 1e-12
                end
            end
        end
    end
end

@testset "StateProjector" begin end

@testset "HydrogenicProjector" begin end
