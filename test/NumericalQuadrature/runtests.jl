using AtomicPotentials.NumericalQuadrature
import AtomicPotentials.NumericalQuadrature:
    is_uniform, compute_log_with_zero, compute_log_without_zero
import SpecialFunctions: erfi

@testset "NumericalQuadrature" begin
    function sin_info(x)
        y = sin.(x)
        I_true = cos(first(x)) - cos(last(x))
        return (; y, I_true)
    end

    function exp_sin_info(x)
        xn = last(x)
        y = exp.(-(x .^ 2)) .* sin.(x)
        I_true = real(
            -(√π * (-2erfi(1 / 2) + erfi(1 / 2 - im * xn) + erfi(1 / 2 + im * xn))) /
            (4 * exp(1 / 4)),
        )
        return (; y, I_true)
    end

    log_x1 = log_b = 1e-6
    lin_x1 = 0.0
    quadrature_methods = (Trapezoid(), Simpson(), QESimpson(), AbinitCorrectedTrapezoid())
    @testset "$(quadrature_method)" for quadrature_method in quadrature_methods
        @testset "$(test_case)" for test_case in (sin_info, exp_sin_info)
            @testset "∫(0,$(xn))[f(x) dx]" for xn in (Float64(2π), Float64(5π / 2))
                @testset "even no. intervals $(isodd(n))" for n in (10001, 10002)
                    log_a = log(xn / log_x1) / (n - 1)
                    @testset "$(mesh_type)" for (mesh_type, x, atol) in (
                        ("linear", range(lin_x1, xn, n), 1e-7),
                        ("log_with_zero", compute_log_with_zero.(log_a, log_b, 1:n), 1e-5),
                        (
                            "log_without_zero",
                            compute_log_without_zero.(log_a, log_b, 1:n),
                            1e-5,
                        ),
                    )
                        # the QE method is _very_ bad in this case
                        if isa(quadrature_method, QESimpson) && !isodd(n)
                            atol = 1e-1
                        end
                        # the ABINIT method doesn't support non-uniform meshes
                        if isa(quadrature_method, AbinitCorrectedTrapezoid) &&
                            !is_uniform(x)
                            continue
                        end
                        y, I_true = test_case(x)
                        weights = integration_weights(x, quadrature_method)
                        @test isapprox(dot(weights, y), I_true; atol)
                    end
                end
            end
        end
    end
end
