@doc raw"""
Composite Simpson's rule quadrature implemented as described in the
[Wikipedia article](https://en.wikipedia.org/wiki/Simpson%27s_rule).

## Non-uniform grid

```math
\int_a^b f(x) dx \approx
\sum_{i=0}{N/2-1} \frac{\Delta x_{2i} + \Delta x_{2i+1}}{6} \left[
    \left( 2 - \frac{\Delta x_{2i+1}}{\Delta x_{2i}} \right) f_{2i} +
    \frac{
        \left( \Delta x_{2i} + \Delta x_{2i+1} \right)^2
        }{
        \Delta x_{2i} \Delta x_{2i+1}
    } f_{2i+1} +
    \left( 2 - \frac{\Delta x_{2i}}{\Delta x_{2i + 1}} \right) f_{2i + 2}
\right]
```
where
```math
f_{k} = f \left( a + \sum_{i=0}^{k-1} \Delta x_{i} \right)
```

In the case of an odd number of subintervals, the above formulae are used up to the second
to last interval, and the last interval is computed seperately as
```math
\frac{
    2 \Delta x_{N-1}^2 + 3 \Delta x_{N-1} \Delta x_{N-2}
    }{
    6\left( \Delta x_{N-1} + \Delta x{N-1} \right)
} f_{N} +
\frac{
    \Delta x_{N-1}^2 + 3 \Delta x_{N-1} \Delta x_{N-1}
    }{
    6 \Delta x_{N-1}
} f_{N-1} +
\frac{
    \Delta x_{N-1}^3
    }{
    6 \Delta x_{N-2} \left( \Delta x_{N-2} + \Delta x_{N-1} \right)
} f_{N-2}
```

## Uniform grid

```math
\int_a^b f(x) dx \approx
\frac{1}{3} \Delta x \left[
    f(x_0) + 4 \sum_{i=0}{N/2} f(x_{2i-1}) + 2 \sum_{i=0}{N/2-1} f(x_{2i}) + f(x_N)
\right]
```

In the case of an odd number of subintervals, the above formula is used for the first N-1
subintervals, and the last subinterval is handled with the Trapezoid method. This approach
is generally more well behaved in the specific case of pseudopotential quantities because
the function is generally close to zero in the last interval and the error made by
the Trapezoid rule w.r.t. Simpson's rule has a smaller effect on the value of the integral.
This approach is also equivalent to the approach taken in the non-uniform case.
"""
struct Simpson <: QuadratureMethod end

function integration_weights!(weights::AbstractVector, x::AbstractVector, method::Simpson)
    length(weights) == length(x) ||
        throw(DimensionMismatch("weights and x arrays must have a common length"))
    length(x) >= 4 || throw(ArgumentError("x must have a length of at least 4"))
    mesh = Mesh(x)
    if typeof(mesh) <: Uniform
        return integration_weights_uniform!(weights, mesh.dx, method)
    else
        return integration_weights_nonuniform!(weights, x, method)
    end
end

function integration_weights_uniform!(weights::AbstractVector, Δx::Real, ::Simpson)
    N = length(weights) - 1  # Number of intervals
    #! format: off
    if !isodd(N)  # Standard Simpson's composite 1/3 rule
        weights[begin                  ]  = 1 / 3 * Δx
        weights[(begin + 1):2:(end - 1)] .= 4 * 1 / 3 * Δx
        weights[(begin + 2):2:(end - 1)] .= 2 * 1 / 3 * Δx
        weights[end                    ]  = 1 / 3 * Δx
    else  # * For functions which decay to zero as r -> ∞, this is a better approximation
        # If the number of intervals is odd, apply Simpsons method to the first N-1 intervals
        # and the Trapezoidal method to the last interval.
        weights[begin                  ]  = 1 / 3 * Δx
        weights[(begin + 1):2:(end - 2)] .= 4 * 1 / 3 * Δx
        weights[(begin + 2):2:(end - 2)] .= 2 * 1 / 3 * Δx
        weights[end - 1                ]  = 5 / 6 * Δx
        weights[end                    ]  = 1 / 2 * Δx
    end
    #! format: on
    # else
    #     # If the number of intervals is odd, apply Simpsons method to the last N-1 intervals
    #     # and the Trapezoidal method to the first interval.
    #     weights[begin] = 1 / 2 * Δx
    #     weights[begin+1] = 5 / 6 * Δx
    #     weights[(begin+2):2:(end-1)] .= 4 * 1 / 3 * Δx
    #     weights[(begin+3):2:(end-1)] .= 2 * 1 / 3 * Δx
    #     weights[end] = 1 / 3 * Δx
    # end
    # else
    #     # If the number of intervals is odd, average the results of applying Simpson's
    #     # composite 1/3 rule to the first N-1 and last N-1 intervals, applying the
    #     # Trapezoidal rule to the last / first interval.
    #     weights[begin] = 5 / 12 * Δx  # (1/3 + 1/2) / 2 = 5/12
    #     weights[begin + 1] = 13 / 12 * Δx  # (4/3 + 1/3 + 1/2) / 2 = 13/12
    #     weights[(begin + 2):(end - 2)] .= 1 * Δx  # (4/3 + 2/3) / 2 = 1
    #     weights[end - 1] = 13 / 12 * Δx  # ''
    #     weights[end] = 5 / 12 * Δx  # ''
    # end

    return weights
end

function integration_weights_nonuniform!(
    weights::AbstractVector, x::AbstractVector, ::Simpson
)
    N = length(weights) - 1  # Number of intervals
    fill!(weights, 0)

    # Skip the last interval if the number of intervals is odd
    istop = isodd(N) ? lastindex(weights) - 3 : lastindex(weights) - 2

    #! format: off
    for i in firstindex(weights):2:istop
        Δx_0 = x[i + 1] - x[i]
        Δx_1 = x[i + 2] - x[i + 1]
        prefac = (Δx_0 + Δx_1) / 6
        weights[i    ] += prefac * (2 - Δx_1 / Δx_0)
        weights[i + 1] += prefac * (Δx_0 + Δx_1)^2 / (Δx_0 * Δx_1)
        weights[i + 2] += prefac * (2 - Δx_0 / Δx_1)
    end

    if isodd(N)  # This handles the last interval when the number of intervals is odd
        Δx_n   = x[end    ] - x[end - 1]
        Δx_nm1 = x[end - 1] - x[end - 2]
        weights[end    ] += (2 * Δx_n^2 + 3 * Δx_n * Δx_nm1) / (6 * (Δx_nm1 + Δx_n))
        weights[end - 1] += (Δx_n^2 + 3 * Δx_n * Δx_nm1) / (6 * Δx_nm1)
        weights[end - 2] -= Δx_n^3 / (6 * Δx_nm1 * (Δx_nm1 + Δx_n))
    end
    #! format: on
    return weights
end
