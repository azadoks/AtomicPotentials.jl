@doc raw"""
Trapezoidal rule quadrature implemented as described in the
[Wikipedia article](https://en.wikipedia.org/wiki/Trapezoidal_rule).

Non-uniform grid
```math
\int_a^b f(x) dx \approx
\sum_{i=0}^{N} \frac{f(x_{i-1}) + f(x_i)}{2} (x_i - x_{i-1})
```

Uniform grid
```math
\int_a^b f(x) dx \approx
\Delta x \left( \sum_{i=1}^{N-1} f(x_i) + \frac{f(x_N) + f(x_0)}{2} \right)
```
"""
struct Trapezoid{U<:Uniformity} <: QuadratureMethod end

function integration_weights!(
    weights::AbstractVector, x::AbstractVector, ::Trapezoid{Uniform}
)
    length(weights) == length(x) ||
        throw(DimensionMismatch("weights and x arrays must have a common length"))
    length(x) >= 2 || throw(ArgumentError("x must have a length of at least 2"))

    Δx = x[begin + 1] - x[begin]

    weights[begin] = Δx / 2
    weights[(begin + 1):(end - 1)] .= Δx
    weights[end] = Δx / 2

    return weights
end

function integration_weights!(
    weights::AbstractVector, x::AbstractVector, ::Trapezoid{NonUniform}
)
    length(weights) == length(x) ||
        throw(DimensionMismatch("weights and x arrays must have a common length"))
    length(x) >= 2 || throw(ArgumentError("x must have a length of at least 2"))

    weights[begin] = (x[begin + 1] - x[begin]) / 2
    for i in (firstindex(weights) + 1):(lastindex(weights) - 1)
        # Δx[i] + Δx[i-1] = (x[i+1] - x[i]) + (x[i] - x[i-1]) = x[i+i] - x[i-1]
        weights[i] = (x[i + 1] - x[i - 1]) / 2
    end
    weights[end] = (x[end] - x[end - 1]) / 2

    return weights
end
