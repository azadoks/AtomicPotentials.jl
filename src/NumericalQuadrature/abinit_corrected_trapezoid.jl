@doc raw"""
ABINIT quadrature -- a collection of composite closed Newton-Cotes rules. Behavior depends
on the number of available points. This method _only_ supports quantities on uniform
meshes!

```math
\int_a^b f(x) dx \approx \begin{cases}
    \frac{\Delta x}{72} \left[
        23.75 f_1 + 95.10 f_2 + 55.20 f_3 + 79.30 f_4 + 70.65 f_5 +
        72 \sum_{i=6}^{N-5} f_i +
        70.65 f_{N-4} 79.30 f_{N-3} + 55.20 f_{N-2} + 95.10 f_{N-1} + 23.75 f_{N}
    \right] & N >= 10 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 49 f_4 +
        48 f_5 +
        49 f_6 + 43 f_7 + 59 f_8 + 17 f_9
    \right] & N = 9 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 49 f_4 +
        49 f_5 + 43 f_6 + 59 f_7 + 17 f_8
    \right] & N = 8 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 50 f_4 + 43 f_5 + 59 f_6 + 17 f_7
    \right] & N = 7 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 44 f_3 + 44 f_4 + 59 f_5 + 17 f_6
    \right] & N = 6 \\
    \frac{\Delta x}{3} \left[
        f_1 + 4 f_2 + 2 f_3 + 4 f_4 + f_5
    \right] & N = 5 \\
    \frac{\Delta x}{8} \left[
        3 f_1 + 9 f_2 + 9 f_3 + 3 f_4
    \right] & N = 4 \\
    \frac{\Delta x}{3} \left[
        1 f_1 + 8 f_2 + 1 f_3
    \right] & N = 3 \\
    \frac{\Delta x}{2} \left[
        f_1 + f_2
    \right] & N = 2 \\
    \frac{\Delta x} \left[
        f_1
    \right] & N = 1
\end{cases}
```
"""
struct AbinitCorrectedTrapezoid <: QuadratureMethod end

function integration_weights!(
    weights::AbstractVector, x::AbstractVector, ::AbinitCorrectedTrapezoid
)
    is_uniform(x) || error("AbinitCorrectedTrapezoid does not support nonuniform meshes")
    Δx = x[begin + 1] - x[begin]
    #! format: off
    if length(weights) >= 10
        weights[begin    ] = 23.75 / 72 * Δx
        weights[begin + 1] = 95.10 / 72 * Δx
        weights[begin + 2] = 55.20 / 72 * Δx
        weights[begin + 3] = 79.30 / 72 * Δx
        weights[begin + 4] = 70.65 / 72 * Δx
        weights[(begin + 5):(end - 5)] .= Δx
        weights[end   - 4] = 70.65 / 72 * Δx
        weights[end   - 3] = 79.30 / 72 * Δx
        weights[end   - 2] = 55.20 / 72 * Δx
        weights[end   - 1] = 95.10 / 72 * Δx
        weights[end      ] = 23.75 / 72 * Δx
    elseif length(weights) == 9
        weights[begin    ] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[begin + 4] = Δx
        weights[end   - 3] = 49 / 48 * Δx
        weights[end   - 2] = 43 / 48 * Δx
        weights[end   - 1] = 59 / 48 * Δx
        weights[end      ] = 17 / 48 * Δx
    elseif length(weights) == 8
        weights[begin    ] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[end   - 3] = 49 / 48 * Δx
        weights[end   - 2] = 43 / 48 * Δx
        weights[end   - 1] = 59 / 48 * Δx
        weights[end      ] = 17 / 48 * Δx
    elseif length(weights) == 7
        weights[begin    ] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 50 / 48 * Δx
        weights[end   - 2] = 43 / 48 * Δx
        weights[end   - 1] = 59 / 48 * Δx
        weights[end      ] = 17 / 48 * Δx
    elseif length(weights) == 6
        weights[begin    ] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 44 / 48 * Δx
        weights[end   - 2] = 44 / 48 * Δx
        weights[end   - 1] = 59 / 48 * Δx
        weights[end      ] = 17 / 48 * Δx
    elseif length(weights) == 5
        weights[begin    ] = 1 / 3 * Δx
        weights[begin + 1] = 4 / 3 * Δx
        weights[begin + 2] = 2 / 3 * Δx
        weights[end   - 1] = 4 / 3 * Δx
        weights[end      ] = 1 / 3 * Δx
    elseif length(weights) == 4
        weights[begin    ] = 3 / 8 * Δx
        weights[begin + 1] = 9 / 8 * Δx
        weights[end   - 1] = 9 / 8 * Δx
        weights[end      ] = 3 / 8 * Δx
    elseif length(weights) == 3
        weights[begin    ] = 1 / 3 * Δx
        weights[begin + 1] = 8 / 3 * Δx
        weights[end      ] = 1 / 3 * Δx
    elseif length(weights) == 2
        weights[begin    ] = 1 / 2 * Δx
        weights[end      ] = 1 / 2 * Δx
    elseif length(weights) == 1
        weights[begin    ] = Δx
    end
    #! format: on
    return weights
end
