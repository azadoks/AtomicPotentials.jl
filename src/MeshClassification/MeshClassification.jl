module MeshClassification

export Mesh
export Uniform
export Log
export ShiftedLog

@doc raw"""
Abstract 1-dimensional mesh.
"""
abstract type Mesh{T<:Real} end

Base.broadcastable(mesh::Mesh) = Ref(mesh)

function Mesh(x::AbstractVector; kwargs...)
    for type in (Uniform, Log, ShiftedLog)
        mesh = type(x)
        x_calculated = mesh.(eachindex(x))
        all(isapprox.(x_calculated, x; kwargs...)) && return mesh
    end
    return error("Could not fit a known mesh type to the data.")
end

@doc raw"""
Evaluate the mesh at index `i`.
"""
function (::Mesh)(::Integer) end

@doc raw"""
Return a callable which evaluates the gradient of the mesh.
"""
function gradient(::Mesh) end

@doc raw"""
Uniform mesh of the form
```math
x\_i = x\_0 + dx (i - i).
```
"""
struct Uniform{T<:Real} <: Mesh{T}
    "Spacing."
    dx::T
    "First mesh point."
    x₀::T
end

function Uniform(x::AbstractVector{T}) where {T<:Real}
    dx = x[2] - x[1]
    x₀ = x[1]
    return Uniform{T}(dx, x₀)
end

function (mesh::Uniform{T})(i::Integer)::T where {T<:Real}
    return mesh.x₀ + mesh.dx * (i - 1)
end

function gradient(mesh::Uniform{T}) where {T<:Real}
    return f(::Integer)::T where {T} = mesh.x₀
end

@doc raw"""
Logarithmic mesh of the form
```math
x\_i = x_c exp(d\mathrm{ln} (i - 1)).
```
"""
struct Log{T<:Real} <: Mesh{T}
    "Spacing."
    dln::T
    "Scaling."
    xc::T
end

function Log(x::AbstractVector{T}) where {T<:Real}
    xc = first(x)
    dln = log(last(x) / first(x)) / (length(x) - 1)
    return Log{T}(dln, xc)
end

function (mesh::Log{T})(i::Integer)::T where {T<:Real}
    return mesh.xc * exp(mesh.dln * (i - 1))
end

function gradient(mesh::Log{T}) where {T<:Real}
    return f(i::Integer)::T where {T} = mesh.dln * mesh.xc * exp(mesh.dln * (i - 1))
end

@doc raw"""
Shifted logarithmic mesh of the form
```math
x\_i = x_c (exp(d\mathrm{ln} (i - 1)) - 1).
```
"""
struct ShiftedLog{T<:Real} <: Mesh{T}
    "Spacing."
    dln::T
    "Scaling."
    xc::T
end

function ShiftedLog(x::AbstractVector{T}) where {T<:Real}
    dln = -x[2]^2 / (2x[2] - x[3])
    xc = log(x[3] / x[2] - 1)
    return ShiftedLog{T}(xc, dln)
end

function (mesh::ShiftedLog{T})(i::Integer)::T where {T<:Real}
    return mesh.xc * (exp(mesh.dln * (i - 1)) - 1)
end

function gradient(mesh::ShiftedLog{T}) where {T<:Real}
    return f(i::Integer)::T where {T} = mesh.dln * mesh.xc * exp(mesh.dln * (i - 1))
end
end
