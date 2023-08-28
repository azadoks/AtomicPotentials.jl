@doc raw"""
Ultrasoft augmentation.

!!! note
    The implementation of parsing, rft, and irft for ultrasoft augmentation is untested
    and likely incorrect.

```math
n(\mathbf{r}) = \sum_i \left[
     |\phi_i(\mathbf{r})|^2
     \sum_{nm,I} Q_{nm}^I(\mathbf{r})
     \bra{\phi_i}\ket{\beta_n^I} \bra{\beta_m^I} \ket{\phi_i}
\right]
```

```math
S = 1 + \sum_{nm,I} q_{nm}^I \ket{\beta_n^I} \bra{\beta_m^I}
```

```math
D^_{nm}^{I} = D_{nm}^{(0)} + \int d\mathbf{r} V_{eff}(\mathbf{r}) Q_{nm}^I(\mathbf{r})
```
where
```math
V_{eff}(\mathbf{r}) = V_{loc}(\mathbf{r}) + V_{xc}(\mathbf{r}) + V_H(\mathbf{r})
```

(CASTEP Documentation)[https://www.tcm.phy.cam.ac.uk/castep/documentation/WebHelp/content/modules/castep/thcastepultrapseudo.htm]
"""
struct Augmentation{
    S<:EvaluationSpace,
    VF<:AbstractVector{<:AbstractMatrix{<:AbstractQuantity{S}}},
    VC<:AbstractVector{<:AbstractMatrix{<:Real}},
}
    "Augmentation functions."
    functions::VF
    "Augmentation charges / coupling constants"
    coupling::VC
end

Base.broadcastable(qty::Augmentation) = Ref(qty)
Base.isempty(qty::Augmentation) = false

function Base.convert(::Type{T}, x::Augmentation) where {T<:Real}
    return Augmentation(
        map(Q_l -> map(Q -> convert(T, Q), Q_l), x.functions),
        map(q_l -> map(q -> convert(T, q), q_l), x.coupling),
    )
end

function Adapt.adapt(to, x::Augmentation)
    return Augmentation(adapt(to, x.functions), adapt(to, x.coupling))
end

function rft(augmentation::Augmentation{RealSpace}, q::AbstractVector; kwargs...)
    functions = map(f_l -> map(f -> rft(f, q; kwargs...), f_l), augmentation.functions)
    return Augmentation(functions, augmentation.coupling)
end

function irft(augmentation::Augmentation{FourierSpace}, r::AbstractVector; kwargs...)
    functions = map(F_l -> map(F -> irft(F, r; kwargs...), F_l), augmentation.functions)
    return Augmentation(functions, augmentation.coupling)
end
