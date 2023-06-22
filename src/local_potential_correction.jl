abstract type AbstractLocalPotentialCorrection{S,Analytical} <: AbstractAtomicQuantity{S,Analytical} end

struct CoulombLocalPotentialCorrection{S,Analytical} <: AbstractLocalPotentialCorrection{S,Analytical}
    Z
end
function CoulombLocalPotentialCorrection(quantity::LocalPotential{S}) where {S}
    return CoulombLocalPotentialCorrection{S,Analytical}(quantity.Z)
end
(correction::CoulombLocalPotentialCorrection{RealSpace})(_) = -correction.Z  #* == (-Z/r) * r
(correction::CoulombLocalPotentialCorrection{FourierSpace})(q) = -correction.Z / q^2

struct ErfLocalPotentialCorrection{S,Analytical} <: AbstractLocalPotentialCorrection{S,Analytical}
    Z
end
(correction::ErfLocalPotentialCorrection{RealSpace})(r) = -correction.Z * erf(r)  #* == (-Z/r * erf(r)) * r
(correction::ErfLocalPotentialCorrection{FourierSpace})(q) = -correction.Z * exp(-q^2 / 4) / q^2
