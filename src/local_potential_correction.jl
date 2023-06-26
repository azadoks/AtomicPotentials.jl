# TODO: figure out how to work this in without causing too much trouble
# TODO: design criteria
# TODO:  (1) selectable at time of fht (i.e. can't be baked-into the LocalPotential or AtomicPotential)
# TODO:  (2) respects fht interface
# TODO: pretty sure can't have both, leaning towards breaking (1) to keep (2)
abstract type AbstractLocalPotentialCorrection{S,Analytical} <:
              AbstractAtomicQuantity{S,Analytical} end

struct CoulombLocalPotentialCorrection{S,Analytical} <:
       AbstractLocalPotentialCorrection{S,Analytical}
    Z
end
function CoulombLocalPotentialCorrection(quantity::LocalPotential{S}) where {S}
    return CoulombLocalPotentialCorrection{S,Analytical}(quantity.Z)
end
(correction::CoulombLocalPotentialCorrection{RealSpace})(_) = -correction.Z  #* == (-Z/r) * r
(correction::CoulombLocalPotentialCorrection{FourierSpace})(q) = -correction.Z / q^2

struct ErfLocalPotentialCorrection{S,Analytical} <:
       AbstractLocalPotentialCorrection{S,Analytical}
    Z
end
(correction::ErfLocalPotentialCorrection{RealSpace})(r) = -correction.Z * erf(r)  #* == (-Z/r * erf(r)) * r
function (correction::ErfLocalPotentialCorrection{FourierSpace})(q)
    return -correction.Z * exp(-q^2 / 4) / q^2
end
