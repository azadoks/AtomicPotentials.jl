function AtomicPotential(file::PseudoPotentialIO.Psp8File)
    r = file.rgrid
    ## Local potential
    Vloc = NumericalLocalPotential{RealSpace}(r, file.v_local, file.header.zion)
    ## Non-local potential
    β = map(enumerate(file.projectors)) do (i_l, proj_l)
        map(enumerate(proj_l)) do (n, proj)
            NumericalQuantity{RealSpace}(r, proj, i_l - 1, 1)
        end
    end
    β = collect(Iterators.flatten(β))
    # Here, we build up a symmetric block-diagonal coupling matrix in which each column
    # corresponds to one of the projectors above
    D = zeros(Float64, length(β), length(β))
    i = 1
    for ekb_l in file.ekb
        ekb_l = diagm(ekb_l)
        n_projectors_l = size(ekb_l, 1)
        j = i + n_projectors_l - 1
        D[i:j, i:j] .= ekb_l
        i += n_projectors_l
    end
    Vnl = NonLocalPotential(β, D, nothing)
    ## Non-linear core correction core charge density
    ρcore = file.rhoc
    ρcore = isnothing(ρcore) ? nothing : NumericalQuantity{RealSpace}(r, ρcore ./ 4π, 0)
    ## Valence charge density for superposition / projection
    ρval = file.rhov
    ρval = isnothing(ρval) ? nothing : NumericalQuantity{RealSpace}(r, ρval ./ 4π, 0)
    ## Atomic potential
    return AtomicPotential(Vloc, Vnl, ρval, ρcore)
end

# function _upf_construct_augmentation_q_with_l(
#     upf_file::PseudoPotentialIO.UpfFile, interpolation_method
# )
#     Q = OffsetVector(
#         [
#             Matrix{AugmentationFunction{RealSpace,Numerical}}(
#                 undef, upf_file.header.number_of_proj, upf_file.header.number_of_proj
#             ) for l in 0:(2upf_file.header.l_max)
#         ],
#         0:(2upf_file.header.l_max),
#     )
#     for l in 0:(2upf_file.header.l_max)
#         # Fill Q with zeroed quantities
#         for i in 1:(upf_file.header.number_of_proj), j in 1:(upf_file.header.number_of_proj)
#             r = upf_file.mesh.r
#             f = zero(r)
#             interpolator = Interpolation.construct_interpolator(r, f, interpolation_method)
#             Q[l][i, j] = AugmentationFunction{RealSpace,Numerical}(
#                 r, f, interpolator, 0, i, j, l
#             )
#         end
#         # Replace the zeroed quantities with data from the UPF where present
#         Q_upf_l = filter(
#             qijl -> qijl.angular_momentum == l, upf_file.nonlocal.augmentation.qijls
#         )
#         for Q_upf in Q_upf_l
#             n1 = Q_upf.first_index
#             n2 = Q_upf.second_index
#             f = Q_upf.qijl
#             r = upf_file.mesh.r[eachindex(f)]
#             interpolator = Interpolation.construct_interpolator(r, f, interpolation_method)
#             Q[l][n1, n2] = AugmentationFunction{RealSpace,Numerical}(
#                 r, f, 0, interpolator, n1, n2, l
#             )
#             Q[l][n2, n1] = AugmentationFunction{RealSpace,Numerical}(
#                 r, f, 0, interpolator, n2, n1, l
#             )
#         end
#     end
#     return Q
# end

# @views function _upf_construct_augmentation_qfcoef(
#     upf_file::PseudoPotentialIO.UpfFile, interpolation_method
# )
#     #TODO check correctness
#     r = upf_file.mesh.r
#     r2 = upf_file.mesh.r .^ 2
#     nqf = upf_file.nonlocal.augmentation.nqf
#     nqlc = 2upf_file.header.l_max + 1

#     Q = OffsetVector(
#         [
#             Matrix{AugmentationFunction{RealSpace,Numerical}}(
#                 undef, upf_file.header.number_of_proj, upf_file.header.number_of_proj
#             ) for l in 0:(2upf_file.header.l_max)
#         ],
#         0:(2upf_file.header.l_max),
#     )
#     for l in 0:(2upf_file.header.l_max),
#         i in 1:(upf_file.header.number_of_proj),
#         j in 1:(upf_file.header.number_of_proj)
#         # Fill Q with zero vectors
#         f = zero(r)
#         interpolator = Interpolation.construct_interpolator(r, f, interpolation_method)
#         Q[l][i, j] = AugmentationFunction{RealSpace,Numerical}(r, f, interpolator, i, j, l)
#     end
#     for (Q_upf, Qfcoef_upf) in
#         zip(upf_file.nonlocal.augmentation.qijs, upf_file.nonlocal.augmentation.qfcoefs)
#         # Replace the zero vectors with datat from the UPF where present
#         # It's not worth the effort to make these into OffsetVectors zero-indexed for l.
#         qfcoef = reshape(Qfcoef_upf.qfcoef, nqf, nqlc)
#         rinner = upf_file.nonlocal.augmentation.rinner

#         i = Q_upf.first_index
#         j = Q_upf.second_index

#         li = upf_file.nonlocal.betas[i].angular_momentum
#         lj = upf_file.nonlocal.betas[j].angular_momentum

#         for l in abs(li - lj):2:(li + lj)
#             qij = copy(Q_upf.qij)
#             ircut = findfirst(i -> r[i] > rinner[l + 1], eachindex(r)) - 1
#             poly = Polynomial(qfcoef[:, l + 1])
#             qij[1:ircut] = r[1:ircut] .^ (l + 2) .* poly.(r2[1:ircut])

#             n1 = Q_upf.first_index
#             n2 = Q_upf.second_index
#             rij = r[eachindex(qij)]
#             interpolator = Interpolation.construct_interpolator(
#                 rij, qij, interpolation_method
#             )

#             Q[l][n1, n2] = AugmentationFunction{RealSpace,Numerical}(
#                 rij, qij, interpolator, 0, n1, n2, l
#             )
#             Q[l][n2, n1] = AugmentationFunction{RealSpace,Numerical}(
#                 rij, qij, interpolator, 0, n2, n1, l
#             )
#         end
#     end
#     return Q
# end

function AtomicPotential(upf_file::PseudoPotentialIO.UpfFile)
    r = upf_file.mesh.r
    ## Local potential
    Vloc = NumericalLocalPotential{RealSpace}(
        r,
        upf_file.local_ ./ 2,  # Ry -> Ha
        upf_file.header.z_valence,
    )
    ## Non-local potential
    β = map(upf_file.nonlocal.betas) do beta
        f = beta.beta  # 1/√a₀ * a₀²
        i_r_cut = beta.cutoff_radius_index
        f = isnothing(i_r_cut) ? f : f[begin:i_r_cut]
        r_beta = isnothing(i_r_cut) ? r : r[begin:i_r_cut]
        return NumericalQuantity{RealSpace}(r_beta, f, beta.angular_momentum, 1)
    end
    D = upf_file.nonlocal.dij ./ 2  # Ry -> Ha
    #TODO: augmentation
    Vnl = NonLocalPotential(β, D)
    ## Core charge density
    if isnothing(upf_file.nlcc)
        ρcore = nothing
    else
        ρcore_f = upf_file.nlcc
        ρcore = NumericalQuantity{RealSpace}(r[eachindex(ρcore_f)], ρcore_f, 0)
    end
    ## Valence charge density
    if isnothing(upf_file.rhoatom)
        ρval = nothing
    else
        ρval_f = upf_file.rhoatom ./ 4π
        ρval = NumericalQuantity{RealSpace}(r[eachindex(ρval_f)], ρval_f, 0, 2)
    end
    ## Pseudo-atomic states
    if !isnothing(upf_file.pswfc)
        χ = map(upf_file.pswfc) do chi
            chi_f = chi.chi
            chi_r = r[eachindex(chi_f)]
            return NumericalQuantity{RealSpace}(chi_r, chi_f, chi.l, 1)
        end
    else
        χ = Nothing[]
    end
    return AtomicPotential(Vloc, Vnl, ρval, ρcore, χ)
end

function AtomicPotential(hgh_file::PseudoPotentialIO.HghFile)
    ## Local potential
    Zval = Float64(sum(hgh_file.zion))
    cloc = hgh_file.cloc
    length(cloc) > 4 && error("length(cloc) > 4 not supported.")
    if length(cloc) < 4
        n_extra = 4 - length(cloc)
        cloc = Tuple([cloc; zeros(Float64, n_extra)])
    end
    Vloc = HghLocalPotential{RealSpace}(hgh_file.rloc, cloc, Zval)

    ## Non-local potential
    lmax = hgh_file.lmax
    # In HGH pseudos, the number of projectors is defined by the maximum angular momentum
    # and the number of coupling constants at each angular momentum.
    β = map(0:lmax) do l
        map(1:size(hgh_file.h[l + 1], 1)) do n
            return HghProjector{RealSpace}(hgh_file.rp[l + 1], n, l)
        end
    end
    β = collect(Iterators.flatten(β))  # Ha / √a₀
    # Here, we build up a symmetric block-diagonal coupling matrix in which each column
    # corresponds to one of the projectors above
    D = Matrix{Float64}(undef, length(β), length(β))  # 1 / Ha
    i = 1
    for h_l in hgh_file.h
        n_projectors_l = size(h_l, 1)
        j = i + n_projectors_l - 1
        D[i:j, i:j] .= h_l
        i += n_projectors_l
    end
    Vnl = NonLocalPotential(β, D)

    ## Atomic potential
    return AtomicPotential(Vloc, Vnl)
end
