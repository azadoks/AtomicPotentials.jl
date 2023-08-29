function AtomicPotential(file::PseudoPotentialIO.Psp8File)
    r = file.rgrid
    ## Local potential
    Vloc = NumericalLocalPotential{RealSpace}(r, file.v_local, file.header.zion)
    ## Non-local potential
    β = map(zip(0:(file.header.lmax), file.projectors)) do (l, proj_l)
        map(proj_l) do proj
            return NumericalQuantity{RealSpace}(r, proj, l, 1)
        end
    end
    D = map(diagm, file.ekb)
    Vnl = NonLocalPotential(β, D, nothing)
    ## Non-linear core correction core charge density
    ρcore = file.rhoc
    ρcore = isnothing(ρcore) ? nothing : NumericalQuantity{RealSpace}(r, ρcore ./ 4π, 0)
    ## Valence charge density for superposition / projection
    ρval = file.rhov
    ρval = isnothing(ρval) ? nothing : NumericalQuantity{RealSpace}(r, ρval ./ 4π, 0)
    ## Atomic potential
    return AtomicPotential(Vloc, Vnl, ρcore, ρval)
end

function _upf_construct_augmentation_q_with_l(upf_file::PseudoPotentialIO.UpfFile)
    T = NumericalQuantity{RealSpace,Numerical,Float64,Vector{Float64}}
    Q = [
        zeros(T, upf_file.header.number_of_proj, upf_file.header.number_of_proj) for
        _ in (0:(2upf_file.header.l_max))
    ]
    for l in 0:(2upf_file.header.l_max)
        # Replace the zeroed quantities with data from the UPF where present
        Q_upf_l = filter(
            qijl -> qijl.angular_momentum == l, upf_file.nonlocal.augmentation.qijls
        )
        for Q_upf in Q_upf_l
            i = Q_upf.first_index
            j = Q_upf.second_index
            f = Q_upf.qijl
            r = upf_file.mesh.r[eachindex(f)]
            Q[l + 1][i, j] = NumericalQuantity{RealSpace}(r, f, l, 2)
            Q[l + 1][j, i] = NumericalQuantity{RealSpace}(r, f, l, 2)
        end
    end
    #TODO would be nice to have this be sparse _and_ symmetric (maybe SparseMatricesCSR)
    return sparse.(Q)
end

function _upf_construct_augmentation_qfcoef(upf_file::PseudoPotentialIO.UpfFile)
    #TODO check correctness
    r = upf_file.mesh.r
    r2 = upf_file.mesh.r .^ 2
    nqf = upf_file.nonlocal.augmentation.nqf
    nqlc = 2upf_file.header.l_max + 1

    T = NumericalQuantity{RealSpace,Numerical,Float64,Vector{Float64}}
    Q = [
        zeros(T, upf_file.header.number_of_proj, upf_file.header.number_of_proj) for
        _ in (0:(2upf_file.header.l_max))
    ]

    for (Q_upf, Qfcoef_upf) in
        zip(upf_file.nonlocal.augmentation.qijs, upf_file.nonlocal.augmentation.qfcoefs)
        qfcoef = reshape(Qfcoef_upf.qfcoef, nqf, nqlc)
        rinner = upf_file.nonlocal.augmentation.rinner

        i = Q_upf.first_index
        j = Q_upf.second_index

        li = upf_file.nonlocal.betas[i].angular_momentum
        lj = upf_file.nonlocal.betas[j].angular_momentum

        for l in abs(li - lj):2:(li + lj)  # parity-constrained (c_CG = 0)
            qij = copy(Q_upf.qij)
            ircut = findfirst(i -> r[i] > rinner[l + 1], eachindex(r)) - 1
            poly = Polynomial(qfcoef[:, l + 1])
            qij[1:ircut] = r[1:ircut] .^ (l + 2) .* poly.(r2[1:ircut])
            rij = r[eachindex(qij)]
            Q[l + 1][i, j] = NumericalQuantity{RealSpace}(rij, qij, l, 2)
            Q[l + 1][j, i] = NumericalQuantity{RealSpace}(rij, qij, l, 2)
        end
    end
    return sparse.(Q)
end

function AtomicPotential(upf_file::PseudoPotentialIO.UpfFile)
    r = upf_file.mesh.r
    ## Local potential
    Vloc = NumericalLocalPotential{RealSpace}(
        r,
        upf_file.local_ ./ 2,  # Ry -> Ha
        upf_file.header.z_valence,
    )
    ## Non-local potential
    β = map(0:(upf_file.header.l_max)) do l
        betas_l = filter(beta -> beta.angular_momentum == l, upf_file.nonlocal.betas)
        map(betas_l) do beta
            i_r_cut = beta.cutoff_radius_index
            f = isnothing(i_r_cut) ? beta.beta : beta.beta[begin:i_r_cut]  # 1/√a₀ * a₀²
            r_beta = isnothing(i_r_cut) ? r : r[begin:i_r_cut]
            return NumericalQuantity{RealSpace}(r_beta, f, l, 1)
        end
    end
    D = map(0:(upf_file.header.l_max)) do l
        betas_l = filter(beta -> beta.angular_momentum == l, upf_file.nonlocal.betas)
        i_betas_l = map(beta -> beta.index, betas_l)
        return upf_file.nonlocal.dij[i_betas_l, i_betas_l] ./ 2  # Ry -> Ha
    end
    if isnothing(upf_file.nonlocal.augmentation)
        augmentation = nothing
    else
        q = map(0:(upf_file.header.l_max)) do l
            betas_l = filter(beta -> beta.angular_momentum == l, upf_file.nonlocal.betas)
            i_betas_l = map(beta -> beta.index, betas_l)
            return upf_file.nonlocal.augmentation.q[i_betas_l, i_betas_l]
        end
        if upf_file.nonlocal.augmentation.q_with_l
            Q = _upf_construct_augmentation_q_with_l(upf_file)
        elseif upf_file.nonlocal.augmentation.nqf > 0
            Q = _upf_construct_augmentation_qfcoef(upf_file)
        else
            error("q_with_l == false and nqf == 0, unsure what to do...")
        end
        augmentation = Augmentation(Q, q)
    end
    Vnl = NonLocalPotential(β, D, augmentation)
    ## Core charge density
    if isnothing(upf_file.nlcc) || all(iszero.(upf_file.rhoatom))
        ρcore = nothing
    else
        ρcore_f = upf_file.nlcc
        ρcore = NumericalQuantity{RealSpace}(r[eachindex(ρcore_f)], ρcore_f, 0)
    end
    ## Valence charge density
    if isnothing(upf_file.rhoatom) || all(iszero.(upf_file.rhoatom))
        ρval = nothing
    else
        ρval_f = upf_file.rhoatom ./ 4π
        ρval = NumericalQuantity{RealSpace}(r[eachindex(ρval_f)], ρval_f, 0, 2)
    end
    ## Pseudo-atomic states
    if !isnothing(upf_file.pswfc)
        χ = map(0:(upf_file.header.l_max)) do l
            chis_l = filter(chi -> chi.l == l, upf_file.pswfc)
            return map(chis_l) do chi
                chi_f = chi.chi
                chi_r = r[eachindex(chi_f)]
                return NumericalQuantity{RealSpace}(chi_r, chi_f, l, 1)
            end
        end
        filter!(!isempty, χ)
    else
        χ = Nothing[]
    end
    ## Atomic potential
    return AtomicPotential(Vloc, Vnl, ρcore, ρval, χ)
end

function AtomicPotential(hgh_file::PseudoPotentialIO.HghFile)
    ## Local potential
    Zval = Float64(sum(hgh_file.zion))
    cloc = hgh_file.cloc
    length(cloc) > 4 && error("length(cloc) > 4 not supported.")
    if length(cloc) < 4
        n_extra = 4 - length(cloc)
        cloc = [cloc; zeros(Float64, n_extra)]
    end
    Vloc = HghLocalPotential{RealSpace}(hgh_file.rloc, cloc, Zval)
    ## Non-local potential
    # In HGH pseudos, the number of projectors is defined by the maximum angular momentum
    # and the number of coupling constants at each angular momentum.
    β = map(0:(hgh_file.lmax)) do l
        map(1:size(hgh_file.h[l + 1], 1)) do i
            return HghProjector{RealSpace}(hgh_file.rp[l + 1], i, l)  # Ha / √a₀
        end
    end
    D = hgh_file.h
    Vnl = NonLocalPotential(β, D)
    ## Atomic potential
    return AtomicPotential(Vloc, Vnl)
end
