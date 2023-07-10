function AtomicPotential(
    psp8_file::PseudoPotentialIO.Psp8File, interpolation_method=Interpolation.Spline(4)
)
    identifier = psp8_file.identifier
    symbol = Symbol(PseudoPotentialIO.element(psp8_file).symbol)
    lmax = psp8_file.header.lmax
    r = psp8_file.rgrid
    j = 0.0

    Vloc = LocalPotential{RealSpace,Numerical}(
        r,
        psp8_file.v_local,
        Interpolation.construct_interpolator(r, psp8_file.v_local, interpolation_method),
        psp8_file.header.zion,
    )

    β = map(enumerate(psp8_file.projectors)) do (i, proj_l)
        l = i - 1
        map(enumerate(proj_l)) do (n, proj)
            f = r .* proj
            interpolator = Interpolation.construct_interpolator(r, f, interpolation_method)
            KleinmanBylanderProjector{RealSpace,Numerical}(r, f, interpolator, n, l, j)
        end
    end
    β = OffsetVector(β, 0:lmax)
    D = OffsetVector(map(l -> diagm(psp8_file.ekb[l + 1]), 0:lmax), 0:lmax)
    Vnl = NonLocalPotential{RealSpace,Numerical,eltype(eltype(β))}(β, D)

    if isnothing(psp8_file.rhoc)
        ρcore = nothing
    else
        ρcore_f = r .^ 2 .* psp8_file.rhoc ./ 4π
        ρcore_interpolator = Interpolation.construct_interpolator(
            r, ρcore_f, interpolation_method
        )
        ρcore = ChargeDensity{RealSpace,Numerical}(r, ρcore_f, ρcore_interpolator)
    end

    if isnothing(psp8_file.rhov)
        ρval = nothing
    else
        ρval_f = r .^ 2 .* psp8_file.rhov ./ 4π
        ρval_interpolator = Interpolation.construct_interpolator(
            r, ρval_f, interpolation_method
        )
        ρval = ChargeDensity{RealSpace,Numerical}(r, ρval_f, ρval_interpolator)
    end

    states = OffsetVector(Vector{Nothing}[], 0:-1)
    aug = nothing

    return AtomicPotential(identifier, symbol, Vloc, Vnl, ρval, ρcore, states, aug)
end

function _upf_construct_augmentation_q_with_l(
    upf_file::PseudoPotentialIO.UpfFile, interpolation_method
)
    Q = OffsetVector(
        [
            Matrix{AugmentationFunction{RealSpace,Numerical}}(
                undef, upf_file.header.number_of_proj, upf_file.header.number_of_proj
            ) for l in 0:(2upf_file.header.l_max)
        ],
        0:(2upf_file.header.l_max),
    )
    for l in 0:(2upf_file.header.l_max)
        # Fill Q with zeroed quantities
        for i in 1:(upf_file.header.number_of_proj), j in 1:(upf_file.header.number_of_proj)
            r = upf_file.mesh.r
            f = zero(r)
            interpolator = Interpolation.construct_interpolator(r, f, interpolation_method)
            Q[l][i, j] = AugmentationFunction{RealSpace,Numerical}(
                r, f, interpolator, i, j, l
            )
        end
        # Replace the zeroed quantities with data from the UPF where present
        Q_upf_l = filter(
            qijl -> qijl.angular_momentum == l, upf_file.nonlocal.augmentation.qijls
        )
        for Q_upf in Q_upf_l
            n1 = Q_upf.first_index
            n2 = Q_upf.second_index
            f = Q_upf.qijl
            r = upf_file.mesh.r[eachindex(f)]
            interpolator = Interpolation.construct_interpolator(r, f, interpolation_method)
            Q[l][n1, n2] = AugmentationFunction{RealSpace,Numerical}(
                r, f, interpolator, n1, n2, l
            )
            Q[l][n2, n1] = AugmentationFunction{RealSpace,Numerical}(
                r, f, interpolator, n2, n1, l
            )
        end
    end
    return Q
end

@views function _upf_construct_augmentation_qfcoef(
    upf_file::PseudoPotentialIO.UpfFile, interpolation_method
)
    #TODO check correctness
    r = upf_file.mesh.r
    r2 = upf_file.mesh.r .^ 2
    nqf = upf_file.nonlocal.augmentation.nqf
    nqlc = 2upf_file.header.l_max + 1

    Q = OffsetVector(
        [
            Matrix{AugmentationFunction{RealSpace,Numerical}}(
                undef, upf_file.header.number_of_proj, upf_file.header.number_of_proj
            ) for l in 0:(2upf_file.header.l_max)
        ],
        0:(2upf_file.header.l_max),
    )
    for l in 0:(2upf_file.header.l_max),
        i in 1:(upf_file.header.number_of_proj),
        j in 1:(upf_file.header.number_of_proj)
        # Fill Q with zero vectors
        f = zero(r)
        interpolator = Interpolation.construct_interpolator(r, f, interpolation_method)
        Q[l][i, j] = AugmentationFunction{RealSpace,Numerical}(r, f, interpolator, i, j, l)
    end
    for (Q_upf, Qfcoef_upf) in
        zip(upf_file.nonlocal.augmentation.qijs, upf_file.nonlocal.augmentation.qfcoefs)
        # Replace the zero vectors with datat from the UPF where present
        # It's not worth the effort to make these into OffsetVectors zero-indexed for l.
        qfcoef = reshape(Qfcoef_upf.qfcoef, nqf, nqlc)
        rinner = upf_file.nonlocal.augmentation.rinner

        i = Q_upf.first_index
        j = Q_upf.second_index

        li = upf_file.nonlocal.betas[i].angular_momentum
        lj = upf_file.nonlocal.betas[j].angular_momentum

        for l in abs(li - lj):2:(li + lj)
            qij = copy(Q_upf.qij)
            ircut = findfirst(i -> r[i] > rinner[l + 1], eachindex(r)) - 1
            poly = Polynomial(qfcoef[:, l + 1])
            qij[1:ircut] = r[1:ircut] .^ (l + 2) .* poly.(r2[1:ircut])

            n1 = Q_upf.first_index
            n2 = Q_upf.second_index
            rij = r[eachindex(qij)]
            interpolator = Interpolation.construct_interpolator(
                rij, qij, interpolation_method
            )

            Q[l][n1, n2] = AugmentationFunction{RealSpace,Numerical}(
                rij, qij, interpolator, n1, n2, l
            )
            Q[l][n2, n1] = AugmentationFunction{RealSpace,Numerical}(
                rij, qij, interpolator, n2, n1, l
            )
        end
    end
    return Q
end

function AtomicPotential(
    upf_file::PseudoPotentialIO.UpfFile, interpolation_method=Interpolation.Spline(4)
)
    identifier = upf_file.identifier
    symbol = Symbol(PseudoPotentialIO.element(upf_file).symbol)
    lmax = upf_file.header.l_max
    r = upf_file.mesh.r
    j = 0.0

    ## Local potential
    Vloc = LocalPotential{RealSpace,Numerical}(
        r,
        upf_file.local_ ./ 2,  # Ry -> Ha
        Interpolation.construct_interpolator(r, upf_file.local_ ./ 2, interpolation_method),
        upf_file.header.z_valence,
    )

    ## Non-local potential
    # Indices in upf.nonlocal.betas for projectors at each angular momentum
    iβ_upf = map(0:lmax) do l
        β_upf = upf_file.nonlocal.betas
        return filter(i -> β_upf[i].angular_momentum == l, eachindex(β_upf))
    end
    iβ_upf = OffsetVector(iβ_upf, 0:lmax)
    # Number of projectors at each angular momentum
    nβ = OffsetArray(length.(iβ_upf), 0:lmax)
    # Find the first/last indices in upf.nonlocal.dij for each angular momentum so the
    # sub-arrays D[l][n,m] can be extracted
    cumul_nβ = [0, cumsum(nβ)...]
    # UPFs store the projectors multiplied by the radial grid, so we multiply again by the
    # grid for consistency.
    # If the cutoff radius index is stored, we use it to truncate the projector for
    # numerical stability.
    β = map(0:lmax) do l
        map(iβ_upf[l]) do n
            βln_data = upf_file.nonlocal.betas[n].beta

            ir_cut = upf_file.nonlocal.betas[n].cutoff_radius_index
            ir_cut = isnothing(ir_cut) ? length(βln) : ir_cut

            βln_f = βln_data[1:ir_cut] .* r[1:ir_cut]  # rβln -> r²βln
            βln_r = r[eachindex(βln_data)]
            βln_interpolator = Interpolation.construct_interpolator(
                βln_r, βln_f, interpolation_method
            )

            return KleinmanBylanderProjector{RealSpace,Numerical}(
                βln_r, βln_f, βln_interpolator, n, l, j
            )
        end
    end
    β = OffsetVector(β, 0:lmax)  # Units are 1/√a₀ * a₀²
    # Extract the blocks from `upf_file.nonlocal.dij` corresponding to each angular momentum
    D = map(1:(length(cumul_nβ) - 1)) do i
        return collect(
            upf_file.nonlocal.dij[
                (cumul_nβ[i] + 1):cumul_nβ[i + 1], (cumul_nβ[i] + 1):cumul_nβ[i + 1]
            ],
        )
    end
    D = OffsetVector(D, 0:lmax) ./ 2  # Ry -> Ha
    Vnl = NonLocalPotential{RealSpace,Numerical,eltype(eltype(β))}(β, D)

    if isnothing(upf_file.nonlocal.augmentation)
        aug = nothing
    else
        q = OffsetVector(
            map(
                i -> collect(
                    upf_file.nonlocal.augmentation.q[
                        (cumul_nβ[i] + 1):cumul_nβ[i + 1],
                        (cumul_nβ[i] + 1):cumul_nβ[i + 1],
                    ],
                ),
                1:(length(cumul_nβ) - 1),
            ),
            0:(lmax),
        )
        if upf_file.nonlocal.augmentation.q_with_l
            Q = _upf_construct_augmentation_q_with_l(upf_file, interpolation_method)
        elseif upf_file.nonlocal.augmentation.nqf > 0
            Q = _upf_construct_augmentation_qfcoef(upf_file, interpolation_method)
        else
            error("q_with_l == false and nqf == 0, unsure what to do...")
        end
        aug = Augmentation{RealSpace,Numerical}(Q, q)
    end

    ## Core charge density
    if isnothing(upf_file.nlcc)
        ρcore = nothing
    else
        ρcore_r = r[eachindex(upf_file.nlcc)]
        ρcore_f = ρcore_r .^ 2 .* upf_file.nlcc
        ρcore_interpolator = Interpolation.construct_interpolator(
            ρcore_r, ρcore_f, interpolation_method
        )
        ρcore = ChargeDensity{RealSpace,Numerical}(ρcore_r, ρcore_f, ρcore_interpolator)
    end

    ## Valence charge density
    if isnothing(upf_file.rhoatom)
        ρval = nothing
    else
        ρval_r = r[eachindex(upf_file.rhoatom)]
        ρval_f = upf_file.rhoatom ./ 4π
        ρval_interpolator = Interpolation.construct_interpolator(
            ρval_r, ρval_f, interpolation_method
        )
        ρval = ChargeDensity{RealSpace,Numerical}(ρval_r, ρval_f, ρval_interpolator)
    end

    ## Pseudo-atomic states
    if !isnothing(upf_file.pswfc)
        # Collect the indices in upf_file.pswfc for projectors at each angular momentum
        iχ_upf = map(0:lmax) do l
            return filter(i -> upf_file.pswfc[i].l == l, eachindex(upf_file.pswfc))
        end
        iχ_upf = OffsetVector(iχ_upf, 0:lmax)

        # UPFs store the wavefunctions multiplied by the radial grid, so we multiply again
        # by r (to get r²χ) for consistency.
        χ = map(0:lmax) do l
            map(iχ_upf[l]) do n
                χln_data = upf_file.pswfc[n].chi
                χln_r = r[eachindex(χln_data)]
                χln_f = χln_data .* χln_r
                χln_interpolator = Interpolation.construct_interpolator(
                    χln_r, χln_f, interpolation_method
                )
                return StateProjector{RealSpace,Numerical}(
                    χln_r, χln_f, χln_interpolator, n, l, j
                )
            end
        end
        states = OffsetVector(χ, 0:lmax)
    else
        states = OffsetVector(Vector{Nothing}[], 0:-1)
    end

    return AtomicPotential(identifier, symbol, Vloc, Vnl, ρval, ρcore, states, aug)
end

function AtomicPotential(hgh_file::PseudoPotentialIO.HghFile)
    lmax = hgh_file.lmax
    symbol = Symbol(PseudoPotentialIO.element(hgh_file).symbol)
    Zval = sum(hgh_file.zion)

    cloc = hgh_file.cloc
    length(cloc) <= 4 || error("length(cloc) > 4 not supported.")
    if length(cloc) < 4
        n_extra = 4 - length(cloc)
        cloc = [cloc; zeros(Float64, n_extra)]
    end

    Vloc = HghLocalPotential{RealSpace,Analytical}(hgh_file.rloc, cloc, Zval)

    β = map(0:lmax) do l
        map(1:size(hgh_file.h[l + 1], 1)) do n
            return HghKleinmanBylanderProjector{RealSpace,Analytical}(hgh_file.rp[l + 1], n, l)
        end
    end
    β = OffsetVector(β, 0:lmax)  # Ha / √a₀
    D = OffsetVector(hgh_file.h, 0:lmax)  # 1 / Ha
    Vnl = NonLocalPotential{
        RealSpace,Analytical,HghKleinmanBylanderProjector{RealSpace,Analytical}
    }(
        β, D
    )

    ρval = nothing
    ρcore = nothing
    states = OffsetVector(Vector{Nothing}[], 0:-1)
    aug = nothing

    return AtomicPotential(hgh_file.identifier, symbol, Vloc, Vnl, ρval, ρcore, states, aug)
end
