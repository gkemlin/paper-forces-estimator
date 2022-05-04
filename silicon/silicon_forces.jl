a = 10.26  # Silicon lattice constant in Bohr
lattice = a / 2 * [[0 1 1.];
                   [1 0 1.];
                   [1 1 0.]]
Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))
atoms = [Si, Si]
positions = [ones(3)/8 + [0.24, -0.33, 0.12] / 20, -ones(3)/8]

model = model_LDA(lattice, atoms, positions)
kgrid = [2, 2, 2]  # k-point grid (Regular Monkhorst-Pack grid)
ss = 4
fft_size_ref = compute_fft_size(model, Ecut_ref; supersampling=ss)
basis_ref = PlaneWaveBasis(model; Ecut=Ecut_ref, kgrid, fft_size=fft_size_ref)
tol = 1e-13

scfres_ref = self_consistent_field(basis_ref; tol,
                                   is_converged=DFTK.ScfConvergenceDensity(tol))
f_ref = compute_forces(scfres_ref)
ψ_ref, occupation = DFTK.select_occupied_orbitals(basis_ref, scfres_ref.ψ, scfres_ref.occupation)

forces_list = []
forces_err_list = []
forces_res_list = []
forces_newton_list = []
diff_forces_err_list = []
diff_forces_res_list = []
diff_forces_newton_list = []
err_list = []
res_list = []
Msqrt_err_list = []
Msqrt_res_list = []
newton_list = []
diff_res_list = []
diff_newton_list = []

for Ecut in Ecut_list
    println("---------------------")
    println("Ecut = $(Ecut)")
    fft_size = compute_fft_size(model, Ecut; supersampling=ss)
    basis = PlaneWaveBasis(model; Ecut, kgrid, fft_size)
    scfres = self_consistent_field(basis; tol,
                                   is_converged=DFTK.ScfConvergenceDensity(tol))
    ψ, _ = DFTK.select_occupied_orbitals(basis_ref, scfres.ψ, scfres.occupation)
    ψr = DFTK.transfer_blochwave(scfres.ψ, basis, basis_ref)
    ρr = compute_density(basis_ref, ψr, scfres.occupation)
    _, ham = energy_hamiltonian(basis_ref, ψr, scfres.occupation; ρ=ρr)
    f = compute_forces(scfres)

    res = DFTK.compute_projected_gradient(basis_ref, ψr, scfres.occupation)
    res, _ = DFTK.select_occupied_orbitals(basis_ref, res, scfres.occupation)
    ψr, _ = DFTK.select_occupied_orbitals(basis_ref, ψr, scfres.occupation)
    err = compute_error(basis_ref, ψr, ψ_ref)

    P = [PreconditionerTPA(basis_ref, kpt) for kpt in basis_ref.kpoints]
    for (ik, ψk) in enumerate(ψr)
        precondprep!(P[ik], ψk)
    end
    Mres = apply_metric(ψr, P, res, apply_inv_M)
    Msqrt_res = apply_metric(ψr, P, res, apply_inv_sqrt_M)
    Msqrt_err = apply_metric(ψr, P, err, apply_sqrt_M)

    ## Rayleigh coefficients
    Λ = map(enumerate(ψr)) do (ik, ψk)
        Hk = ham.blocks[ik]
        Hψk = Hk * ψk
        ψk'Hψk
    end

    ## schur 1st step
    resLF = DFTK.transfer_blochwave(res, basis_ref, basis)
    resHF = res - DFTK.transfer_blochwave(resLF, basis, basis_ref)
    e2 = apply_metric(ψr, P, resHF, apply_inv_T)
    ΩpKe2 = apply_Ω(e2, ψr, ham, Λ) .+ apply_K(basis_ref, e2, ψr, ρr, occupation)
    ΩpKe2 = DFTK.transfer_blochwave(ΩpKe2, basis_ref, basis)
    rhs = resLF - ΩpKe2

    e1 = DFTK.solve_ΩplusK(basis, ψ, rhs, occupation; tol_cg=tol).δψ
    e1 = DFTK.transfer_blochwave(e1, basis, basis_ref)
    e_newton = e1 + Mres

    f_err = δforces(basis_ref, occupation, ψr, proj_tangent(err, ψr), ρr)
    f_res = δforces(basis_ref, occupation, ψr, Mres, ρr)
    f_newton = δforces(basis_ref, occupation, ψr, e_newton, ρr)

    # update lists
    append!(forces_list, norm(f-f_ref))
    append!(forces_err_list, norm(f_err))
    append!(forces_res_list, norm(f_res))
    append!(forces_newton_list, norm(f_newton))
    append!(diff_forces_err_list, norm(f-f_ref-f_err))
    append!(diff_forces_res_list, norm(f-f_ref-f_res))
    append!(diff_forces_newton_list, norm(f-f_ref-f_newton))
    append!(err_list, norm(err))
    append!(res_list, norm(Mres))
    append!(Msqrt_err_list, norm(Msqrt_err))
    append!(Msqrt_res_list, norm(Msqrt_res))
    append!(newton_list, norm(e_newton))
    append!(diff_res_list, norm(compute_error(basis_ref, err, Mres)))
    append!(diff_newton_list, norm(compute_error(basis_ref, err, e_newton)))
end

h5open("forces_silicon.h5", "w") do file
    println("writing h5 file")
    T = eltype(basis_ref)
    file["E_ref"] = scfres_ref.energies.total
    file["f_ref"] = Array{T}(hcat(f_ref...))
    file["Ecut_ref"] = Ecut_ref
    file["Ecut_list"] = collect(Ecut_list)
    file["forces_list"] = Array{T}(forces_list)
    file["forces_err_list"] = Array{T}(forces_err_list)
    file["forces_res_list"] = Array{T}(forces_res_list)
    file["forces_newton_list"] = Array{T}(forces_newton_list)
    file["diff_forces_err_list"] = Array{T}(diff_forces_err_list)
    file["diff_forces_res_list"] = Array{T}(diff_forces_res_list)
    file["diff_forces_newton_list"] = Array{T}(diff_forces_newton_list)
    file["err_list"] = Array{T}(err_list)
    file["res_list"] = Array{T}(res_list)
    file["Msqrt_err_list"] = Array{T}(Msqrt_err_list)
    file["Msqrt_res_list"] = Array{T}(Msqrt_res_list)
    file["newton_list"] = Array{T}(newton_list)
    file["diff_res_list"] = Array{T}(diff_res_list)
    file["diff_newton_list"] = Array{T}(diff_newton_list)
end
