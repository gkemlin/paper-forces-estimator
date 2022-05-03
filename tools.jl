using LinearMaps
using IterativeSolvers
using Statistics

############################# ERROR AND RESIDUAL ###############################

# compute the error on the orbitals by aligning the eigenvectors
# this is done by solving min |ϕ - ψ*U| for U unitary matrix of size NxN
# whose solution is U = M(M^*M)^-1/2 where M = ψ^*ϕ
function compute_error(basis, ϕ, ψ)

    # necessary quantites
    Nk = length(basis.kpoints)

    # compute error
    err = similar(ϕ)
    for ik = 1:Nk
        ϕk = ϕ[ik]
        ψk = ψ[ik]
        # compute overlap matrix
        M = ψk'ϕk
        U = M*(M'M)^(-1/2)
        err[ik] = ϕk - ψk*U
    end
    err
end

############################## CHANGES OF NORMS ################################

## T = -1/2 Δ + t

# applies per kpoint and per band

apply_T(φk, Pk, δφnk, n) = (Pk.mean_kin[n] .+ Pk.kin) .* δφnk
apply_sqrt_T(φk, Pk, δφnk, n) = sqrt.(Pk.mean_kin[n] .+ Pk.kin) .* δφnk

apply_inv_T(φk, Pk, δφnk, n) = δφnk ./ (Pk.mean_kin[n] .+ Pk.kin)
apply_inv_sqrt_T(φk, Pk, δφnk, n) = δφnk ./ sqrt.(Pk.mean_kin[n] .+ Pk.kin)

function apply_M(φk, Pk, δφnk, n)
    δφnk = proj_tangent_kpt(δφnk, φk)
    δφnk = apply_sqrt_T(φk, Pk, δφnk, n)
    δφnk = proj_tangent_kpt(δφnk, φk)
    δφnk = apply_sqrt_T(φk, Pk, δφnk, n)
    δφnk = proj_tangent_kpt(δφnk, φk)
end

function apply_sqrt_M(φk, Pk, δφnk, n)
    δφnk = proj_tangent_kpt(δφnk, φk)
    δφnk = apply_sqrt_T(φk, Pk, δφnk, n)
    δφnk = proj_tangent_kpt(δφnk, φk)
end

function apply_inv_M(φk, Pk, δφnk, n)
    δφnk = proj_tangent_kpt(δφnk, φk)
    function op(x)
        apply_M(φk, Pk, x, n)
    end
    function f_ldiv!(x, y)
        x .= proj_tangent_kpt(y, φk)
        x .= apply_inv_T(φk, Pk, x, n)
        proj_tangent_kpt!(x, φk)
    end
    J = LinearMap{eltype(φk)}(op, size(δφnk, 1))
    δφnk = cg(J, δφnk, Pl=FunctionPreconditioner(f_ldiv!),
              verbose=false, reltol=0, abstol=1e-15)
    proj_tangent_kpt(δφnk, φk)
end

function apply_inv_sqrt_M(φk, Pk, δφnk, n)
    δφnk = proj_tangent_kpt(δφnk, φk)
    function op(x)
        apply_sqrt_M(φk, Pk, x, n)
    end
    function f_ldiv!(x, y)
        x .= proj_tangent_kpt(y, φk)
        x .= apply_inv_sqrt_T(φk, Pk, x, n)
        proj_tangent_kpt!(x, φk)
    end
    J = LinearMap{eltype(φk)}(op, size(δφnk, 1))
    δφnk = cg(J, δφnk, Pl=FunctionPreconditioner(f_ldiv!),
              verbose=false, reltol=0, abstol=1e-15)
    proj_tangent_kpt(δφnk, φk)
end

# applies for full orbitals

# /!\ A(φk, Pk, δφnk, n) should be one of the eight functions above /!\
function apply_metric(φ, P, δφ, A::Function)
    map(enumerate(δφ)) do (ik, δφk)
        Aδφk = similar(δφk)
        φk = φ[ik]
        for n = 1:size(δφk,2)
            Aδφk[:,n] = A(φk, P[ik], δφk[:,n], n)
        end
        Aδφk
    end
end

"""
Extension of Diagonal to nonlocal operators
"""
function LinearAlgebra.Diagonal(opnl::DFTK.NonlocalOperator)
    [dot(p, opnl.D * p) for p in eachrow(opnl.P)]
end

"""
Compute the average of the local potential and the nonlocal potential
"""
function compute_avg(basis::PlaneWaveBasis, H::Hamiltonian)

    # average of the local part of the potential of the Hamiltonian
    avg_local_pot = mean(DFTK.total_local_potential(H))

    # adding the average on the nonlocal part of the potential depending on the
    # k point
    total_pot_avg = []
    for (ik, kpt) in enumerate(basis.kpoints)
        non_local_op = [op for op in H.blocks[ik].operators
                        if (op isa DFTK.NonlocalOperator)][1]
        avg_non_local_op = Diagonal(non_local_op)

        # compute potential average if used in the perturbation
        total_pot_avgk = avg_local_pot .+ avg_non_local_op
        push!(total_pot_avg, total_pot_avgk)
    end
    total_pot_avg
end
