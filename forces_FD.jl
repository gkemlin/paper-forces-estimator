using ForwardDiff

function δforces(basis, occupation, ψ, δψ, ρ)
    δρ = DFTK.compute_δρ(basis, ψ, δψ, occupation)
    function f(ε)
        compute_forces(basis, ψ .+ ε .* δψ, occupation; ρ=ρ .+ ε .* δρ)
    end
    ForwardDiff.derivative(f, 0)
end
