using DFTK
using LinearMaps
using IterativeSolvers
using HDF5

import DFTK: proj_tangent, proj_tangent!, proj_tangent_kpt, proj_tangent_kpt!
import DFTK: pack_ψ, unpack_ψ, reinterpret_real, reinterpret_complex
import DFTK: precondprep!, FunctionPreconditioner
import DFTK: apply_K, apply_Ω

# the computations were run on a 32 cores cluster
using LinearAlgebra
using FFTW
BLAS.set_num_threads(32)
FFTW.set_num_threads(32)

include("tools.jl")
include("forces_FD.jl")

# reference Ecut and list of Ecut's for defining coarse grids
Ecut_ref = 125
Ecut_list = 5:5:90

@time begin
    # silicon computations
    cd("silicon/")
    include("silicon/silicon_forces.jl")
    include("silicon/generate_data_silicon_forces.jl")
    cd("../")

    # GaAs computations
    cd("GaAs/")
    include("GaAs/GaAs_forces.jl")
    include("GaAs/generate_data_GaAs_forces.jl")
    cd("../")

    # TiO2 computations
    cd("TiO2/")
    include("TiO2/TiO2_forces.jl")
    include("TiO2/generate_data_TiO2_forces.jl")
    cd("../")
end
