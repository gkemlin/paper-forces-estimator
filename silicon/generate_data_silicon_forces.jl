using HDF5
using CSV, Tables

h5open("forces_silicon_ss4.h5", "r") do file
    # read h5
    E_ref = read(file["E_ref"])
    println(E_ref)
    f_ref = read(file["f_ref"])
    display(f_ref)
    Ecut_ref = read(file["Ecut_ref"])
    Ecut_list = read(file["Ecut_list"])
    forces_list = read(file["forces_list"])
    forces_err_list = read(file["forces_err_list"])
    forces_res_list = read(file["forces_res_list"])
    forces_newton_list = read(file["forces_newton_list"])
    diff_forces_err_list = read(file["diff_forces_err_list"])
    diff_forces_res_list = read(file["diff_forces_res_list"])
    diff_forces_newton_list = read(file["diff_forces_newton_list"])
    err_list = read(file["err_list"])
    res_list = read(file["res_list"])
    Msqrt_err_list = read(file["Msqrt_err_list"])
    Msqrt_res_list = read(file["Msqrt_res_list"])

    # write CSV
    CSV.write("forces_silicon.csv", Tables.table([Ecut_list forces_list forces_err_list forces_res_list forces_newton_list]), writeheader=false)
    CSV.write("error_silicon.csv", Tables.table([Ecut_list Msqrt_err_list Msqrt_res_list]), writeheader=false)
    CSV.write("diff_forces_silicon.csv", Tables.table([Ecut_list diff_forces_err_list diff_forces_res_list diff_forces_newton_list]), writeheader=false)
end
