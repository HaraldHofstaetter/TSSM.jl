# Usage instruction:
# julia quantum_dot.jl groundstate # calculates groundstate at t0 = 0
# julia quantum_dot.jl ref [tol] # tol::Float64, calc ref solution with specified tol for AdaptiveAdamsLawson
# julia quantum_dot.jl {conv | stab} [method] [num_steps] # run experiment for equidistant method as listed below
# conv compares the result with a reference solution and prints the norm distance function
# stab calculates until t = 70 without external potential
#
# Parameters are changed as global variables otherwise
# Experiment results are being saved to subdirectories, please review those lines and mkdir before running!
# Saved data includes the wave_function as .h5 for ref and conv,
# the time, norm, energy (, stepsize) for every timestep in a .jld file.
# For conv, the error and the number of B_Calls is additionally saved.
using TSSM
using JLD
include("mctdhf1d.jl")
include("check.jl")
include("../examples/time_propagators.jl")

V1(x) = 1. / 32 * (x^2)
V2(x,y) = 1/sqrt((x-y)^2+1/16)
V(x,y) = V1(x) + V1(y) + V2(x,y)

const N = 4
const interval_bound = 20
const nx = Int(512*interval_bound/10)
const t0 = 0
const groundstate_tol=1e-9

println("N = ", N, ", nx = ", nx, " on [-", interval_bound, ",", interval_bound, "]")

function method_from_string(variant::String)
    if variant == "yoshida"
        g = [1.351207191959657634, -1.702414383919315268, 1.351207191959657634] # Yoshida
        a, b = get_coeffs_composition(g)
        method = SplittingRK4BMethod(a,b)
    end
    if variant == "suzuki"
        g = [1/(4-4^(1/3)),1/(4-4^(1/3)),-4^(1/3)/(4-4^(1/3)), 1/(4-4^(1/3)), 1/(4-4^(1/3))] # Suzuki
        a, b = get_coeffs_composition(g)
        method = SplittingRK4BMethod(a,b)
    end
    if variant == "strang"
        a, b = [0.5, 0.5], [1.0, 0.0]
        method = SplittingRK4BMethod(a, b, secondorder=true)
    end
    if variant == "krogstad"
        method = ExponentialRungeKutta(:krogstad)
    end
    if variant == "rk4"
        method = ExponentialRungeKutta(:rk4)
    end
    if variant == "rklawson"
        method = ExponentialRungeKutta(:lawson)
    end
    if variant == "strehmel_weiner"
        method = ExponentialRungeKutta(:strehmel_weiner)
    end
    if variant == "adamslawson3"
        method = ExponentialMultistep(3, version=2, iters=1, starting_method=ExponentialRungeKutta(:lawson), starting_subdivision=50)
    end
    if variant == "expmultistep3"
        method = ExponentialMultistep(3, version=1, iters=1, starting_method=ExponentialRungeKutta(:lawson), starting_subdivision=50)
    end
    if variant == "adamslawson5"
        method = ExponentialMultistep(5, version=2, iters=1, starting_method=ExponentialRungeKutta(:lawson), starting_subdivision=50)
        println("adamslawson5")
    end
    if variant == "expmultistep5"
        method = ExponentialMultistep(5, version=1, iters=1, starting_method=ExponentialRungeKutta(:lawson), starting_subdivision=50)
        println("exponential multistep5")
    end
    if variant == "adamslawson3_nocorr"
        method = ExponentialMultistep(3, version=2, iters=0, starting_method=ExponentialRungeKutta(:lawson), starting_subdivision=50)
        println("adamslawson3 without correction step")
    end
    if variant == "expmultistep3_nocorr"
        method = ExponentialMultistep(3, version=1, iters=0, starting_method=ExponentialRungeKutta(:lawson), starting_subdivision=50)
        println("exponential multistep 3 without correction step")
    end
    return method
end

function ref(tol::Float64=1e-5)
    E(t) = sin(2*t)
    tend = 12

    V1_t(x, t) = E(t)*x
    V_t(x,y,t) = E(t)*(x+y)
    
    m = MCTDHF1D(2, N, nx, -interval_bound, interval_bound, potential1=V1, potential1_t=V1_t, potential2=V2, spin_restricted=true);
    psi = wave_function(m);
    TSSM.load!(psi, "wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_t0_1e-9.h5")
    set_time!(psi, t0)

    times, energies, norms, stepsizes = propagate_adaptive!(psi, t0, tend, tol)

    TSSM.save(psi, "results_tolerance/wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_t"*string(tend)*"_"*string(tol)*".h5")
    JLD.save("results_tolerance/wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_"*string(tol)*".jld", "times", times, "energies", energies, "norms", norms, "stepsizes", stepsizes)
end

function conv(variant::String, steps::Int)
    #laser field
    E(t) = sin(2*t)
    tend = 12
    dt=tend/steps

    V1_t(x, t) = E(t)*x
    V_t(x,y,t) = E(t)*(x+y)

    m = MCTDHF1D(2, N, nx, -interval_bound, interval_bound, potential1=V1, potential1_t=V1_t, potential2=V2, spin_restricted=true)
    psi = wave_function(m)
    
    TSSM.load!(psi, "wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_t0_1e-9.h5")
    set_time!(psi, t0)
    
    psi_ref = wave_function(m);
    TSSM.load!(psi_ref, "results_tolerance/wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_t"*string(tend)*"_1.0e-9.h5")
    
    method = method_from_string(variant)

    times, energies, norms = propagate_equidistant!(method, psi, t0, dt, steps)

    err = distance(psi, psi_ref)
    println("err:", err)
    TSSM.save(psi, "results_conv/wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_t"*string(tend)*"_"*variant*"_"*string(steps)*".h5")
    JLD.save("results_conv/wfquantumdot_conv_f2_n"*string(N)*"_nx"*string(nx)*"_t"*string(tend)*"_"*variant*"_"*string(steps)*".jld", "times", times, "energies", energies, "dt", dt, "err", err, "B_calls", COUNT_B)
    
end

function stab(variant::String, steps::Int)
    #laser field

    E(t) = 0
    tend = 70
    dt=tend/steps

    V1_t(x, t) = E(t)*x
    V_t(x,y,t) = E(t)*(x+y)
    
    m = MCTDHF1D(2, N, nx, -interval_bound, interval_bound, potential1=V1, potential1_t=V1_t, potential2=V2, spin_restricted=true);
    psi = wave_function(m);

    TSSM.load!(psi, "wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_t0_1e-9.h5")
    set_time!(psi, t0)

    method = method_from_string(variant)
    
    times, energies, norms = propagate_equidistant!(method, psi, t0, dt, steps)

    JLD.save("results_stab/wfquantumdot_stab_f2_n"*string(N)*"_nx"*string(nx)*"_t"*string(tend)*"_"*variant*"_"*string(steps)*".jld", "times", times, "energies", energies, "dt", dt, "norms", norms)
end

function propagate_equidistant!(method, psi, t0, dt, steps)
    times = Vector{Float64}()
    energies = Vector{Float64}()
    norms = Vector{Float64}()
    step_counter = 1
    time0 = time()
    for (step, tsi) in EquidistantTimeStepper(method, psi, t0, dt, steps)
        n = norm(psi)
        E_pot = potential_energy(psi)
        E_kin = kinetic_energy(psi)
        E_tot = E_pot+E_kin
        if step_counter % 50 == 0 || step == steps
            @printf("%5i  %14.10f  %14.10f  %14.10f  %14.10f  %14.10f  %14.10f %6i %10.2f\n", 
                step_counter, get_time(psi), E_pot, E_kin, E_tot, dt, n, COUNT_B, time()-time0)
        end
        append!(times, get_time(psi))
        append!(energies, E_tot)
        append!(norms, n)
        step_counter += 1
    end
    return times, energies, norms
end

function propagate_adaptive!(psi, t0, tend, tol)
    method = AdaptiveAdamsLawson(6)
    times = Vector{Float64}()
    energies = Vector{Float64}()
    stepsizes = Vector{Float64}()
    norms = Vector{Float64}()
    step_counter = 1
    time0 = time()
    told = t0
    for (t, tsi) in AdaptiveTimeStepper(method, psi, t0, tend, tol, 0.00001)
        n = norm(psi)
        E_pot = potential_energy(psi)
        E_kin = kinetic_energy(psi)
        E_tot = E_pot+E_kin
        stepsize = t - told
        if step_counter % 50 == 0
            @printf("%5i  %14.10f  %14.10f  %14.10f  %14.10f  %14.10f  %14.10f  %10.2f\n", 
                step_counter, get_time(psi), n, E_pot, E_kin, E_tot, stepsize, time()-time0)
        end
        append!(times, get_time(psi))
        append!(energies, E_tot)
        append!(stepsizes, stepsize)
        append!(norms, n)
        told = t
        step_counter += 1
    end
    return times, energies, norms, stepsizes
end

if ARGS[1] == "groundstate"
    println("calculate groundstate with tolerance ", groundstate_tol)
    include("propagators.jl")
    m = MCTDHF1D(2, N, nx, -interval_bound, interval_bound, potential1=V1, potential2=V2, spin_restricted=true)
    psi = wave_function(m)
    groundstate!(psi, dt=0.1, max_iter=500000, output_step=100, tol=groundstate_tol)
    TSSM.save(psi, "wfquantumdot_f2_n"*string(N)*"_nx"*string(nx)*"_t0_1e-9.h5")
elseif ARGS[1] == "ref"
    println("calculate reference solution with tolerance ", ARGS[2])
    ref(parse(Float64, ARGS[2]))
elseif ARGS[1] == "conv"
    println("calculate solution using method ", ARGS[2], " with ", ARGS[3], " steps")
    conv(ARGS[2], parse(Int, ARGS[3]))
elseif ARGS[1] == "stab"
    stab(ARGS[2], parse(Int, ARGS[3]))
else
    throw(ErrorException("first cmdline argument is not 'ref', 'conv' or 'stab'"))
end
    
