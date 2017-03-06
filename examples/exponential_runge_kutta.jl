type ExponentialRungeKuttaEquidistantTimeStepperIterator
    psi::WaveFunction
    t0::Real
    tend::Real
    dt::Real
    steps::Int
    s::Int
    G::WaveFunction
    U::Vector{WaveFunction}
end


function exponential_runge_kutta_equidistant_time_stepper(psi::WaveFunction, 
         t0::Real, tend::Real, dt::Real; steps::Int=-1,)
    m = psi.m
    set_time!(psi, t0)
    s = 4
    G = wave_function(m)
    U = [wave_function(m) for j=1:4]
    ExponentialRungeKuttaEquidistantTimeStepperIterator(psi, t0, tend, dt, steps, s, G, U)
end

Base.start(tsi::ExponentialRungeKuttaEquidistantTimeStepperIterator) = (tsi.steps<0 ? tsi.t0 : 0)


Base.done(tsi::ExponentialRungeKuttaEquidistantTimeStepperIterator, t) = (tsi.steps<0 ? t >= tsi.tend : t>=tsi.steps )


function ERK4_step!(psi::WaveFunction, dt::Number; 
    G::WaveFunction=wave_function(psi.m),
    U::Vector{WaveFunction}=[wave_function(psi.m) for j=1:4])
    
    s = 4
    c = [0.0, 0.5, 0.5, 1.0]*dt
    t = get_time(psi)
    for j=1:s
        copy!(U[j], psi)
        set_time!(U[j], t+c[j])
    end
    set_time!(psi, t+dt)
    
    gen_rhs!(G, U[1])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, U[2], c[2], 1, 0.5*dt)
    add_phi_A!(G, U[3], c[3], 1, 0.5*dt)
    add_phi_A!(G, U[3], c[3], 2, -1.0*dt)
    add_phi_A!(G, U[4], c[4], 1, 1.0*dt)
    add_phi_A!(G, U[4], c[4], 2, -2.0*dt)
    add_phi_A!(G, psi, dt, 1, 1.0*dt)
    add_phi_A!(G, psi, dt, 2, -3.0*dt)
    add_phi_A!(G, psi, dt, 3, 4.0*dt)
    
    gen_rhs!(G, U[2])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, U[3], c[3], 2, 1.0*dt)
    add_phi_A!(G, psi, dt, 2, 2.0*dt)
    add_phi_A!(G, psi, dt, 3, -4.0*dt)
    
    gen_rhs!(G, U[3])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, U[4], c[4], 2, 2.0*dt)
    add_phi_A!(G, psi, dt, 2, 2.0*dt)
    add_phi_A!(G, psi, dt, 3, -4.0*dt)

    gen_rhs!(G, U[4])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, psi, dt, 2, -1.0*dt)
    add_phi_A!(G, psi, dt, 3, 4.0*dt)
    psi
end


function Base.next(tsi::ExponentialRungeKuttaEquidistantTimeStepperIterator, t)
    ERK4_step!(tsi.psi, tsi.dt, G=tsi.G, U=tsi.U)    
    if tsi.steps<0 
        t1 = t + tsi.dt < tsi.tend ? t + tsi.dt : tsi.tend
    else
        t1 = t+1
    end
    return t1, t1
 end
 
 
 function global_orders(psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real; rows=8)
    @assert psi.m==reference_solution.m
    tab = Array(Float64, rows, 3)

    wf_save_initial_value = clone(psi)
    copy!(wf_save_initial_value, psi)

    steps = Int(floor((tend-t0)/dt))
    dt1 = dt
    err_old = 0.0
    println("             dt         err      p")
    println("-----------------------------------")
    for row=1:rows
        for t in exponential_runge_kutta_equidistant_time_stepper(psi, t0, tend, dt1, steps=steps) 
        end    
        err = distance(psi, reference_solution)
        if (row==1) then
            @printf("%3i%12.3e%12.3e\n", row, Float64(dt1), Float64(err))
            tab[row,1] = dt1
            tab[row,2] = err
            tab[row,3] = 0 
        else
            p = log(err_old/err)/log(2.0);
            @printf("%3i%12.3e%12.3e%7.2f\n", row, Float64(dt1), Float64(err), Float64(p))
            tab[row,1] = dt1
            tab[row,2] = err
            tab[row,3] = p 
        end
        err_old = err
        dt1 = 0.5*dt1
        steps = 2*steps
        copy!(psi,wf_save_initial_value)
    end
end
