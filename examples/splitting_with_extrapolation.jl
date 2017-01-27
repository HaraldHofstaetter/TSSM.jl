function gen_interpolation_matrix(x::Vector{Float64}, z::Vector{Float64})
    n = length(x)
    m = length(z)    
    L = zeros(n, m)
    for j=1:m
        for i=1:n
            L[i,j] = prod([(z[j]-x[k])/(x[i]-x[k]) for k=1:i-1])*prod([(z[j]-x[k])/(x[i]-x[k]) for k=i+1:n])
        end
    end
    L
end


function gen_extrapolated_wf!(psi::WaveFunction, L::Vector{Float64}, psi_back::Vector{WaveFunction}, ptr::Int; 
                              order::Int=length(psi_back)-1)
    n = length(psi_back)
    set!(psi, 0)
    for j=1:order+1
        k = mod(j-order+ptr-2, n)+1
        axpy!(psi, psi_back[k], L[j])
    end
end



function gen_starting_values(psi::WaveFunction, dt::Number, N::Int, a::Vector{Float64}, b::Vector{Float64};
    psi_ex::WaveFunction=wave_function(psi.m), rhs::WaveFunction=wave_function(psi.m), nonlinear_potential::Bool=false)
    s = length(a)
    z = cumsum(a)
    psi_back = WaveFunction[wave_function(m) for j=1:N]
    copy!(psi_back[1], psi)   
    t0 = get_time(psi)
    for K=2:N
        copy!(psi, psi_back[1])
        for J=2:K
            L = gen_interpolation_matrix(collect(0.0:(K-2.0)), J-2.0+z)
            for j = 1:s
                if a[j]!=0.0
                    propagate_A!(psi, a[j]*dt)            
                end
                if b[j]!=0.0
                    gen_extrapolated_wf!(psi_ex, L[:,j], psi_back[1:K-1], 2*K-2)
                    set_time!(psi_ex, get_time(psi))
                    if nonlinear_potential
                        gen_nonlinear_potential!(rhs, psi_ex)
                        to_real_space!(psi)
                        to_real_space!(rhs)
                        u = get_data(psi, true)
                        f = get_data(rhs, true)
                        u[:] .*= exp(b[j]*dt*f)
                    else
                        gen_rhs!(rhs, psi_ex)
                        axpy!(psi, rhs, b[j]*dt)
                    end
                end
            end            
            copy!(psi_back[J], psi)
            to_real_space!(psi_back[J])
        end
    end 
    psi_back
end



type SplittingWithExtrapolationEquidistantTimeStepperIterator
    psi::WaveFunction
    t0::Real
    tend::Real
    dt::Real
    steps::Int
    a::Vector{Float64}
    b::Vector{Float64}
    L::Matrix{Float64}
    ptr::Int
    psi_back::Vector{WaveFunction} # storage for previous solution values
    psi_ex::WaveFunction  # storage for the extrapolated values
    rhs::WaveFunction     # storage for the righthandside of the B equation
    nonlinear_potential::Bool
end


function splitting_with_extrapolation_equidistant_time_stepper(psi::WaveFunction, 
         t0::Real, tend::Real, dt::Real, order::Int, a::Vector{Float64}, b::Vector{Float64}; steps::Int=-1, nonlinear_potential::Bool=false) 
    m = psi.m
    psi_ex = wave_function(m)
    rhs = wave_function(m)
    set_propagate_time_together_with_A!(m, true)
    set_time!(psi, t0)
    psi_back = gen_starting_values(psi, dt, order+1, a, b, psi_ex=psi_ex, rhs=rhs, nonlinear_potential=nonlinear_potential)
    copy!(psi, psi_back[order+1])
    z = cumsum(a)
    L = gen_interpolation_matrix(collect((-order-0.0):0.0), z)
    ptr = order+1
    SplittingWithExtrapolationEquidistantTimeStepperIterator(psi, t0, tend, dt, steps, a, b, L, ptr, psi_back, psi_ex, rhs, nonlinear_potential)
end


Base.start(tsi::SplittingWithExtrapolationEquidistantTimeStepperIterator) = (tsi.steps<0 ? 
    tsi.t0 + (length(tsi.psi_back)-1)*tsi.dt : length(tsi.psi_back)-1)


Base.done(tsi::SplittingWithExtrapolationEquidistantTimeStepperIterator, t) = (tsi.steps<0 ? t >= tsi.tend : t>=tsi.steps )


function Base.next(tsi::SplittingWithExtrapolationEquidistantTimeStepperIterator, t)
    for j = 1:length(tsi.a)
        if tsi.a[j]!=0.0
            propagate_A!(tsi.psi, tsi.a[j]*tsi.dt)            
        end
        if tsi.b[j]!=0.0
            gen_extrapolated_wf!(tsi.psi_ex, tsi.L[:,j], tsi.psi_back, tsi.ptr)
            set_time!(tsi.psi_ex, get_time(tsi.psi))
            if tsi.nonlinear_potential
                gen_nonlinear_potential!(tsi.rhs, tsi.psi_ex)
                to_real_space!(tsi.psi)
                to_real_space!(tsi.rhs)
                u = get_data(tsi.psi, true)
                f = get_data(tsi.rhs, true)
                u[:] .*= exp(tsi.b[j]*tsi.dt*f)
            else
                gen_rhs!(tsi.rhs, tsi.psi_ex)
                axpy!(tsi.psi, tsi.rhs, tsi.b[j]*tsi.dt)
            end
        end
    end      
    n = length(tsi.psi_back)
    tsi.ptr = mod(tsi.ptr, n) + 1
    copy!(tsi.psi_back[tsi.ptr], tsi.psi)
    to_real_space!(tsi.psi_back[tsi.ptr])
    if tsi.steps<0 
        t1 = t + tsi.dt < tsi.tend ? t + tsi.dt : tsi.tend
    else
        t1 = t+1
    end
    return t1, t1
end


function global_orders(psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real, order::Int, a::Vector{Float64}, b::Vector{Float64}; 
                        operator_sequence="AB", rows=8, nonlinear_potential::Bool=false)
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
        for t in splitting_with_extrapolation_equidistant_time_stepper(psi, t0, tend, dt1, order, a, b, steps=steps,
                nonlinear_potential=nonlinear_potential)
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





type SplittingWithExtrapolationAdaptiveTimeStepperIterator
    psi::WaveFunction
    t::Float64
    t0::Float64
    tend::Float64
    dt::Float64
    tol::Float64
    dt_old::Float64
    err_old::Float64
    order::Int
    order1::Int
    a1::Vector{Float64}
    b1::Vector{Float64}
    z1::Vector{Float64}
    order2::Int
    a2::Vector{Float64}
    b2::Vector{Float64}
    z2::Vector{Float64}
    N::Int
    ptr::Int
    psi_back::Vector{WaveFunction} # storage for previous solution values
    t_back::Vector{Float64}  # previous times
    psi_ex::WaveFunction  # storage for the extrapolated values
    rhs::WaveFunction     # storage for the righthandside of the B equation
    psi2::WaveFunction
    psi0::WaveFunction    
    nonlinear_potential::Bool
end



function splitting_with_extrapolation_adaptive_time_stepper(psi::WaveFunction, 
        t0::Float64, tend::Float64,  dt::Float64, tol::Float64,order::Int,
        order1::Int, a1::Vector{Float64}, b1::Vector{Float64},
        order2::Int, a2::Vector{Float64}, b2::Vector{Float64}; nonlinear_potential::Bool=false)
    m =psi.m
    psi_ex = wave_function(m)
    rhs = wave_function(m)
    psi0 = wave_function(m)
    psi2 = wave_function(m)
    
    set_propagate_time_together_with_A!(m, true)
    set_time!(psi, t0)
    N = max(order1, order2)+1
    psi_back = gen_starting_values(psi, dt, N, a1, b1, psi_ex=psi_ex, rhs=rhs, nonlinear_potential=nonlinear_potential)
    t_back = [dt*k for k=0:N-1] 
    t = t_back[N]
    ptr = order+1
    copy!(psi, psi_back[ptr])
    z1 = cumsum(a1)
    z2 = cumsum(a2)
    SplittingWithExtrapolationAdaptiveTimeStepperIterator(psi,t,t0,tend,dt,tol,-1.0,-1.0,order,
        order1,a1,b1,z1,order2,a2,b2,z2, N, ptr, psi_back, t_back, psi_ex, rhs, psi2, psi0, nonlinear_potential)
end


Base.start(tsi::SplittingWithExtrapolationAdaptiveTimeStepperIterator) = tsi.t 

function Base.done(tsi::SplittingWithExtrapolationAdaptiveTimeStepperIterator, state_t::Float64) 
   tsi.t >= tsi.tend
end  

function Base.next(tsi::SplittingWithExtrapolationAdaptiveTimeStepperIterator, state_t::Float64 ) 
    const facmin = 0.25
    const facmax = 4.0
    const fac = 0.9
    const beta1=0.25
    const beta2=0.25
    const alpha2=0.25
    dt = tsi.dt
    dt0 = dt
    dt_old=dt
    copy!(tsi.psi0, tsi.psi)
    err = 2.0 #error/tol
    while err>=1.0
        dt = min(dt, tsi.tend-tsi.t)
        dt0 = dt
        tt = Float64[ tsi.t_back[mod(j-tsi.order1+tsi.ptr-2, tsi.N)+1] for j=1:tsi.order1+1]
        L1 = gen_interpolation_matrix(tt, tsi.t + dt*tsi.z1)
        L2 = gen_interpolation_matrix(tt, tsi.t + dt*tsi.z2)
        
        #step_embedded!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.scheme2, tsi.operator_sequence)
        k_start2 = -1
        start2_with_A = false
        
        for k = 1:length(tsi.a1)
            if k_start2<0.0 && tsi.a1[k]!=tsi.a2[k]
                k_start2 = k                
                start2_with_A = true                
                copy!(tsi.psi2, tsi.psi)
            end
            if tsi.a1[k]!=0.0
                propagate_A!(tsi.psi, tsi.a1[k]*tsi.dt)            
            end
            if k_start2<0.0 && b1[k]!=tsi.b2[k]
                k_start2 = k                
                start2_with_A = false
                copy!(tsi.psi2, tsi.psi)
            end
            if tsi.b1[k]!=0.0
                gen_extrapolated_wf!(tsi.psi_ex, L1[:,k], tsi.psi_back, tsi.ptr, order=tsi.order1)
                set_time!(tsi.psi_ex, get_time(tsi.psi))
                if tsi.nonlinear_potential
                    gen_nonlinear_potential!(tsi.rhs, tsi.psi_ex)
                    to_real_space!(tsi.psi)
                    to_real_space!(tsi.rhs)
                    u = get_data(tsi.psi, true)
                    f = get_data(tsi.rhs, true)
                    u[:] .*= exp(tsi.b1[k]*tsi.dt*f)                    
                else
                    gen_rhs!(tsi.rhs, tsi.psi_ex)
                    axpy!(tsi.psi, tsi.rhs, tsi.b1[k]*tsi.dt)
                end
            end
        end     
        
        if k_start2<0
            k_start2 = length(tsi.a1)+1
            start2_with_A = true
            copy!(tsi.psi2, tsi.psi)
        end     
        
        for k = k_start2:length(tsi.a2)
            if tsi.a2[k]!=0.0 &&(k!=k_start2||start2_with_A)
                propagate_A!(tsi.psi2, tsi.a2[k]*tsi.dt)            
            end
            if tsi.b2[k]!=0.0
                gen_extrapolated_wf!(tsi.psi_ex, L2[:,k], tsi.psi_back, tsi.ptr, order=tsi.order1)
                set_time!(tsi.psi_ex, get_time(tsi.psi2))
                if tsi.nonlinear_potential
                    gen_nonlinear_potential!(tsi.rhs, tsi.psi_ex)
                    to_real_space!(tsi.psi2)
                    to_real_space!(tsi.rhs)
                    u = get_data(tsi.psi2, true)
                    f = get_data(tsi.rhs, true)
                    u[:] .*= exp(tsi.b2[k]*tsi.dt*f)                                       
                else               
                    gen_rhs!(tsi.rhs, tsi.psi_ex)
                    axpy!(tsi.psi2, tsi.rhs, tsi.b2[k]*tsi.dt)
                end
            end
        end        
            
        err = distance(tsi.psi, tsi.psi2)/tsi.tol
        if tsi.dt_old<0.0
           dt = dt*min(facmax, max(facmin, fac*(1.0/err)^(1.0/(tsi.order+1))))
        else
            dt = dt*min(facmax, max(facmin, fac*((1.0/err)^(beta1/(tsi.order+1))*(1.0/tsi.err_old)^(beta2/(tsi.order+1))*(dt/tsi.dt_old)^(-alpha2))))         
        end
        if err>=1.0
           copy!(tsi.psi, tsi.psi0)
            @printf("t=%17.9e  err=%17.8e  dt=%17.8e  rejected...\n", Float64(tsi.t), Float64(err), Float64(dt))
        end   
    end
    
    tsi.t += dt0
    tsi.dt = dt
    tsi.dt_old = dt_old
    tsi.err_old = err
    
    tsi.ptr = mod(tsi.ptr, tsi.N) + 1
    copy!(tsi.psi_back[tsi.ptr], tsi.psi)
    to_real_space!(tsi.psi_back[tsi.ptr])   
    tsi.t_back[tsi.ptr] = tsi.t      
    
    tsi.t, tsi.t 
end

