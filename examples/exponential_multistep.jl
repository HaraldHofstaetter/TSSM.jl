
function gen_exponential_multistep_starting_values(psi::WaveFunction, dt::Number, N::Int; 
          final_iteration=false, combine_first::Bool=false)    
    acc = wave_function(psi.m)
    psi0 = wave_function(psi.m)
    rhs_back = WaveFunction[wave_function(psi.m) for j=1:N]
    gen_rhs!(rhs_back[1], psi)
    t0 = get_time(psi)
    copy!(psi0, psi)
    for K=2:N
        for R=1:(final_iteration&&K==N?2:1)
            copy!(psi, psi0)            
            for J=2:K
                V  = Matrix{Float64}(inv(Rational{Int}[k^m//factorial(m) for k=-J+2:K-J, m=0:K-2])) 
                if combine_first
                    copy!(acc, rhs_back[J-1])
                    add_apply_A!(psi, acc, 1.0)
                    add_phi_A!(acc, psi, dt, 1, dt)
                    set_time!(psi, get_time(psi)+dt)
                else
                    propagate_A!(psi, dt)
                end
                for m=(combine_first?2:1):K-1                
                    set!(acc, 0.0)
                    for k=1:K-1
                        axpy!(acc, rhs_back[k], V[m,k])
                    end
                    add_phi_A!(acc, psi, dt, m, dt)        
                end 
                gen_rhs!(rhs_back[J], psi)   
            end
        end
    end
    rhs_back
end




type ExponentialMultistepEquidistantTimeStepperIterator
    psi::WaveFunction
    t0::Real
    tend::Real
    dt::Real
    steps::Int
    V1::Matrix{Float64}
    V2::Matrix{Float64}
    N::Int
    N1::Int
    N2::Int
    iter::Int    
    ptr::Int
    rhs_back::Vector{WaveFunction} # storage for previous solution values
    psi0::WaveFunction
    acc::WaveFunction   
    combine_first::Bool
end    


function exponential_multistep_equidistant_time_stepper(psi::WaveFunction, 
         t0::Real, tend::Real, dt::Real, N1::Int; 
         iter::Int=0, N2::Int=(iter>0?N1+1:0), steps::Int=-1,
         psi_back=nothing,  final_iteration=false, combine_first::Bool=true) 
    m = psi.m
    set_propagate_time_together_with_A!(m, true)
    set_time!(psi, t0)
    N = max(N1, N2)
    if psi_back!=nothing 
        rhs_back = psi_back
        @assert length(psi_back)==N
        for k=1:N
            gen_rhs!(rhs, rhs_back[k])
            copy!(rhs_back[k], rhs)
        end        
    else
        rhs_back = gen_exponential_multistep_starting_values(psi, dt, N, 
                       final_iteration=final_iteration, combine_first=combine_first)
    end    
    V1 = Matrix{Float64}(inv(Rational{Int}[n^m//factorial(m) for n=-N1+1:0, m=0:N1-1]))
    V2 = Matrix{Float64}(inv(Rational{Int}[n^m//factorial(m) for n=-N2+2:1, m=0:N2-1]))
    
    ptr = N
    acc = wave_function(m)
    psi0 = wave_function(m)
    ExponentialMultistepEquidistantTimeStepperIterator(psi, t0, tend, dt, steps, V1, V2, 
                 N, N1, N2, iter, ptr, rhs_back, psi0, acc, combine_first)
end


Base.start(tsi::ExponentialMultistepEquidistantTimeStepperIterator) = (tsi.steps<0 ? 
    tsi.t0 + (tsi.N-1)*tsi.dt : tsi.N-1)


Base.done(tsi::ExponentialMultistepEquidistantTimeStepperIterator, t) = (tsi.steps<0 ? t >= tsi.tend : t>=tsi.steps )


function Base.next(tsi::ExponentialMultistepEquidistantTimeStepperIterator, t)
    if tsi.iter>=1
        copy!(tsi.psi0, tsi.psi)
    end
    
    #predictor
    if tsi.combine_first
        k1 = mod(tsi.ptr-1, tsi.N)+1
        copy!(tsi.acc, tsi.rhs_back[k1])
        add_apply_A!(tsi.psi, tsi.acc, 1.0)
        add_phi_A!(tsi.acc, tsi.psi, tsi.dt, 1, tsi.dt)
        set_time!(tsi.psi, get_time(tsi.psi)+tsi.dt)
    else
        propagate_A!(tsi.psi, tsi.dt)
    end
    for m=(tsi.combine_first?2:1):tsi.N1
        set!(tsi.acc, 0.0)
        for k=1:tsi.N1
            k1 = mod(k-tsi.N1+tsi.ptr-1, tsi.N)+1
            axpy!(tsi.acc, tsi.rhs_back[k1], tsi.V1[m,k])
        end
        add_phi_A!(tsi.acc, tsi.psi, tsi.dt, m, tsi.dt)
    end
   
    tsi.ptr = mod(tsi.ptr, tsi.N) + 1
    gen_rhs!(tsi.rhs_back[tsi.ptr], tsi.psi)    
        
    for iter=1:tsi.iter
        copy!(tsi.psi, tsi.psi0)
        
        #corrector
        if tsi.combine_first
            k1 = mod(tsi.ptr-2, tsi.N)+1
            copy!(tsi.acc, tsi.rhs_back[k1])
            add_apply_A!(tsi.psi, tsi.acc, 1.0)
            add_phi_A!(tsi.acc, tsi.psi, tsi.dt, 1, tsi.dt)
            set_time!(tsi.psi, get_time(tsi.psi)+tsi.dt)
        else
            propagate_A!(tsi.psi, tsi.dt)
        end        
        for m=(tsi.combine_first?2:1):tsi.N2
            set!(tsi.acc, 0.0)
            for k=1:tsi.N2
                k1 = mod(k-tsi.N2+tsi.ptr-1, tsi.N)+1
                axpy!(tsi.acc, tsi.rhs_back[k1], tsi.V2[m,k])
            end
            add_phi_A!(tsi.acc, tsi.psi, tsi.dt, m, tsi.dt)
        end        
        
        gen_rhs!(tsi.rhs_back[tsi.ptr], tsi.psi)
    end
    
    if tsi.steps<0 
        t1 = t + tsi.dt < tsi.tend ? t + tsi.dt : tsi.tend
    else
        t1 = t+1
    end
    return t1, t1
end


function global_orders(psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real, N1::Int; 
                       iter::Int=0, N2::Int=(iter>0?N1+1:0), 
                       rows=8, final_iteration=false, combine_first::Bool=true)
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
        for t in exponential_multistep_equidistant_time_stepper(psi, t0, tend, dt1, N1, N2=N2, steps=steps, 
                         iter=iter, final_iteration=final_iteration, combine_first=combine_first)
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


