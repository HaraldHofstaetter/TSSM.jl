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


function gen_extrapolated_wf!(psi::WaveFunction, L::Vector{Float64}, psi_back::Vector{WaveFunction}, first::Int)
    n = length(psi_back)
    set!(psi, 0)
    for j=1:n
        k = mod(j+first-2, n)+1
        axpy!(psi, psi_back[k], L[j])
    end
end




function gen_starting_values(psi::WaveFunction, dt::Number, N::Int, a::Vector{Float64}, b::Vector{Float64};
         psi_ex::WaveFunction=wave_function(psi.m), rhs::WaveFunction=wave_function(psi.m))
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
                    gen_extrapolated_wf!(psi_ex, L[:,j], psi_back[1:K-1], K)
                    set_time!(psi_ex, get_time(psi))
                    gen_rhs!(rhs, psi_ex)
                    axpy!(psi, rhs, b[j]*dt)
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
   first::Int
   psi_back::Vector{WaveFunction} # storage for previous solution values
   psi_ex::WaveFunction  # storage for the extrapolated values
   rhs::WaveFunction     # storage for the righthandside of the B equation
end


function splitting_with_extrapolation_equidistant_time_stepper(psi::WaveFunction, 
         t0::Real, tend::Real, dt::Real, order::Int, a::Vector{Float64}, b::Vector{Float64}; steps::Int=-1) 
    m = psi.m
    psi_ex = wave_function(m)
    rhs = wave_function(m)
    set_propagate_time_together_with_A!(m, true)
    set_time!(psi, t0)
    psi_back = gen_starting_values(psi, dt, order+1, a, b, psi_ex=psi_ex, rhs=rhs)
    copy!(psi, psi_back[order+1])
    z = cumsum(a)
    L = gen_interpolation_matrix(collect((-order-0.0):0.0), z)
    first = 1
    SplittingWithExtrapolationEquidistantTimeStepperIterator(psi, t0, tend, dt, steps, a, b, L, first, psi_back, psi_ex, rhs)
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
            gen_extrapolated_wf!(tsi.psi_ex, tsi.L[:,j], tsi.psi_back, tsi.first)
            set_time!(tsi.psi_ex, get_time(tsi.psi))
            gen_rhs!(tsi.rhs, tsi.psi_ex)
            axpy!(tsi.psi, tsi.rhs, tsi.b[j]*tsi.dt)
        end
    end            
    copy!(tsi.psi_back[tsi.first], tsi.psi)
    to_real_space!(tsi.psi_back[tsi.first])
    n = length(tsi.psi_back)
    tsi.first = mod(tsi.first, n) + 1
    if tsi.steps<0 
        t1 = t + tsi.dt < tsi.tend ? t + tsi.dt : tsi.tend
    else
        t1 = t+1
    end
    return t1, t1
end


function global_orders(psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real, order::Int, a::Vector{Float64}, b::Vector{Float64}; 
                       operator_sequence="AB", rows=8)
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
        for t in splitting_with_extrapolation_equidistant_time_stepper(psi, t0, tend, dt1, order, a, b, steps=steps)
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
