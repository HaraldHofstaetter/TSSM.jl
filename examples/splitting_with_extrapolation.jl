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
   a::Vector{Float64}
   b::Vector{Float64}
   L::Matrix{Float64}
   first::Int
   psi_back::Vector{WaveFunction} # storage for previous solution values
   psi_ex::WaveFunction  # storage for the extrapolated values
   rhs::WaveFunction     # storage for the righthandside of the B equation
end


function splitting_with_extrapolation_equidistant_time_stepper(psi::WaveFunction, 
         t0::Real, tend::Real, dt::Real, order::Int, a::Vector{Float64}, b::Vector{Float64}) 
    m = psi.m
    psi_ex = wave_function(m)
    rhs = wave_function(m)
    set_propagate_time_together_with_A!(m, true)
    psi_back = gen_starting_values(psi, dt, order+1, a, b, psi_ex=psi_ex, rhs=rhs)
    set_time!(psi, t0)
    copy!(psi, psi_back[order+1])
    z = cumsum(a)
    L = gen_interpolation_matrix(collect((-order-0.0):0.0), [z; 1.0])
    first = 1
    SplittingWithExtrapolationEquidistantTimeStepperIterator(psi, t0, tend, dt, a, b, L, first, psi_back, psi_ex, rhs)
end


Base.start(tsi::SplittingWithExtrapolationEquidistantTimeStepperIterator) = tsi.t0 + (length(tsi.psi_back)-1)*tsi.dt


Base.done(tsi::SplittingWithExtrapolationEquidistantTimeStepperIterator, t) = (t >= tsi.tend)


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
      
    t1 = t + tsi.dt < tsi.tend ? t + tsi.dt : tsi.tend
    t1, t1
end
