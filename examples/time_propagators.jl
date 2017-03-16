abstract TimePropagationMethod

type EquidistantTimeStepper
    method::TimePropagationMethod

    psi::WaveFunction
    t0::Real
    dt::Real
    steps::Int
end

#default initializer 
function initialize!(method::TimePropagationMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int)
    0 # iteration  starts with step 0
end

#default finalizer
function finalize!(method::TimePropagationMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
end     

function Base.start(tsi::EquidistantTimeStepper) 
    set_propagate_time_together_with_A!(tsi.psi.m, true)
    set_time!(tsi.psi, tsi.t0)
    initialize!(tsi.method, tsi.psi, tsi.t0, tsi.dt, tsi.steps)
end

function Base.done(tsi::EquidistantTimeStepper, step::Int) 
    f = (step>=tsi.steps )
    if f 
        finalize!(tsi.method, tsi.psi, tsi.t0, tsi.dt, tsi.steps, step)
    end
    f    
end
    
function Base.next(tsi::EquidistantTimeStepper, step::Int)    
    step!(tsi.method, tsi.psi, tsi.t0, tsi.dt, tsi.steps, step)
    (step+1, tsi), step+1
end

function global_orders(method::TimePropagationMethod, 
                       psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real; rows=8)
    @assert psi.m==reference_solution.m
    tab = Array(Float64, rows, 3)

    wf_save_initial_value = clone(psi)
    copy!(wf_save_initial_value, psi)

    steps = Int(floor((tend-t0)/dt))
    dt1 = dt
    err_old = 0.0
    println("             dt         err           C      p")
    println("-----------------------------------------------")
    for row=1:rows
        for t in EquidistantTimeStepper(method, psi, t0, dt1, steps) end
        err = distance(psi, reference_solution)
        if (row==1) then
            @printf("%3i%12.3e%12.3e\n", row, Float64(dt1), Float64(err))
            tab[row,1] = dt1
            tab[row,2] = err
            tab[row,3] = 0 
        else
            p = log(err_old/err)/log(2.0)
            C = err/dt1^p
            @printf("%3i%12.3e%12.3e%12.3e%7.2f\n", row, Float64(dt1), Float64(err), Float64(C), Float64(p))
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

#######################################################################################
#Splitting methods


type SplittingMethod <: TimePropagationMethod
    s::Int
    a::Vector{Float64}
    b::Vector{Float64}
    function SplittingMethod(a::Vector{Float64}, b::Vector{Float64})         
         s = length(a)
         @assert s==length(b)
         new(s, a, b)
    end     
end

function step!(m::SplittingMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    for j = 1:m.s
        if m.a[j]!=0.0
            propagate_A!(psi, m.a[j]*dt)            
        end
        if m.b[j]!=0.0
            propagate_B!(psi, m.b[j]*dt) 
        end    
    end         
end  




#######################################################################################
#Composition methods


function get_coeffs_composition(g::Vector{Float64})
    n = length(g)
    a = zeros(n+1)
    b = vcat(g, 0)
    a[1] = g[1]/2
    a[n+1] = a[1]
    for j=2:div(n,2)+1
        a[j] = (g[j-1]+g[j])/2
        a[n-j+2] = a[j]
    end
    a,b
end

type CompositionMethod <: TimePropagationMethod
    s::Int
    c::Vector{Float64}    
    a::Vector{Float64}
    b::Vector{Float64}
    iters::Int
    psi0 #::WaveFunction
    k1 #::WaveFunction
    k2 #::WaveFunction
    function CompositionMethod(c::Vector{Float64}, iters::Int)
         s = length(c)+1
         a,b = get_coeffs_composition(c)
         new(s, c, a, b, iters, nothing, nothing, nothing)
    end     
end


function initialize!(method::CompositionMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int)
    method.k1 = wave_function(psi.m)
    method.k2 = wave_function(psi.m)
    method.psi0 = wave_function(psi.m)    
    0 # iteration  starts with step 0
end


function finalize!(method::CompositionMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    method.k1 = nothing
    method.k2 = nothing
    method.psi0 = nothing
end     


function step!(m::CompositionMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    for j = 1:m.s
        if m.a[j]!=0.0
            propagate_A!(psi, m.a[j]*dt)            
        end
        if m.b[j]!=0.0
            if m.iters>0
                copy!(m.psi0, psi)
            end            
            
            #Euler step:
            gen_rhs!(m.k1, psi)
            axpy!(psi, m.k1, m.b[j]*dt)            
            
            #fixed point iteration for implicit midpoint rule
            for it=1:m.iters
                copy!(m.k1, m.psi0)
                scale!(m.k1, 0.5)
                axpy!(m.k1, psi, 0.5)
                copy!(psi, m.psi0)
                set_time!(m.k1, get_time(psi))
                gen_rhs!(m.k2, m.k1)
                axpy!(psi, m.k2, m.b[j]*dt)                    
            end
        end 
    end         
end  


#######################################################################################
#(Exponential) Runge-Kutta


type ExponentialRungeKutta <: TimePropagationMethod
    scheme::Int
    s::Int
    G #::WaveFunction
    U #::Vector{WaveFunction}
    function ExponentialRungeKutta(scheme::Symbol=:krogstad)
        scheme1 = 42
        s = 4
        if scheme==:rk4
            scheme1=4
            s = 4
        elseif scheme==:erk2
            scheme1=39
            s = 2   
        elseif scheme==:erk2a
            scheme1=40
            s = 2           
        elseif scheme==:krogstad
            scheme1=42
            s = 4    
    #    elseif scheme==:cox_mathews
    #        scheme1=41
    #        s = 4
        elseif scheme==:strehmel_weiner    
            scheme1=43
            s = 4
        else
            warn("scheme not known, I use :krogstad")
        end
        new(scheme1, s, nothing, nothing)
    end         
end


function initialize!(method::ExponentialRungeKutta, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int)
    method.G = wave_function(psi.m)
    method.U = WaveFunction[wave_function(psi.m) for j=1:method.s]
    0 # iteration  starts with step 0
end

function finalize!(method::ExponentialRungeKutta, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    method.G = nothing
    method.U = nothing
end     




function ERK2_step!(psi::WaveFunction, dt::Number; 
    #eq. (2.39) in Hochbruck/Ostermann, c2=1
    G::WaveFunction=wave_function(psi.m),
    U::Vector{WaveFunction}=[wave_function(psi.m) for j=1:2])

    s = 2
    c = [0.0, 1.0]*dt
    t = get_time(psi)
    for j=1:s
        copy!(U[j], psi)
        set_time!(U[j], t+c[j])
    end
    set_time!(psi, t+dt)

    gen_rhs!(G, U[1])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, U[2], c[2], 1, dt)
    add_phi_A!(G, psi, dt, 1, dt)
    add_phi_A!(G, psi, dt, 2, -dt)  
    
    gen_rhs!(G, U[2])
    add_apply_A!(U[1], G, 1.0)   
    add_phi_A!(G, psi, dt, 2, dt)    
end    
    
function ERK2A_step!(psi::WaveFunction, dt::Number; 
    #eq. (2.40) in Hochbruck/Ostermann, c2=0.5
    G::WaveFunction=wave_function(psi.m),
    U::Vector{WaveFunction}=[wave_function(psi.m) for j=1:2])

    s = 2
    c = [0.0, 0.5]*dt
    t = get_time(psi)
    for j=1:s
        copy!(U[j], psi)
        set_time!(U[j], t+c[j])
    end
    set_time!(psi, t+dt)

    gen_rhs!(G, U[1])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, U[2], c[2], 1, 0.5*dt)
        
    gen_rhs!(G, U[2])
    add_apply_A!(U[1], G, 1.0)   
    add_phi_A!(G, psi, dt, 1, dt)    
end    


function ERK4_step_krogstad!(psi::WaveFunction, dt::Number; 
    #Krogstad, eq. (2.42) in Hochbruck/Ostermann
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

function RK4_step!(psi::WaveFunction, dt::Number; 
    #classical explicit Runge-Kutta method (NOT exponential!)
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
    axpy!(U[2], G, 0.5*dt)
    axpy!(psi, G, 1/6*dt)
    
    gen_rhs!(G, U[2])
    add_apply_A!(U[2], G, 1.0)
    axpy!(U[3], G, 0.5*dt)
    axpy!(psi, G, 1/3*dt)

    gen_rhs!(G, U[3])
    add_apply_A!(U[3], G, 1.0)
    axpy!(U[4], G, dt)
    axpy!(psi, G, 1/3*dt)
    
    gen_rhs!(G, U[4])
    add_apply_A!(U[4], G, 1.0)
    axpy!(psi, G, 1/6*dt)
    psi
end    



function ERK4_step_strehmel_weiner!(psi::WaveFunction, dt::Number; 
    #Strehmel/Weiner, eq. (2.43) in Hochbruck/Ostermann
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
    add_phi_A!(G, U[3], c[3], 2, -0.5*dt)
    add_phi_A!(G, U[4], c[4], 1, 1.0*dt)
    add_phi_A!(G, U[4], c[4], 2, -2.0*dt)
    add_phi_A!(G, psi, dt, 1, 1.0*dt)
    add_phi_A!(G, psi, dt, 2, -3.0*dt)
    add_phi_A!(G, psi, dt, 3, 4.0*dt)
    
    gen_rhs!(G, U[2])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, U[3], c[3], 2, 0.5*dt)
    add_phi_A!(G, U[4], c[4], 2, -2.0*dt)
    
    gen_rhs!(G, U[3])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, U[4], c[4], 2, 4.0*dt)
    add_phi_A!(G, psi, dt, 2, 4.0*dt)
    add_phi_A!(G, psi, dt, 3, -8.0*dt)

    gen_rhs!(G, U[4])
    add_apply_A!(U[1], G, 1.0)
    add_phi_A!(G, psi, dt, 2, -1.0*dt)
    add_phi_A!(G, psi, dt, 3, 4.0*dt)
    psi
end

function step!(m::ExponentialRungeKutta, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    if m.scheme==4
        RK4_step!(psi, dt, G=m.G, U=m.U) 
    elseif m.scheme==39
        ERK2_step!(psi, dt, G=m.G, U=m.U)    
    elseif m.scheme==40
        ERK2A_step!(psi, dt, G=m.G, U=m.U)    
    elseif m.scheme==42
        ERK4_step_krogstad!(psi, dt, G=m.G, U=m.U)    
    elseif m.scheme==43
        ERK4_step_strehmel_weiner!(psi, dt, G=m.G, U=m.U)    
    elseif m.scheme==41
        ERK4_step_cox_mathews!(psi, dt, G=m.G, U=m.U)    
    end
 end

########################################################################################
# Exponential multistep

immutable QuadratureRule
    c::Vector{Float64} # nodes normed to interval [0,1]
    b::Vector{Float64} # weights 
end


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


type ExponentialMultistep <: TimePropagationMethod
    N::Int
    N1::Int
    N2::Int
    iters::Int    
    ptr::Int
    version::Int
    C1::Union{Matrix{Float64}, Vector{Float64}}
    C2::Union{Matrix{Float64}, Vector{Float64}}
    a::Vector{Float64}
    b::Vector{Float64}
    final_iteration::Bool
    combine_first::Bool    
    rhs_back #::Vector{WaveFunction} # storage for previous solution values
    psi0 #::WaveFunction
    acc #::WaveFunction    
    starting_method #::TimePropagationMethod
     
    function ExponentialMultistep(N1::Int; version::Int=1, iters::Int=0, N2::Int=(iters>0?N1+1:0),
        final_iteration::Bool=false, combine_first::Bool=true, starting_method::Union{Void,TimePropagationMethod}=nothing, 
        quadrature::Union{Void,QuadratureRule}=nothing)
        N = max(N1, N2)
        @assert version>=0 && version<=2
        a = Float64[]
        b = Float64[]
        if quadrature!=nothing
            version=3
            a = vcat(quadrature.c,1)
            a = vcat(a[1], a[2:end]-a[1:end-1])
            b = vcat(quadrature.b,0)
            C1 = gen_interpolation_matrix(collect((-N1+1.0):0.0), quadrature.c)
            C2 = gen_interpolation_matrix(collect((-N2+2.0):1.0), quadrature.c)
        elseif version==1
            C1 = Matrix{Float64}(inv(Rational{Int}[n^m//factorial(m) for n=-N1+1:0, m=0:N1-1]))
            C2 = Matrix{Float64}(inv(Rational{Int}[n^m//factorial(m) for n=-N2+2:1, m=0:N2-1]))
        else 
            C1 = Vector{Float64}(Rational{Int}[n^m for m=0:N1-1, n=-N1+1:0]\Rational{Int}[1//(m+1) for m=0:N1-1])
            C2 = Vector{Float64}(Rational{Int}[n^m for m=0:N2-1, n=-N2+2:1]\Rational{Int}[1//(m+1) for m=0:N2-1])            
        end    
        ptr = N        
        new(N, N1, N2, iters, ptr, version, C1, C2, a, b, final_iteration, combine_first, nothing, nothing, nothing, starting_method)
    end
end    

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
                C  = Matrix{Float64}(inv(Rational{Int}[k^m//factorial(m) for k=-J+2:K-J, m=0:K-2])) 
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
                        axpy!(acc, rhs_back[k], C[m,k])
                    end
                    add_phi_A!(acc, psi, dt, m, dt)        
                end 
                gen_rhs!(rhs_back[J], psi)   
            end
        end
    end
    rhs_back
end

function gen_exponential_multistep2_starting_values(psi::WaveFunction, dt::Number, N::Int; 
          final_iteration=false)    
    psi0 = wave_function(psi.m)
    rhs_back = WaveFunction[wave_function(psi.m) for j=1:N]
    gen_rhs!(rhs_back[1], psi)
    t0 = get_time(psi)
    copy!(psi0, psi)
    for K=2:N
        for R=1:(final_iteration&&K==N?2:1)
            copy!(psi, psi0)            
            for k=1:K-1
                propagate_A!(rhs_back[k], -(K-2)*dt) 
            end 
            for J=2:K
                C = Vector{Float64}(Rational{Int}[n^m for m=0:K-2, n=-J+2:K-J]\Rational{Int}[1//(m+1) for m=0:K-2])
                for k=1:K-1
                   axpy!(psi, rhs_back[k], C[k]*dt)
                end
                propagate_A!(psi, dt)
                gen_rhs!(rhs_back[J], psi)   
                propagate_A!(rhs_back[J], -dt) 
                for k=1:K
                    propagate_A!(rhs_back[k], dt) 
                end
            end            
        end
    end
    rhs_back
end

function gen_extrapolated_wf!(psi::WaveFunction, L::Vector{Float64}, psi_back::Vector{WaveFunction}, ptr::Int; 
                              order::Int=length(psi_back))
    n = length(psi_back)
    set!(psi, 0)
    for j=1:order
        k = mod(j-order+ptr-1, n)+1
        axpy!(psi, psi_back[k], L[j])
    end
end

function gen_exponential_multistep3_starting_values(psi::WaveFunction, dt::Number, N::Int, 
    a::Vector{Float64}, b::Vector{Float64};final_iteration=false) 
    rhs = wave_function(psi.m)
    s = length(a)
    aa = cumsum(a)
    psi0 = wave_function(psi.m)
    rhs_back = WaveFunction[wave_function(psi.m) for j=1:N]
    copy!(psi0, psi)   
    gen_rhs!(rhs_back[1], psi)
    t0 = get_time(psi)
    for K=2:N
        for R=1:(final_iteration&&K==N?2:1)
            copy!(psi, psi0)
            for J=2:K
                L = gen_interpolation_matrix(collect(0.0:(K-2.0)), J-2.0+aa) # TODO: check order->order+/-1            
                for j = 1:s
                    if a[j]!=0.0
                        propagate_A!(psi, a[j]*dt)            
                    end
                    if b[j]!=0.0
                        gen_extrapolated_wf!(rhs, L[:,j], rhs_back[1:K-1], 2*K-2) # TODO: check order->order+/-1 
                        axpy!(psi, rhs, b[j]*dt)
                    end
                end            
                gen_rhs!(rhs_back[J], psi)
            end
        end
    end 
    rhs_back
end

function initialize!(method::ExponentialMultistep, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int)
    if method.starting_method!=nothing
        method.rhs_back = WaveFunction[wave_function(psi.m) for j=1:method.N]              
        gen_rhs!(method.rhs_back[1], psi)
        k=1  
        for I in EquidistantTimeStepper(method.starting_method, psi, t0, dt, method.N-1) 
            k += 1
            gen_rhs!(method.rhs_back[k], psi)
        end
        if method.version==2
            for k=0:method.N-1
                propagate_A!(method.rhs_back[method.N-k], dt*k)
            end         
        end            
    else
        if method.version==1
            method.rhs_back = gen_exponential_multistep_starting_values(psi, dt, method.N, 
                           final_iteration=method.final_iteration, combine_first=method.combine_first)
        elseif method.version==2
            method.rhs_back = gen_exponential_multistep2_starting_values(psi, dt, method.N, 
                           final_iteration=method.final_iteration)
        elseif method.version==3
            method.rhs_back = gen_exponential_multistep3_starting_values(psi, dt, method.N, method.a, method.b,
                           final_iteration=method.final_iteration)                   
        end
    end
    if method.iters>=1
        method.psi0 = wave_function(psi.m)
    end    
    if method.version==1 || method.version==3
        method.acc = wave_function(psi.m)
    end    
    method.ptr = method.N
    method.N-1 
end


function finalize!(method::ExponentialMultistep, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    method.psi0 = nothing
    method.acc = nothing
    method.rhs_back = nothing
end     


function step!(m::ExponentialMultistep, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)

    if m.iters>=1
        copy!(m.psi0, psi)
    end
    
    #predictor
    if m.version==1
        if m.combine_first
            k1 = mod(m.ptr-1, m.N)+1
            copy!(m.acc, m.rhs_back[k1])
            add_apply_A!(psi, m.acc, 1.0)
            add_phi_A!(m.acc, psi, dt, 1, dt)
            set_time!(psi, get_time(psi)+dt)
        else
            propagate_A!(psi, dt)
        end
        for n=(m.combine_first?2:1):m.N1
            set!(m.acc, 0.0)
            for k=1:m.N1
                k1 = mod(k-m.N1+m.ptr-1, m.N)+1
                axpy!(m.acc, m.rhs_back[k1], m.C1[n,k])
            end
            add_phi_A!(m.acc, psi, dt, n, dt)
        end
    elseif m.version==3            
        for j = 1:length(m.a)
            if m.a[j]!=0.0
                propagate_A!(psi, m.a[j]*dt)            
            end
            if m.b[j]!=0.0
                gen_extrapolated_wf!(m.acc, m.C1[:,j], m.rhs_back, m.ptr, order=m.N1)
                axpy!(psi, m.acc, m.b[j]*dt)
            end    
        end                            
    elseif m.version==2
        for k=1:m.N1
            k1 = mod(k-m.N1+m.ptr-1, m.N)+1
            axpy!(psi, m.rhs_back[k1], m.C1[k]*dt)
        end
        propagate_A!(psi, dt)
    end
   
    m.ptr = mod(m.ptr, m.N) + 1
    gen_rhs!(m.rhs_back[m.ptr], psi)    
    if m.version==2 && m.iters>=1
        propagate_A!(m.rhs_back[m.ptr], -dt)
    end
        
    for iter=1:m.iters
        copy!(psi, m.psi0)
        
        #corrector
        if m.version==1
            if m.combine_first
                k1 = mod(m.ptr-2, m.N)+1
                copy!(m.acc, m.rhs_back[k1])
                add_apply_A!(psi, m.acc, 1.0)
                add_phi_A!(m.acc, psi, dt, 1, dt)
                set_time!(psi, get_time(psi)+dt)
            else
                propagate_A!(psi, dt)
            end        
            for n=(m.combine_first?2:1):m.N2
                set!(m.acc, 0.0)
                for k=1:m.N2
                    k1 = mod(k-m.N2+m.ptr-1, m.N)+1
                    axpy!(m.acc, m.rhs_back[k1], m.C2[n,k])
                end
                add_phi_A!(m.acc, psi, dt, n, dt)
            end  
        elseif m.version==3
            for j = 1:length(m.a)
                if m.a[j]!=0.0
                    propagate_A!(psi, m.a[j]*dt)            
                end
                if m.b[j]!=0.0
                    gen_extrapolated_wf!(m.acc, m.C2[:,j], m.rhs_back, m.ptr, order=m.N2)
                    axpy!(psi, m.acc, m.b[j]*dt)
                end
            end                                                                   
        elseif m.version==2
            for k=1:m.N2
                k1 = mod(k-m.N2+m.ptr-1, m.N)+1
                axpy!(psi, m.rhs_back[k1], m.C2[k]*dt)
            end
            propagate_A!(psi, dt)                        
        end
        gen_rhs!(m.rhs_back[m.ptr], psi)
    end
    
    if m.version==2
        for k=1:m.N
            if k==m.ptr continue end
            propagate_A!(m.rhs_back[k], dt)
        end    
    end
end
    
    
