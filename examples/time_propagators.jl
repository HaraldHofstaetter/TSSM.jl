using LinearAlgebra
using Printf

abstract type TimePropagationMethod end

COUNT_B = 0

mutable struct EquidistantTimeStepper
    method::TimePropagationMethod
    psi::WaveFunction
    t0::Real
    dt::Real
    steps::Int
    function EquidistantTimeStepper(method::TimePropagationMethod, psi::WaveFunction,
        t0::Real, dt::Real, steps::Int)
        global COUNT_B = 0
        new(method, psi, t0, dt, steps)
    end
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

function Base.iterate(tsi::EquidistantTimeStepper, step=0)
    if step == 0
        set_propagate_time_together_with_A!(tsi.psi.m, true)
        set_time!(tsi.psi, tsi.t0)
        step = initialize!(tsi.method, tsi.psi, tsi.t0, tsi.dt, tsi.steps)
    end
    
    if step >= tsi.steps
        finalize!(tsi.method, tsi.psi, tsi.t0, tsi.dt, tsi.steps, step)
        return nothing
    end
    
    step!(tsi.method, tsi.psi, tsi.t0, tsi.dt, tsi.steps, step)
    return (step+1, tsi), step+1
end

function global_orders(method::TimePropagationMethod, 
                       psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real; rows=8)
    @assert psi.m==reference_solution.m
    tab = zeros(Float64, rows, 3)

    wf_save_initial_value = clone(psi)
    copy!(wf_save_initial_value, psi)

    steps = Int(floor((tend-t0)/dt))
    dt1 = dt
    err_old = 0.0
    println("             dt         err           C      p      B calls")
    println("-----------------------------------------------------------")
    for row=1:rows
        for t in EquidistantTimeStepper(method, psi, t0, dt1, steps) end
        err = distance(psi, reference_solution)
        if (row==1) 
            @Printf.printf("%3i%12.3e%12.3e                   %13i\n", row, Float64(dt1), Float64(err), Int64(COUNT_B))
            tab[row,1] = dt1
            tab[row,2] = err
            tab[row,3] = 0 
        else
            p = log(err_old/err)/log(2.0)
            C = err/dt1^p
            @Printf.printf("%3i%12.3e%12.3e%12.3e%7.2f%13i\n", row, Float64(dt1), Float64(err),
                                                       Float64(C), Float64(p), Int64(COUNT_B))
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

#######################################
#Splitting methods
#######################################

mutable struct SplittingMethod <: TimePropagationMethod
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
            global COUNT_B += 1
        end    
    end         
end


# Splitting with RK4 for B step
mutable struct SplittingRK4BMethod <: TimePropagationMethod
    s::Int
    a::Vector{Float64}
    b::Vector{Float64}
    G #::WaveFunction
    U #::Vector{WaveFunction}
    function SplittingRK4BMethod(a::Vector{Float64}, b::Vector{Float64})         
         s = length(a)
         @assert s==length(b)
         new(s, a, b)
    end     
end

function initialize!(method::SplittingRK4BMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int)
    method.G = wave_function(psi.m)
    method.U = WaveFunction[wave_function(psi.m) for j=1:4]
    0 # iteration  starts with step 0
end

function finalize!(method::SplittingRK4BMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    method.G = nothing
    method.U = nothing
end     

function step!(m::SplittingRK4BMethod, psi::WaveFunction, 
         t0::Real, dt::Real, steps::Int, step::Int)
    for j = 1:m.s
        if m.a[j]!=0.0
            propagate_A!(psi, m.a[j]*dt)
        end
        if m.b[j]!=0.0
            RK4_stepB!(psi, m.b[j]*dt, G=m.G, U=m.U)
        end    
    end         
end

function RK4_stepB!(psi::WaveFunction, dt::Number; 
    # classical explicit Runge-Kutta method for propagating B by dt, using time at t+dt as implied
    # by autonomization for splitting methods.
    G::WaveFunction=wave_function(psi.m),
    U::Vector{WaveFunction}=[wave_function(psi.m) for j=1:4])

    t = get_time(psi)
    set_time!(psi, t+dt)
    for j=1:4
        copy!(U[j], psi)
        set_time!(U[j], t)
    end
    set_time!(psi, t)

    gen_rhs!(G, U[1])
    axpy!(U[2], G, 0.5*dt)
    axpy!(psi, G, 1/6*dt)

    gen_rhs!(G, U[2])
    axpy!(U[3], G, 0.5*dt)
    axpy!(psi, G, 1/3*dt)

    gen_rhs!(G, U[3])
    axpy!(U[4], G, dt)
    axpy!(psi, G, 1/3*dt)

    gen_rhs!(G, U[4])
    axpy!(psi, G, 1/6*dt)
    psi
end

#######################################
#Composition methods
#######################################
function gen_rhs!(rhs::WaveFunction, psi::WaveFunction)
    set!(rhs, 0)
    add_apply_B!(psi,rhs)
    global COUNT_B += 1
end


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

mutable struct CompositionMethod <: TimePropagationMethod
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


#######################################
#(Exponential) Runge-Kutta
#######################################

mutable struct ExponentialRungeKutta <: TimePropagationMethod
    scheme::Int
    s::Int
    iters::Int    
    G #::WaveFunction
    U #::Vector{WaveFunction}
    function ExponentialRungeKutta(scheme::Symbol=:krogstad; iters::Int=(scheme==:gauss_lawson ? 4 : 0))
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
        elseif scheme==:lawson
            scheme1=14
            s = 4
        elseif scheme==:krogstad
            scheme1=42
            s = 4    
    #    elseif scheme==:cox_mathews
    #        scheme1=41
    #        s = 4
        elseif scheme==:strehmel_weiner    
            scheme1=43
            s = 4
        elseif scheme==:gauss_lawson
            scheme1=24
            s = 3     
        else
            warn("scheme not known, I use :krogstad")
        end
        new(scheme1, s, iters, nothing, nothing)
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

function ERK4_step_lawson!(psi::WaveFunction, dt::Number; 
    #Lawson Runge-Kutta method 
    G::WaveFunction=wave_function(psi.m),
    U::Vector{WaveFunction}=[wave_function(psi.m) for j=1:4])
    
    s = 4
    c = [0.0, 0.5, 0.5, 1.0]*dt
    t = get_time(psi)
    copy!(U[1], psi)
    copy!(U[2], U[1])
    propagate_A!(U[2], 0.5*dt)
    copy!(U[3], U[2])
    copy!(U[4], U[3])
    propagate_A!(U[4], 0.5*dt)
    copy!(psi, U[4])
    for j=1:s #seems to be not necessary: propagate_A! already yields right time...
        set_time!(U[j], t+c[j])
    end
    set_time!(psi, t+dt)
    
    gen_rhs!(G, U[1])
    propagate_A!(G, 0.5*dt)
    axpy!(U[2], G, 0.5*dt)
    propagate_A!(G, 0.5*dt)
    axpy!(psi, G, 1/6*dt)  
    
    gen_rhs!(G, U[2])
    axpy!(U[3], G, 0.5*dt)
    propagate_A!(G, 0.5*dt)
    axpy!(psi, G, 1/3*dt)  
    
    gen_rhs!(G, U[3])
    propagate_A!(G, 0.5*dt)
    axpy!(U[4], G, dt)
    axpy!(psi, G, 1/3*dt)      

    gen_rhs!(G, U[4])
    axpy!(psi, G, 1/6*dt)
    psi
end    


function ERK4_step_gauss_lawson!(psi::WaveFunction, dt::Number; 
    iters::Int=4,
    #Lawson Runge-Kutta method 
    G::WaveFunction=wave_function(psi.m),
    U::Vector{WaveFunction}=[wave_function(psi.m) for j=1:3])
    
    c1 = 1/2-sqrt(3)/6
    c2 = 1/2+sqrt(3)/6
    a11 = 1/4
    a12 = 1/4-sqrt(3)/6
    a21 = 1/4+sqrt(3)/6
    a22 = 1/4
    G1 = G
    G2 = U[3]
    copy!(U[1], psi)
    copy!(U[2], psi)
    propagate_A!(U[1], c1*dt)    
    propagate_A!(U[2], c2*dt)        
    
    gen_rhs!(G1, U[1])    
    gen_rhs!(G2, U[2])
 
    for iter = 1:iters #fixed point iteration
        if iter>1
            copy!(U[1], psi)
            copy!(U[2], psi)        
            propagate_A!(U[1], c1*dt)    
            propagate_A!(U[2], c2*dt)                
        end
    
        axpy!(U[1], G1, a11*dt)
        axpy!(U[2], G2, a22*dt)
    
        propagate_A!(G1, (c2-c1)*dt)    
        propagate_A!(G2, (c1-c2)*dt)        

        axpy!(U[1], G2, a12*dt)
        axpy!(U[2], G1, a21*dt)

        gen_rhs!(G1, U[1])    
        gen_rhs!(G2, U[2])
    end
    
    propagate_A!(psi, dt)
    
    propagate_A!(G1, (1.0-c1)*dt)    
    propagate_A!(G2, (1.0-c2)*dt) 
    
    axpy!(psi, G1, 0.5*dt)
    axpy!(psi, G2, 0.5*dt)
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
    elseif m.scheme==14
        ERK4_step_lawson!(psi, dt, G=m.G, U=m.U) 
    elseif m.scheme==24
        ERK4_step_gauss_lawson!(psi, dt, G=m.G, U=m.U, iters=m.iters)         
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

#######################################
# Exponential multistep
#######################################

struct QuadratureRule
    c::Vector{Float64} # nodes normed to interval [0,1]
    b::Vector{Float64} # weights 
end


function gen_interpolation_matrix(x::Vector{Float64}, z::Vector{Float64})
    n = length(x)
    m = length(z)    
    L = zeros(n, m)
    for j=1:m
        for i=1:n
            L[i,j] = prod(Float64[(z[j]-x[k])/(x[i]-x[k]) for k=1:i-1])*prod(Float64[(z[j]-x[k])/(x[i]-x[k]) for k=i+1:n])
        end
    end
    L
end


mutable struct ExponentialMultistep <: TimePropagationMethod
    N::Int
    N1::Int
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
    starting_subdivision::Int # subdivision of dt for calculation of starting values
     
    function ExponentialMultistep(N1::Int; version::Int=1, iters::Int=0,
        final_iteration::Bool=false, combine_first::Bool=true, starting_method::Union{Nothing,TimePropagationMethod}=nothing, 
        quadrature::Union{Nothing,QuadratureRule}=nothing, starting_subdivision::Int=1)
        N = iters==0 ? N1 : N1+1
        @assert version>=0 && version<=2
        a = Float64[]
        b = Float64[]
        if quadrature!=nothing
            version=3
            a = vcat(quadrature.c,1)
            a = vcat(a[1], a[2:end]-a[1:end-1])
            b = vcat(quadrature.b,0)
            C1 = gen_interpolation_matrix(collect((-N1+1.0):0.0), quadrature.c)
            C2 = gen_interpolation_matrix(collect((-N1+1.0):1.0), quadrature.c)
        elseif version==1
            C1 = Matrix{Float64}(inv(Rational{Int}[n^m//factorial(m) for n=-N1+1:0, m=0:N1-1]))
            C2 = Matrix{Float64}(inv(Rational{Int}[n^m//factorial(m) for n=-N1+1:1, m=0:N1]))
        else 
            C1 = Vector{Float64}(Rational{Int}[n^m for m=0:N1-1, n=-N1+1:0]\Rational{Int}[1//(m+1) for m=0:N1-1])
            C2 = Vector{Float64}(Rational{Int}[n^m for m=0:N1,   n=-N1+1:1]\Rational{Int}[1//(m+1) for m=0:N1])            
        end    
        ptr = N        
        new(N, N1, iters, ptr, version, C1, C2, a, b, final_iteration, combine_first, nothing, nothing, nothing, starting_method, starting_subdivision)
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
        for R=1:(final_iteration&&K==N ? 2 : 1)
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
                for m=(combine_first ? 2 : 1):K-1                
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
        for R=1:(final_iteration && K==N ? 2 : 1)
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
        for R=1:(final_iteration && K==N ? 2 : 1)
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
        for I in EquidistantTimeStepper(method.starting_method, psi, t0, dt/method.starting_subdivision, (method.N1-1)*method.starting_subdivision)
            if I[1] % method.starting_subdivision == 0
                k += 1
                gen_rhs!(method.rhs_back[k], psi)
            end
        end
        if method.version==2
            for k=0:method.N1-1
                propagate_A!(method.rhs_back[method.N1-k], dt*k)
            end         
        end            
    else
        if method.version==1
            method.rhs_back = gen_exponential_multistep_starting_values(psi, dt, method.N1, 
                           final_iteration=method.final_iteration, combine_first=method.combine_first)
        elseif method.version==2
            method.rhs_back = gen_exponential_multistep2_starting_values(psi, dt, method.N1, 
                           final_iteration=method.final_iteration)
        elseif method.version==3
            method.rhs_back = gen_exponential_multistep3_starting_values(psi, dt, method.N1, method.a, method.b,
                           final_iteration=method.final_iteration)                   
        end
        if method.iters>=1
            push!(method.rhs_back, wave_function(psi.m))
        end
    end
    if method.iters>=1
        method.psi0 = wave_function(psi.m)
    end    
    if method.version==1 || method.version==3
        method.acc = wave_function(psi.m)
    end    
    method.ptr = method.N1
    method.N1-1 
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
        for n=(m.combine_first ? 2 : 1):m.N1
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
            for n=(m.combine_first ? 2 : 1):m.N1+1
                set!(m.acc, 0.0)
                for k=1:m.N1+1
                    k1 = mod(k-m.N1+m.ptr-2, m.N)+1
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
                    gen_extrapolated_wf!(m.acc, m.C2[:,j], m.rhs_back, m.ptr, order=m.N1+1)
                    axpy!(psi, m.acc, m.b[j]*dt)
                end
            end                                                                   
        elseif m.version==2
            for k=1:m.N1+1
                k1 = mod(k-m.N1+m.ptr-2, m.N)+1
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
    
###############################################################################
# Adaptive Propagators
###############################################################################

abstract type AdaptiveTimePropagationMethod end

mutable struct AdaptiveTimeStepper
    method::AdaptiveTimePropagationMethod

    psi::WaveFunction
    t0::Real
    tend::Real
    tol::Real
    dt::Real
end

#default initializer 
function initialize!(method::AdaptiveTimePropagationMethod, psi::WaveFunction, 
         t0::Real, tend::Real, tol::Real, dt::Real)
    t0 # iteration  starts with step 0
end

#default finalizer
function finalize!(method::AdaptiveTimePropagationMethod, psi::WaveFunction, 
         t0::Real, tend::Real, tol::Real, dt::Real, t::Real)
end     

function Base.start(tsi::AdaptiveTimeStepper) 
    set_propagate_time_together_with_A!(tsi.psi.m, true)
    set_time!(tsi.psi, tsi.t0)
    initialize!(tsi.method, tsi.psi, tsi.t0, tsi.tend, tsi.tol, tsi.dt)
    tsi.t0
end

function Base.done(tsi::AdaptiveTimeStepper, t::Real) 
    f = (t>=tsi.tend )
    if f 
        finalize!(tsi.method, tsi.psi, tsi.t0, tsi.tend, tsi.tol, tsi.dt, t)
    end
    f    
end
    
function Base.next(tsi::AdaptiveTimeStepper, t::Real)    
    (tnew, dtnew) = step!(tsi.method, tsi.psi, tsi.t0, tsi.tend, tsi.tol, tsi.dt, t)
    tsi.dt = dtnew
    (tnew, tsi), tnew 
end

########################################################################################

function solve_vander_trans(x::Vector, b::Vector)
    # Algorithm 4.6.2 from Golub/van Loan
    n = length(x)
    for k=1:n-1
        for i=n:-1:k+1
            b[i] -= x[k]*b[i-1]
        end
    end
    for k=n-1:-1:1
        for i=k+1:n
            b[i] /= (x[i]-x[i-k])
        end
        for i=k:n-1
            b[i] -= b[i+1]
        end
    end
end

#c = Float64[dt^m/(m+1) for m=0:N1-1]
#solve_vander_trans(t_back, c)
#c


mutable struct AdaptiveAdamsLawson <: AdaptiveTimePropagationMethod
    N::Int
    N1::Int
    ptr::Int
    t_back::Vector{Float64}
    rhs_back #::Vector{WaveFunction} # storage for previous solution values
    psi0 #::WaveFunction
    psi1 #::WaveFunction
    acc  #::WaveFunction
    starting_method #::TimePropagationMethod
    bootstrap_mode::Bool
    N1_final::Int

    combine_first::Bool
    version::Int
    
     
    function AdaptiveAdamsLawson(N1::Int; starting_method::Union{Nothing,TimePropagationMethod}=nothing,
                                combine_first::Bool=true, version::Int=2)
        new(N1+1, N1, N1, Float64[], nothing, nothing, nothing, nothing, starting_method, starting_method==nothing, N1, combine_first, version)
    end
end    

function initialize!(m::AdaptiveAdamsLawson, psi::WaveFunction, 
         t0::Real, tend::Real, tol::Real, dt::Real)    
    m.bootstrap_mode = m.starting_method==nothing
    N_final = m.N1_final+1
    m.t_back = zeros(Float64, N_final)
    m.rhs_back = WaveFunction[wave_function(psi.m) for j=1:N_final]              
    m.psi0 = wave_function(psi.m)
    m.psi1 = wave_function(psi.m)
    if m.version==1 
        m.acc = wave_function(psi.m)
    end
    gen_rhs!(m.rhs_back[1], psi)
    if m.bootstrap_mode
        m.N = 2
        m.N1 = 1
    else    
        k=1   
        for I in EquidistantTimeStepper(m.starting_method, psi, t0, dt, m.N1-1) 
            k += 1
            gen_rhs!(m.rhs_back[k], psi)
        end
        m.t_back[:] = Float64[dt*k for k=0:m.N1]
        if m.version==2 #Lawson
            for k=0:m.N1-1
                propagate_A!(m.rhs_back[m.N1-k], dt*k)
            end         
        end
    end
    m.ptr = m.N1
    (m.N1-1)*dt 
end    

function finalize!(method::AdaptiveAdamsLawson, psi::WaveFunction, 
         t0::Real, tend::Real, tol::Real, dt::Real, t::Real)
    method.psi0 = nothing
    method.psi1 = nothing
    method.acc = nothing
    method.rhs_back = nothing
end     

function step!(m::AdaptiveAdamsLawson, psi::WaveFunction, 
         t0::Real, tend::Real, tol::Real, dt::Real, t::Real)
    facmin = 0.25
    facmax = 4.0
    fac = 0.9
    
    copy!(m.psi0, psi)
    copy!(m.psi1, psi)

    ptr0 = m.ptr
    dt0 = dt
    
    err = 2.0 #error/t    
    while err>=1.0
        dt = min(dt, tend-t)
        dt0 = dt
        
        #predictor 
        tt = Float64[m.t_back[mod(k-m.N1+m.ptr-1, m.N)+1] for k=1:m.N1]
        tt = tt .- tt[end]
        if m.version==2 # Lawson
            C = Float64[dt^m/(m+1) for m=0:m.N1-1]
            solve_vander_trans(tt, C)
            for k=1:m.N1
                k1 = mod(k-m.N1+m.ptr-1, m.N)+1
                axpy!(m.psi1, m.rhs_back[k1], C[k]*dt)
            end
            propagate_A!(m.psi1, dt)    
        elseif m.version==1
            C = inv([(t/dt)^m/factorial(m) for t in tt, m=0:m.N1-1])
            if m.combine_first
                k1 = mod(m.ptr-1, m.N)+1
                copy!(m.acc, m.rhs_back[k1])
                add_apply_A!(m.psi1, m.acc, 1.0)
                add_phi_A!(m.acc, m.psi1, dt, 1, dt)
                set_time!(m.psi1, get_time(m.psi1)+dt)
            else
                propagate_A!(m.psi1, dt)
            end
            for n=(m.combine_first ? 2 : 1):m.N1
                set!(m.acc, 0.0)
                for k=1:m.N1
                    k1 = mod(k-m.N1+m.ptr-1, m.N)+1
                    axpy!(m.acc, m.rhs_back[k1], C[n,k])
                end
                add_phi_A!(m.acc, m.psi1, dt, n, dt)
            end
        end
        
        if m.bootstrap_mode
            m.ptr += 1
        else
            m.ptr = mod(m.ptr, m.N) + 1
        end
        gen_rhs!(m.rhs_back[m.ptr], m.psi1)    
        m.t_back[m.ptr] = t+dt
        if m.version==2 
            propagate_A!(m.rhs_back[m.ptr], -dt)
        end
        
        #corrector
        tt = Float64[m.t_back[mod(k-m.N1+m.ptr-2, m.N)+1] for k=1:m.N1+1]
        tt = tt .- tt[end-1]
        if m.version==2 # Lawson
            C = Float64[dt^m/(m+1) for m=0:m.N1]       
            solve_vander_trans(tt, C)
            for k=1:m.N1+1
                k1 = mod(k-m.N1+m.ptr-2, m.N)+1
                axpy!(psi, m.rhs_back[k1], C[k]*dt)
            end
            propagate_A!(psi, dt)                        
        elseif m.version==1
            C = inv([(t/dt)^m/factorial(m) for t in tt, m=0:m.N1])
            if m.combine_first
                k1 = mod(m.ptr-2, m.N)+1
                copy!(m.acc, m.rhs_back[k1])
                add_apply_A!(psi, m.acc, 1.0)
                add_phi_A!(m.acc, psi, dt, 1, dt)
                set_time!(psi, get_time(psi)+dt)
            else
                propagate_A!(psi, dt)
            end        
            for n=(m.combine_first ? 2 : 1):m.N1+1
                set!(m.acc, 0.0)
                for k=1:m.N1+1
                    k1 = mod(k-m.N1+m.ptr-2, m.N)+1
                    axpy!(m.acc, m.rhs_back[k1], C[n,k])
                end
                add_phi_A!(m.acc, psi, dt, n, dt)
            end  

        end
        
        err = distance(psi, m.psi1)/tol
        dt = dt*min(facmax, max(facmin, fac*(1.0/err)^(1.0/(m.N1+1)))) # CHECK m.N+1=p+1 !!!
        
        if err>=1.0
            copy!(m.psi1, m.psi0)
            copy!(psi, m.psi0)
            m.ptr = ptr0
            @Printf.printf("t=%17.9e  err=%17.8e  dt=%17.8e  rejected...\n", 
                    Float64(t), Float64(err), Float64(dt))
        end
    end    
    
    gen_rhs!(m.rhs_back[m.ptr], psi)    
    t += dt0
    m.t_back[m.ptr] = t
    
    if m.version==2
        for k=1:m.N
            if k==m.ptr continue end
            propagate_A!(m.rhs_back[k], dt0)
        end   
    end
    if m.bootstrap_mode 
        if m.N1_final>m.N1
            m.N1 += 1
            m.N += 1
        else
            m.bootstrap_mode = false
        end
    end    
    (t, dt)
end    

