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
    tsi, step+1
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
    println("             dt         err      p")
    println("-----------------------------------")
    for row=1:rows
        for t in EquidistantTimeStepper(method, psi, t0, dt1, steps) end
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

    
