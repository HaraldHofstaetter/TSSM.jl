
Strang = (0.5, 1.0, 0.5)

function step!(psi::WaveFunction, dt::Real, scheme=Strang, operator_sequence="AB")
    for k = 1:length(scheme)
        which_operator = operator_sequence[mod(k-1, length(operator_sequence)) + 1]
        if which_operator == 'A'
            propagate_A!(psi, dt*scheme[k]) 
        elseif which_operator == 'B'
            propagate_B!(psi, dt*scheme[k]) 
        elseif which_operator == 'C'
            propagate_C!(psi, dt*scheme[k]) 
        elseif which_operator == 'T'
            propagate_time!(psi, dt*scheme[k]) 
        else
            # TODO! Error or Warning
        end
    end 
end

struct EquidistantTimeStepperIterator
   psi::WaveFunction
   t0::Real
   tend::Real
   dt::Real
   scheme
   operator_sequence
end

equidistant_time_stepper(psi::WaveFunction, t0::Real, tend::Real, dt::Real, 
                         scheme::Tuple=Strang, operator_sequence="AB") =
	EquidistantTimeStepperIterator(psi, t0, tend, dt, scheme, operator_sequence)

Base.start(tsi::EquidistantTimeStepperIterator) = tsi.t0

Base.done(tsi::EquidistantTimeStepperIterator, t) = t >= tsi.tend

function Base.next(tsi::EquidistantTimeStepperIterator, t)
    step!(tsi.psi, tsi.dt, tsi.scheme, tsi.operator_sequence)
    t1 = t + tsi.dt < tsi.tend ? t + tsi.dt : tsi.tend
    t1, t1
end

function step_embedded!(psi1::WaveFunction, psi2::WaveFunction, dt::Real, scheme1, scheme2, operator_sequence="AB")
    kk = -1
    for k = 1:length(scheme1)
        if kk<0 && scheme1[k] != scheme2[k]
             kk = k
             copy!(psi2, psi1)
        end     
        which_operator = operator_sequence[mod(k-1, length(operator_sequence)) + 1]
        if which_operator == 'A'
            propagate_A!(psi1, dt*scheme1[k]) 
        elseif which_operator == 'B'
            propagate_B!(psi1, dt*scheme1[k]) 
        elseif which_operator == 'C'
            propagate_C!(psi1, dt*scheme1[k]) 
        elseif which_operator == 'T'
            propagate_time!(psi1, dt*scheme[k]) 
        else
            # TODO! Error or Warning
        end
    end 
    if kk<0
       kk = length(scheme1)+1
       copy!(psi2, psi1)
    end   
    for k = kk:length(scheme2)
        which_operator = operator_sequence[mod(k-1, length(operator_sequence)) + 1]
        if which_operator == 'A'
            propagate_A!(psi2, dt*scheme2[k]) 
        elseif which_operator == 'B'
            propagate_B!(psi2, dt*scheme2[k]) 
        elseif which_operator == 'C'
            propagate_C!(psi2, dt*scheme2[k]) 
        elseif which_operator == 'T'
            propagate_time!(psi2, dt*scheme[k]) 
        else
            # TODO! Error or Warning
        end
    end 
end

function step_palindromic!(psi1::WaveFunction, psi2::WaveFunction, dt::Real, scheme, operator_sequence="AB")
    copy!(psi2, psi1)
    step!(psi1, dt, scheme, operator_sequence)
    step!(psi2, dt, scheme, reverse(operator_sequence))
end

function step_defect_based!(psi::WaveFunction, h::WaveFunction, dt::Real, scheme, 
                            operator_sequence="AB"; symmetrized_defect::Bool=false )
    set!(h, 0.0) #TODO: set is_real_space=false for AB, true for BA    
    set_time!(h, 0.0)
    s = length(scheme)
    m =  length(operator_sequence)

    f = 1.0
    if symmetrized_defect
        f = 0.5
        which_operator = operator_sequence[1]
        if 'A' in operator_sequence && which_operator != 'A' 
            add_apply_A!(psi, h, -f)
        end
        if 'B' in operator_sequence && which_operator != 'B' 
            add_apply_B!(psi, h, -f)
        end
        if 'C' in operator_sequence && which_operator != 'C' 
            add_apply_C!(psi, h, -f)
        end
    end
    which_operator = 'X'
    for k = 1:s
        which_operator = operator_sequence[mod(k-1, m) + 1]

        y = scheme[k]
        if (symmetrized_defect && k==1) || (k==s)
            y-=f
        end

        if which_operator == 'A'
            add_apply_A!(psi, h, y)
        elseif which_operator == 'B'
            add_apply_B!(psi, h, y)
        elseif which_operator == 'C'
            add_apply_C!(psi, h, y)
        end

        if which_operator == 'A'
            propagate_A_derivative!(psi, h, dt*scheme[k])
        elseif which_operator == 'B'
            propagate_B_derivative!(psi, h, dt*scheme[k])
        elseif which_operator == 'C'
            propagate_C_derivative!(psi, h, dt*scheme[k])
        end    

    end
    if 'A' in operator_sequence && which_operator != 'A' 
        add_apply_A!(psi, h, -f)
    end
    if 'B' in operator_sequence && which_operator != 'B' 
        add_apply_B!(psi, h, -f)
    end
    if 'C' in operator_sequence && which_operator != 'C' 
        add_apply_C!(psi, h, -f)
    end
end
 

struct AdaptiveTimeStepperIterator
   psi::WaveFunction
   t0::Real
   tend::Real
   dt::Real
   tol::Real
   scheme1
   scheme2
   order::Integer
   operator_sequence
   psi2::WaveFunction
   psi0::WaveFunction
end

struct AdaptiveTimeStepperState
   t::Real
   dt::Real
end   

function adaptive_time_stepper(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, 
                         scheme1, scheme2, order::Integer, operator_sequence="AB")
    psi2=clone(psi)                     
    psi0=clone(psi)                     
    AdaptiveTimeStepperIterator(psi, t0, tend, dt, tol,
        scheme1, scheme2, order, operator_sequence, psi2, psi0)
end        

Base.start(tsi::AdaptiveTimeStepperIterator) = AdaptiveTimeStepperState(tsi.t0, tsi.dt)

function Base.done(tsi::AdaptiveTimeStepperIterator, state::AdaptiveTimeStepperState)
  state.t >= tsi.tend
end  

using Printf

function Base.next(tsi::AdaptiveTimeStepperIterator, state::AdaptiveTimeStepperState)
    facmin = 0.25
    facmax = 4.0
    fac = 0.9

    dt = state.dt
    dt0 = dt
    copy!(tsi.psi0, tsi.psi)
    err = 2.0
    while err>=1.0
        dt = min(dt, tsi.tend-state.t)
        dt0 = dt
        if tsi.scheme2=="palindromic"
           step_palindromic!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.operator_sequence)
           err = 0.5*distance(tsi.psi, tsi.psi2)/tsi.tol
        elseif tsi.scheme2=="defect_based"
           step_defect_based!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.operator_sequence)
           err = dt*norm(tsi.psi2)/(tsi.order+1)/tsi.tol
        elseif tsi.scheme2=="symmetrized_defect_based"
           step_defect_based!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.operator_sequence, symmetrized_defect=true)
           err = dt*norm(tsi.psi2)/(tsi.order+1)/tsi.tol
        else
           step_embedded!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.scheme2, tsi.operator_sequence)
           err = distance(tsi.psi, tsi.psi2)/tsi.tol
        end   
        dt = dt*min(facmax, max(facmin, fac*(1.0/err)^(1.0/(tsi.order+1))))
        if err>=1.0
           copy!(tsi.psi, tsi.psi0)
           @printf("t=%17.9e  err=%17.8e  dt=%17.8e  rejected...\n", Float64(state.t), Float64(err), Float64(dt))
        end   
    end
    state.t + dt0, AdaptiveTimeStepperState(state.t+dt0, dt)
end

################################################################################
struct AdaptiveTimeStepper2Iterator
   psi::WaveFunction
   t0::Real
   tend::Real
   dt::Real
   tol::Real
   scheme1
   scheme2
   order::Integer
   operator_sequence
   psi2::WaveFunction
   psi0::WaveFunction
end

struct AdaptiveTimeStepper2State
   t::Real
   dt::Real
   dt_old::Real
   err_old::Real
end   

function adaptive_time_stepper2(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, 
                         scheme1, scheme2, order::Integer, operator_sequence="AB")
    psi2=clone(psi)                     
    psi0=clone(psi)                     
    AdaptiveTimeStepper2Iterator(psi, t0, tend, dt, tol,
        scheme1, scheme2, order, operator_sequence, psi2, psi0)
end        

Base.start(tsi::AdaptiveTimeStepper2Iterator) = AdaptiveTimeStepper2State(tsi.t0, tsi.dt, -1.0, -1.0)

function Base.done(tsi::AdaptiveTimeStepper2Iterator, state::AdaptiveTimeStepper2State)
  state.t >= tsi.tend
end  

function Base.next(tsi::AdaptiveTimeStepper2Iterator, state::AdaptiveTimeStepper2State)
    facmin = 0.25
    facmax = 4.0
    fac = 0.9
    beta1=0.25
    beta2=0.25
    alpha2=0.25

    dt = state.dt
    dt0 = dt
    dt_old=dt
    copy!(tsi.psi0, tsi.psi)
    err = 2.0 #error/tol
    while err>=1.0
        dt = min(dt, tsi.tend-state.t)
        dt0 = dt
        if tsi.scheme2=="palindromic"
           step_palindromic!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.operator_sequence)
           err = 0.5*distance(tsi.psi, tsi.psi2)/tsi.tol
        elseif tsi.scheme2=="defect_based"
           step_defect_based!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.operator_sequence)
           err = dt*norm(tsi.psi2)/(tsi.order+1)/tsi.tol
        elseif tsi.scheme2=="symmetrized_defect_based"
           step_defect_based!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.operator_sequence, symmetrized_defect=true)
           err = dt*norm(tsi.psi2)/(tsi.order+1)/tsi.tol
        else
           step_embedded!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.scheme2, tsi.operator_sequence)
           err = distance(tsi.psi, tsi.psi2)/tsi.tol
        end   
        if state.dt_old<0.0
           dt = dt*min(facmax, max(facmin, fac*(1.0/err)^(1.0/(tsi.order+1))))
        else
           dt = dt*min(facmax, max(facmin, fac*((1.0/err)^(beta1/(tsi.order+1))*(1.0/state.err_old)^(beta2/(tsi.order+1))*(dt/state.dt_old)^(-alpha2))))         
        end
        if err>=1.0
           copy!(tsi.psi, tsi.psi0)
           @printf("t=%17.9e  err=%17.8e  dt=%17.8e  rejected...\n", Float64(state.t), Float64(err), Float64(dt))
        end   
    end
    state.t + dt0, AdaptiveTimeStepper2State(state.t+dt0, dt, dt_old, err)
end



################################################################################
struct EmbeddedScheme
    scheme1
    scheme2
    order :: Integer
end

struct PalindromicScheme
    scheme
    order :: Integer
end

struct DefectBasedScheme
    scheme
    order :: Integer
    symmetrized :: Bool
end

function adaptive_time_stepper(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, es::EmbeddedScheme, 
                         operator_sequence="AB")
    adaptive_time_stepper(psi, t0, tend, dt, tol, es.scheme1, es.scheme2, es.order,
                         operator_sequence)
end                          

function adaptive_time_stepper(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, ps::PalindromicScheme, 
                         operator_sequence="AB")
    adaptive_time_stepper(psi, t0, tend, dt, tol, ps.scheme, "palindromic", ps.order,
                         operator_sequence)
end           

function adaptive_time_stepper(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, ds::DefectBasedScheme, 
                         operator_sequence="AB")
    adaptive_time_stepper(psi, t0, tend, dt, tol, ds.scheme, (ds.symmetrized ? "symmetrized_defect_based" : "defect_based"), ds.order,
                         operator_sequence)
end    


function adaptive_time_stepper2(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, es::EmbeddedScheme, 
                         operator_sequence="AB")
    adaptive_time_stepper2(psi, t0, tend, dt, tol, es.scheme1, es.scheme2, es.order,
                         operator_sequence)
end                          

function adaptive_time_stepper2(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, ps::PalindromicScheme, 
                         operator_sequence="AB")
    adaptive_time_stepper2(psi, t0, tend, dt, tol, ps.scheme, "palindromic", ps.order,
                         operator_sequence)
end         

function adaptive_time_stepper2(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, tol::Real, ds::DefectBasedScheme, 
                         operator_sequence="AB")
    adaptive_time_stepper2(psi, t0, tend, dt, tol, ds.scheme, (ds.symmetrized ? "symmetrized_defect_based" : "defect_based"), ds.order,
                         operator_sequence)
end    
                 

########################################################################


function global_orders(psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real;
                       scheme::Tuple=Strang, operator_sequence="AB", rows=8, corrector::String="none", order::Integer=2)
    @assert psi.m==reference_solution.m
    tab = zeros(Float64, rows, 3)

    wf_save_initial_value = clone(psi)
    copy!(wf_save_initial_value, psi)

    psi2 = clone(psi)

    steps = floor((tend-t0)/dt)
    dt1 = dt
    err_old = 0.0
    println("             dt         err      p")
    println("-----------------------------------")
    for row=1:rows
        for k=1:steps
            if corrector=="palindromic"
                step_palindromic!(psi, psi2, dt1, scheme, operator_sequence)
                axpy!(psi, psi2, 1.0)
                scale!(psi, 0.5)
            elseif corrector=="defect_based"
                step_defect_based!(psi, psi2, dt1, scheme, operator_sequence)
                scale!(psi2, -dt1/(order+1))
                axpy!(psi, psi2, +1.0)
            elseif corrector=="symmetrized_defect_based"
                step_defect_based!(psi, psi2, dt1, scheme, operator_sequence, symmetrized_defect=true)
                scale!(psi2, -dt1/(order+1))
                axpy!(psi, psi2, +1.0)
            else
                step!(psi, dt1, scheme, operator_sequence)
            end  
        end
        err = distance(psi, reference_solution)
        if (row==1) 
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
    tab
end

function local_orders(psi::WaveFunction, t0::Real, dt::Real; 
                      scheme::Tuple=Strang, operator_sequence="AB", 
                      reference_scheme::Tuple=scheme, reference_operator_sequence=operator_sequence,
                      reference_steps=10,
                      rows=8)
    tab = zeros(Float64, rows, 3)

    wf_save_initial_value = clone(psi)
    psi_ref = clone(psi)
    copy!(wf_save_initial_value, psi)

    dt1 = dt
    err_old = 0.0
    println("             dt         err      p")
    println("-----------------------------------")
    for row=1:rows
        step!(psi, dt1, scheme, operator_sequence)
        copy!(psi_ref,wf_save_initial_value)
        for k=1:reference_steps
            step!(psi_ref, dt1/reference_steps, reference_scheme, reference_operator_sequence)
        end    
        err = distance(psi, psi_ref)
        if (row==1) 
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
        copy!(psi,wf_save_initial_value)
    end
    tab
end


function local_orders(psi::WaveFunction, get_reference_solution::Function, 
                       t0::Real, dt::Real; 
                       scheme::Tuple=Strang, operator_sequence="AB", rows=8)
    tab = zeros(Float64, rows, 3)

    wf_save_initial_value = clone(psi)
    reference_solution = clone(psi)
    copy!(wf_save_initial_value, psi)

    dt1 = dt
    println("             dt         err      p")
    println("-----------------------------------")
    for row=1:rows
        step!(psi, dt1, scheme, operator_sequence)
        set!(reference_solution, get_reference_solution, t0+dt1)
        err = distance(psi, reference_solution)
        if (row==1) 
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
        copy!(psi,wf_save_initial_value)
    end
    tab
end


function local_orders_0(psi::WaveFunction, get_reference_solution::Function, 
                       t0::Real, dt::Real, 
                       scheme1, scheme2; operator_sequence="AB", rows=8, order=0)
    tab = zeros(Float64, rows, 5)

    wf_save_initial_value = clone(psi)
    psi2 = clone(psi)
    reference_solution = clone(psi)
    copy!(wf_save_initial_value, psi)

    dt1 = dt
    err_old1 = 0.0 
    err_old2 = 0.0
    println("             dt         err      p         err      p")
    println("------------------------------------------------------")
    for row=1:rows
        if scheme2=="palindromic"
            step_palindromic!(psi, psi2, dt1, scheme1, operator_sequence)
            axpy!(psi2, psi, 1.0)
            scale!(psi2, 0.5)
        elseif scheme2=="defect_based"
            step_defect_based!(psi, psi2, dt1, scheme1, operator_sequence)
            scale!(psi2, -dt1/(order+1))
            axpy!(psi2, psi, +1.0)
        elseif scheme2=="symmetrized_defect_based"
            step_defect_based!(psi, psi2, dt1, scheme1, operator_sequence, symmetrized_defect=true)
            scale!(psi2, -dt1/(order+1))
            axpy!(psi2, psi, +1.0)
        else
            step_embedded!(psi, psi2, dt1, scheme1, scheme2 , operator_sequence)
        end
        set!(reference_solution, get_reference_solution, t0+dt1)
        err1 = distance(psi, reference_solution)
        err2 = distance(psi2, reference_solution)
        if (row==1) 
            @printf("%3i%12.3e%12.3e       %12.3e\n", row, Float64(dt1), Float64(err1), Float64(err2))
            tab[row,1] = dt1
            tab[row,2] = err1
            tab[row,3] = 0 
            tab[row,4] = err2
            tab[row,5] = 0 
        else
            p1 = log(err_old1/err1)/log(2.0);
            p2 = log(err_old2/err2)/log(2.0);
            @printf("%3i%12.3e%12.3e%7.2f%12.3e%7.2f\n", row, Float64(dt1), Float64(err1), Float64(p1), Float64(err2), Float64(p2))
            tab[row,1] = dt1
            tab[row,2] = err1
            tab[row,3] = p1
            tab[row,4] = err2
            tab[row,5] = p2 
 
        end
        err_old1 = err1
        err_old2 = err2
        dt1 = 0.5*dt1
        copy!(psi,wf_save_initial_value)
    end
    tab
end

function local_orders(psi::WaveFunction, get_reference_solution::Function, 
                       t0::Real, dt::Real, 
                       embedded_scheme::EmbeddedScheme; operator_sequence="AB", rows=8)
    local_orders_0(psi, get_reference_solution, t0, dt,  embedded_scheme.scheme1, embedded_scheme.scheme2,
                   operator_sequence=operator_sequence, rows=rows)
end   

function local_orders(psi::WaveFunction, get_reference_solution::Function, 
                       t0::Real, dt::Real, 
                       palindromic_scheme::PalindromicScheme; operator_sequence="AB", rows=8)
    local_orders_0(psi, get_reference_solution, t0, dt,  palindromic_scheme.scheme, "palindromic",
                   operator_sequence=operator_sequence, rows=rows)
end    

function local_orders(psi::WaveFunction, get_reference_solution::Function, 
                       t0::Real, dt::Real, 
                       defect_based_scheme::DefectBasedScheme; operator_sequence="AB", rows=8)
    local_orders_0(psi, get_reference_solution, t0, dt, defect_based_scheme.scheme, 
                   (defect_based_scheme.symmetrized ? "symmetrized_defect_based" : "defect_based"),
                   operator_sequence=operator_sequence, rows=rows, order=defect_based_scheme.order)
end    
