
Strang = (0.5, 10.0, 0.5)

function step!(psi::WaveFunction, dt::Real, scheme=Strang, operator_sequence="AB")
    for k = 1:length(scheme)
        which_operator = operator_sequence[mod(k-1, length(operator_sequence)) + 1]
        if which_operator == 'A'
            propagate_A!(psi, dt*scheme[k]) 
        elseif which_operator == 'B'
            propagate_B!(psi, dt*scheme[k]) 
        elseif which_operator == 'C'
            propagate_C!(psi, dt*scheme[k]) 
        else
            # TODO! Error or Warning
        end
    end 
end

immutable EquidistantTimeStepperIterator
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
 

immutable AdaptiveTimeStepperIterator
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

immutable AdaptiveTimeStepperState
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

function Base.next(tsi::AdaptiveTimeStepperIterator, state::AdaptiveTimeStepperState)
    const facmin = 0.25
    const facmax = 4.0
    const fac = 0.9

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
        else
           step_embedded!(tsi.psi, tsi.psi2, dt, tsi.scheme1, tsi.scheme2, tsi.operator_sequence)
           err = distance(tsi.psi, tsi.psi2)/tsi.tol
        end   
        dt = dt*min(facmax, max(facmin, fac*(1.0/err)^(1.0/(tsi.order+1))))
        if err>=1.0
           copy!(tsi.psi, tsi.psi0)
           @printf("t=%17.9e  err=%17.8e  dt=%17.8e  rejected...\n", state.t, err, dt)
        end   
    end
    state.t + dt0, AdaptiveTimeStepperState(state.t+dt0, dt)
end

################################################################################
immutable AdaptiveTimeStepper2Iterator
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

immutable AdaptiveTimeStepper2State
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
    const facmin = 0.25
    const facmax = 4.0
    const fac = 0.9
    const beta1=0.25
    const beta2=0.25
    const alpha2=0.25
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
           @printf("t=%17.9e  err=%17.8e  dt=%17.8e  rejected...\n", state.t, err, dt)
        end   
    end
    state.t + dt0, AdaptiveTimeStepper2State(state.t+dt0, dt, dt_old, err)
end



################################################################################
immutable EmbeddedScheme
    scheme1
    scheme2
    order :: Integer
end

immutable PalindromicScheme
    scheme
    order :: Integer
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

########################################################################


function global_orders(psi::WaveFunction, reference_solution::WaveFunction, 
                       t0::Real, tend::Real, dt::Real, 
                       scheme::Tuple=Strang, operator_sequence="AB", rows=8)
    @assert psi.m==reference_solution.m
    tab = Array(Float64, rows, 3)

    wf_save_initial_value = clone(psi)
    copy!(wf_save_initial_value, psi)

    steps = floor((tend-t0)/dt)
    dt1 = dt
    println("             dt         err      p")
    println("-----------------------------------")
    for row=1:rows
        for k=1:steps
           step!(psi, dt1, scheme, operator_sequence)
        end   
        err = distance(psi, reference_solution)
        if (row==1) then
            @printf("%3i%12.3e%12.3e\n", row, dt1, err)
            tab[row,1] = dt1
            tab[row,2] = err
            tab[row,3] = 0 
        else
            p = log(err_old/err)/log(2.0);
            @printf("%3i%12.3e%12.3e%7.2f\n", row, dt1, err, p)
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

function local_orders(psi::WaveFunction, get_reference_solution, 
                       t0::Real, dt::Real, 
                       scheme::Tuple=Strang, operator_sequence="AB", rows=8)
    tab = Array(Float64, rows, 3)

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
        if (row==1) then
            @printf("%3i%12.3e%12.3e\n", row, dt1, err)
            tab[row,1] = dt1
            tab[row,2] = err
            tab[row,3] = 0 
        else
            p = log(err_old/err)/log(2.0);
            @printf("%3i%12.3e%12.3e%7.2f\n", row, dt1, err, p)
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


function local_orders_0(psi::WaveFunction, get_reference_solution, 
                       t0::Real, dt::Real, 
                       scheme1, scheme2, operator_sequence="AB", rows=8)
    tab = Array(Float64, rows, 5)

    wf_save_initial_value = clone(psi)
    psi2 = clone(psi)
    reference_solution = clone(psi)
    copy!(wf_save_initial_value, psi)

    dt1 = dt
    println("             dt         err       p         err      p")
    println("-------------------------------------------------------")
    for row=1:rows
        if scheme2=="palindromic"
            step_palindromic!(psi, psi2, dt1, scheme1, operator_sequence)
            axpy!(psi2, psi, 1.0)
            scale!(psi2, 0.5)
        else
            step_embedded!(psi, psi2, dt1, scheme1, scheme2 , operator_sequence)
        end
        set!(reference_solution, get_reference_solution, t0+dt1)
        err1 = distance(psi, reference_solution)
        err2 = distance(psi2, reference_solution)
        if (row==1) then
            @printf("%3i%12.3e%12.3e        %12.3e\n", row, dt1, err1, err2)
            tab[row,1] = dt1
            tab[row,2] = err1
            tab[row,3] = 0 
            tab[row,4] = err2
            tab[row,5] = 0 
        else
            p1 = log(err_old1/err1)/log(2.0);
            p2 = log(err_old2/err2)/log(2.0);
            @printf("%3i%12.3e%12.3e%7.2fe%12.3e%7.2f\n", row, dt1, err1, p1, err2, p2)
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

function local_orders(psi::WaveFunction, get_reference_solution, 
                       t0::Real, dt::Real, 
                       embedded_scheme::EmbeddedScheme, operator_sequence="AB", rows=8)
    local_orders_0(psi, get_reference_solution, t0, dt,  embedded_scheme.scheme1, embedded_scheme.scheme2,
                 operator_sequence, rows)
end   

function local_orders(psi::WaveFunction, get_reference_solution, 
                       t0::Real, dt::Real, 
                       palindromic_scheme::PalindromicScheme, operator_sequence="AB", rows=8)
    local_orders_0(psi, get_reference_solution, t0, dt,  palindromic_scheme.scheme, "palindromic",
                 operator_sequence, rows)
end    
