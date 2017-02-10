function RK2_step!(psi::WfMCTDHF1D, dt::Number; include_kinetic_part::Bool=false, 
    include_one_particle_potential_part::Bool=true, freeze_time::Bool=!include_kinetic_part)
    # Usually time is frozen if the kinetic part is not included, because in this case
    # time is propagated by the A (=kinetic) part of the time splitting method.
    m = psi.m
    t = get_time(psi)
    t1 = freeze_time ? t : t+0.5*dt
    t2 = freeze_time ? t : t+dt
    gen_rhs!(m.k1, psi, include_kinetic_part=include_kinetic_part,
                        include_one_particle_potential_part=include_one_particle_potential_part)
    scale!(m.k1, 0.5*dt)
    axpy!(m.k1, psi, 1.0)
    set_time!(m.k1, t) 
    gen_rhs!(m.k2, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    axpy!(psi, m.k2, dt)
    set_time!(psi, t2) 
    orthonormalize_orbitals!(psi)
end

function RK4_step!(psi::WfMCTDHF1D, dt::Number; include_kinetic_part::Bool=false, 
    include_one_particle_potential_part::Bool=true, freeze_time::Bool=!include_kinetic_part)
    # Usually time is frozen if the kinetic part is not included, because in this case
    # time is propagated by the A (=kinetic) part of the time splitting method.
    m = psi.m
    t = get_time(psi)
    t1 = freeze_time ? t : t
    t2 = freeze_time ? t : t+0.5*dt
    t3 = freeze_time ? t : t+0.5*dt
    t4 = freeze_time ? t : t+dt
    #For this storage-efficient implementation, see
    #E. K. Blum: A  Modification of the  Runge-Kutta Fourth-Order Method
    gen_rhs!(m.k3, psi, include_kinetic_part=include_kinetic_part,
                        include_one_particle_potential_part=include_one_particle_potential_part)
    scale!(m.k3, dt)
    copy!(m.k2, m.k3)
    copy!(m.k1, psi)
    axpy!(m.k1, m.k3, 0.5)
    set_time!(m.k1, t2) 
    gen_rhs!(m.k3, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    scale!(m.k3, dt)
    axpy!(m.k1, m.k3, 0.5)
    axpy!(m.k1, m.k2, -0.5)
    scale!(m.k2,1.0/6)
    scale!(m.k3, -0.5)
    set_time!(m.k1, t3) 
    gen_rhs!(m.k4, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    axpy!(m.k3, m.k4, dt)
    axpy!(m.k1, m.k3, 1.0)
    axpy!(m.k2, m.k3, -1.0)
    scale!(m.k3, 2.0)
    set_time!(m.k1, t4)     
    gen_rhs!(m.k4, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    axpy!(m.k3, m.k4, dt)
    axpy!(m.k1, m.k2, 1.0)
    axpy!(m.k1, m.k3, 1.0/6)
    copy!(psi,m.k1)
    set_time!(psi, t4) 
    orthonormalize_orbitals!(psi)
end




function groundstate!(psi::WfMCTDHF1D;  dt::Float64=0.05, max_iter::Int=2000, output_step::Int=20, tol=1e-5,
                      keep_initial_value::Bool=false)
    m = psi.m
    m.k1 = wave_function(m)
    m.k2 = wave_function(m)
    psi1 = wave_function(m)
    

    if !keep_initial_value
        to_frequency_space!(psi)
        nx = get_nx(m.m)
        for k=1:m.N
            u=get_data(psi.o[k].phi,true)
            u[:]=zeros(Complex{Float64}, length(u))
            u[div(k+1,2)]=1
            if k>=3
                u[nx-div(k+1,2)+2]=1
            end
        end
        to_real_space!(psi)
        psi.a[:] = zeros(m.lena)
        for k=1:m.lena #generate singlet state (works only for f=2!!!)
            (i1,i2) = m.slater_indices[k]            
            if isodd(i1) && i2==i1+1
                psi.a[k] = 1.0
            end
        end
    end

    orthonormalize_orbitals!(psi)
    psi.a[:] = psi.a[:]/Base.norm(psi.a)
    time0 = time()

    for k=1:max_iter
        if mod(k,output_step)==0
            copy!(psi1, psi)
            imaginary_time_propagate_A!(psi1, 0.25*dt)
            orthonormalize_orbitals!(psi1)
            RK2_step!(psi1, -0.5im*dt)
            orthonormalize_orbitals!(psi1)
            imaginary_time_propagate_A!(psi1, 0.25*dt)
            normalize!(psi1)
            E1, dev1 = get_energy_expectation_deviation(psi1)
            err1 = dev1/max(abs(E1),1.0)
        end
        imaginary_time_propagate_A!(psi, 0.5*dt)
        orthonormalize_orbitals!(psi)
        RK2_step!(psi, -1im*dt)
        orthonormalize_orbitals!(psi)
        imaginary_time_propagate_A!(psi, 0.5*dt)
        normalize!(psi)
        
        if mod(k,output_step)==0
            #E_pot = potential_energy(psi)
            #E_kin = kinetic_energy(psi)
            #E = E_pot + E_kin
            #ctime = time() - time0
            #@printf("step=%4i  E_pot=%14.10f  E_kin=%14.10f  E=%14.10f ctime=%10.2f\n", k, E_pot, E_kin, E, ctime)               
            E, dev = get_energy_expectation_deviation(psi)
            err = dev/max(abs(E),1.0)
            ctime = time() - time0
            @printf("step=%4i  E=%14.10f  err=%12.3e  E1=%14.10f  err1=%12.3e  ctime=%10.2f\n", k, E, err, E1, err1, ctime)
            if err1<err 
               @printf("changed step size, old:%24.15e  new:%24.15e\n", dt, 0.5*dt)                 
               dt = 0.5*dt
               copy!(psi,psi1)
            end 
            if (err<tol)||(err1<tol)
                break;
            end
        end
        
    end
    
    m.k1 = nothing
    m.k2 = nothing
end


function strang_step!(psi::WfMCTDHF1D, dt::Real; include_kinetic_part_in_B::Bool=false, 
    include_one_particle_potential_part_in_B::Bool=true)
    if !include_kinetic_part_in_B
        propagate_A!(psi, 0.25*dt)
    end    
    if !include_one_particle_potential_part_in_B
        propagate_B!(psi, 0.5*dt)
    end    
    if !include_kinetic_part_in_B
        propagate_A!(psi, 0.25*dt)
    end       
    RK2_step!(psi, dt, include_kinetic_part=include_kinetic_part_in_B,
        include_one_particle_potential_part=include_one_particle_potential_part_in_B)    
    if !include_kinetic_part_in_B
        propagate_A!(psi, 0.25*dt)
    end   
    if !include_one_particle_potential_part_in_B
        propagate_B!(psi, 0.5*dt)
    end        
    if !include_kinetic_part_in_B
        propagate_A!(psi, 0.25*dt)
    end
end

    
function run!(psi::WfMCTDHF1D, dt::Real, n::Int, method::Function; output_step::Int=1) 
    m = psi.m
    m.k1 = wave_function(m)
    m.k2 = wave_function(m)
    m.k3 = wave_function(m)
    m.k4 = wave_function(m)
    time0 = time()
    set_propagate_time_together_with_A!(m, true)

    orthonormalize_orbitals!(psi)

    for k=1:n
        method(psi, dt)
        
        if mod(k,output_step)==0
            t = get_time(psi)
            nn = norm(psi)
            E_pot = potential_energy(psi)
            E_kin = kinetic_energy(psi)
            E = E_pot + E_kin
            ctime = time() - time0
            @printf("step=%5i  t=%14.10f  norm=%14.10f  E_pot=%14.10f  E_kin=%14.10f  E=%14.10f  ctime=%10.2f\n", k, t, nn, E_pot, E_kin, E, ctime)            
        end
    end
    
    m.k1 = nothing
    m.k2 = nothing
    m.k3 = nothing
    m.k4 = nothing
end


function local_orders(psi::WfMCTDHF1D, dt::Real, method::Function; method_ref::Function=method,                      
                      reference_steps=10,
                      rows=8)
    m = psi.m
    m.k1 = wave_function(m)
    m.k2 = wave_function(m)
    m.k3 = wave_function(m)
    m.k4 = wave_function(m)   
    tab = Array(Float64, rows, 7)
    set_propagate_time_together_with_A!(m, true)
    wf_save_initial_value = wave_function(m)    
    psi_ref = wave_function(m)
    copy!(wf_save_initial_value, psi)

    dt1 = dt
    err_a_old = 0.0
    err_phi_old = 0.0
    err_old = 0.0
    println("             dt       err_a      p      err_phi      p          err      p")
    println("----------------------------------------------------------------------------")
    for row=1:rows
        method(psi, dt1)
        copy!(psi_ref,wf_save_initial_value)
        for k=1:reference_steps
            method_ref(psi_ref, dt1/reference_steps)
        end
        err_phi = -1.0
        for j=1:m.N
            err_phi=max(err_phi, TSSM.distance(psi.o[j].phi, psi_ref.o[j].phi))
        end
        err_a = norm(psi.a-psi_ref.a)
        err = abs(distance(psi, psi_ref))
        if (row==1) then
            @printf("%3i%12.3e%12.3e %19.3e %19.3e\n", row, Float64(dt1), Float64(err_a),Float64(err_phi),Float64(err))
            tab[row,1] = dt1
            tab[row,2] = err_a
            tab[row,3] = 0
            tab[row,4] = err_phi
            tab[row,5] = 0
            tab[row,6] = err
            tab[row,7] = 0            
        else
            p_a = log(err_a_old/err_a)/log(2.0);
            p_phi = log(err_phi_old/err_phi)/log(2.0);
            p = log(err_old/err)/log(2.0);            
            @printf("%3i%12.3e%12.3e%7.2f %12.3e%7.2f %12.3e%7.2f\n", row, Float64(dt1), 
                Float64(err_a), Float64(p_a), Float64(err_phi), Float64(p_phi), Float64(err), Float64(p))
            tab[row,1] = dt1
            tab[row,2] = err_a
            tab[row,3] = p_a
            tab[row,4] = err_phi
            tab[row,5] = p_phi
            tab[row,6] = err
            tab[row,7] = p                        
        end
        err_a_old = err_a
        err_phi_old = err_phi
        err_old = err
        dt1 = 0.5*dt1
        copy!(psi,wf_save_initial_value)
    end
    m.k1 = nothing
    m.k2 = nothing
    m.k3 = nothing
    m.k4 = nothing   
    tab
end



function run!(psi::WfMCTDHF1D, dt::Real, steps::Int, a::Vector{Float64}, b::Vector{Float64};
    method::Symbol=:RK2, iters::Int=0, output_step::Int=1, return_solutions::Bool=false) 
    time0 = time()    
    s = length(a)
    m = psi.m
    if method==:RK2 || method==:RK4
        m.k1 = wave_function(m)
        m.k2 = wave_function(m)
    end
    if method==:RK4
        m.k3 = wave_function(m)
        m.k4 = wave_function(m)
    end
    if iters>0
        psi0 = wave_function(m)
    end
    if return_solutions
        psi_back = WaveFunction[wave_function(m) for j=1:steps+1]
        copy!(psi_back[1],psi)
    end
    set_propagate_time_together_with_A!(m, true)

    orthonormalize_orbitals!(psi)

    for step=1:steps
        for j = 1:s
            if a[j]!=0.0
                propagate_A!(psi, a[j]*dt)            
            end
            if b[j]!=0.0
                if iters>0
                    copy!(psi0, psi)
                end            
                if method == :RK2
                    RK2_step!(psi, b[j]*dt)
                elseif method ==:RK4
                    RK4_step!(psi, b[j]*dt)
                end
                
                for it=1:iters
                    copy!(m.k1, psi0)
                    scale!(m.k1, 0.5)
                    axpy!(m.k1, psi, 0.5)
                    copy!(psi, psi0)
                    set_time!(m.k1, get_time(psi))
                    gen_rhs!(m.k2, m.k1)
                    axpy!(psi, m.k2, b[j]*dt)                    
                end
            end 
            if return_solutions
                copy!(psi_back[step+1], psi) 
            end
        end
        
        if mod(step,output_step)==0
            t = get_time(psi)
            nn = norm(psi)
            E_pot = potential_energy(psi)
            E_kin = kinetic_energy(psi)
            E = E_pot + E_kin
            ctime = time() - time0
            @printf("step=%5i  t=%14.10f  norm=%14.10f  E_pot=%14.10f  E_kin=%14.10f  E=%14.10f  ctime=%10.2f\n", 
                     step, t, nn, E_pot, E_kin, E, ctime)            
        end
    end
    
    m.k1 = nothing
    m.k2 = nothing
    m.k3 = nothing
    m.k4 = nothing
    if return_solutions
        return psi_back
    else
        return nothing
    end
end
