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
    scale!(m.k1, -0.5im*dt)
    axpy!(m.k1, psi, 1.0)
    set_time!(m.k1, t) 
    gen_rhs!(m.k2, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    axpy!(psi, m.k2, -1im*dt)
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
    t4 = freeze_time ? t : t+t
    #For this storage-efficient implementation, see
    #E. K. Blum: A  Modification of the  Runge-Kutta Fourth-Order Method
    gen_rhs!(m.k3, psi, include_kinetic_part=include_kinetic_part,
                        include_one_particle_potential_part=include_one_particle_potential_part)
    scale!(m.k3, -1im*dt)
    copy!(m.k2, m.k3)
    copy!(m.k1, psi)
    axpy!(m.k1, m.k3, 0.5)
    set_time!(m.k1, t2) 
    gen_rhs!(m.k3, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    scale!(m.k3, -1im*dt)
    axpy!(m.k1, m.k3, 0.5)
    axpy!(m.k1, m.k2, -0.5)
    scale!(m.k2,1.0/6)
    scale!(m.k3, -0.5)
    set_time!(m.k1, t3) 
    gen_rhs!(m.k4, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    axpy!(m.k3, m.k4, -1im*dt)
    axpy!(m.k1, m.k3, 1.0)
    axpy!(m.k2, m.k3, -1.0)
    scale!(m.k3, 2.0)
    set_time!(m.k1, t4)     
    gen_rhs!(m.k4, m.k1, include_kinetic_part=include_kinetic_part,
                         include_one_particle_potential_part=include_one_particle_potential_part)
    axpy!(m.k3, m.k4, -1im*dt)
    axpy!(m.k1, m.k2, 1.0)
    axpy!(m.k1, m.k3, 1.0/6)
    copy!(psi,m.k1)
    set_time!(psi, t4) 
    orthonormalize_orbitals!(psi)
end




function groundstate!(psi::WfMCTDHF1D, dt::Real, n::Int; output_step::Int=1, 
                      keep_initial_value::Bool=false)
    m = psi.m
    m.k1 = wave_function(m)
    m.k2 = wave_function(m)

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
        psi.a[:] = ones(m.lena)
    end

    orthonormalize_orbitals!(psi)
    psi.a[:] = psi.a[:]/Base.norm(psi.a)
    time0 = time()

    for k=1:n
        imaginary_time_propagate_A!(psi, 0.5*dt)
        orthonormalize_orbitals!(psi)
        RK2_step!(psi, -1im*dt)
        orthonormalize_orbitals!(psi)
        imaginary_time_propagate_A!(psi, 0.5*dt)
        orthonormalize_orbitals!(psi)
        norm_psi = norm(psi)
        psi.a[:] *= 1/norm_psi
        
        if mod(k,output_step)==0
            E_pot = potential_energy(psi)
            E_kin = kinetic_energy(psi)
            E = E_pot + E_kin
            ctime = time() - time0
            @printf("step=%4i  E_pot=%14.10f  E_kin=%14.10f  E=%14.10f ctime=%10.2f\n", k, E_pot, E_kin, E, ctime)            
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

    
function run!(psi::WfMCTDHF1D, dt::Real, n::Int; output_step::Int=1, include_kinetic_part_in_B::Bool=false, 
    include_one_particle_potential_part_in_B::Bool=true)
    m = psi.m
    m.k1 = wave_function(m)
    m.k2 = wave_function(m)
    time0 = time()
    set_propagate_time_together_with_A!(m, true)

    orthonormalize_orbitals!(psi)

    for k=1:n
        strang_step!(psi, dt, include_kinetic_part_in_B=include_kinetic_part_in_B,
            include_one_particle_potential_part_in_B=include_one_particle_potential_part_in_B)
        
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
