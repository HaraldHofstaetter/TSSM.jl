Strang=(0.5, 1, 0.5)

function step_imaginary_time!(this::WaveFunction, dt::Real, scheme=Strang, operator_sequence="AB")
    for k = 1:length(scheme)
        which_operator = operator_sequence[mod(k, length(operator_sequence)) + 1]
        if which_operator == 'A'
            imaginary_time_propagate_A!(this, dt*scheme[k]) 
        elseif which_operator == 'B'
            imaginary_time_propagate_B!(this, dt*scheme[k]) 
        elseif which_operator == 'C'
            imaginary_time_propagate_C!(this, dt*scheme[k]) 
        else
            # TODO! ERror or Warning
        end
    end 
end

function step_imaginary_time_extrapolated!(this::WaveFunction, dt::Real, extrapolation_order, operator_sequence, psi0, psi1)
    if extrapolation_order==1
        c = ( 1.0, )
    elseif extrapolation_order==2
        c = ( -1.0/3.0, 4.0/3.0 )
    elseif extrapolation_order==3
        c = ( 1.0/24.0, -16.0/15.0, 81.0/40.0 )
    elseif extrapolation_order==4
        c = ( -1.0/360.0,  16.0/45.0, -729.0/280.0, 1024.0/315.0 )
    elseif extrapolation_order==5
        c = ( 1.0/8640.0, -64.0/945.0, 6561.0/4480.0, -16384.0/2835.0, 390625.0/72576.0 )
    elseif extrapolation_order==6
        c = ( -1.0/302400.0, 8.0/945.0, -2187.0/4480.0, 65536.0/14175.0, -9765625.0/798336.0, 17496.0/1925.0 )
    else
        # TODO: Error or Warning
    end    

    if operator_sequence[1]=='A'
       to_frequency_space!(this)
       copy!(psi0, this)
       copy!(psi1, this)

       imaginary_time_propagate_A!(this, 0.5*dt)
       imaginary_time_propagate_B!(this, dt, 2)
       imaginary_time_propagate_A!(this, 0.5*dt)
       
       scale!(this, c[1])

       for k = 2:extrapolation_order
           dt1 = dt/k
           copy!(psi1, psi0)
           imaginary_time_propagate_A!(psi1, 0.5*dt1)
           for j = 1:k-1
               imaginary_time_propagate_B!(psi1, dt1, 2)
               imaginary_time_propagate_A!(psi1, dt1)
           end 
           imaginary_time_propagate_B!(psi1, dt1, 2)
           imaginary_time_propagate_A!(psi1, 0.5*dt1)

           axpy!(this, psi1, c[k])
       end 
    else
       to_real_space!(this)
       copy!(psi0, this)
       copy!(psi1, this)

       imaginary_time_propagate_B!(this, 0.5*dt, 2)
       imaginary_time_propagate_A!(this, dt)
       imaginary_time_propagate_B!(this, 0.5*dt, 2)
       
       scale!(this, c[1])

       for k = 2:extrapolation_order
           dt1 = dt/k
           copy!(psi1, psi0)
           imaginary_time_propagate_B!(psi1, 0.5*dt1, 2)
           for j = 1:k-1
               imaginary_time_propagate_A!(psi1, dt1)
               imaginary_time_propagate_B!(psi1, dt1, 2)
           end 
           imaginary_time_propagate_A!(psi1, dt1, 2)
           imaginary_time_propagate_B!(psi1, 0.5*dt1, 2)

           axpy!(this, psi1, c[k])
           this%u = this%u + c(k)*psi1%u
       end 
    end   
end

#########################################

function  groundstate!(this::WaveFunction; dt::Real=0.05, 
                       tol::Real=1e-8, max_iters::Integer=10000, 
                       extrapolation_order::Integer=2, scheme=Strang, operator_sequence="AB")
    time0 =time()
    if dim(this)==1    
        init(x) = 1.0
    elseif dim(this)==2    
        init(x,y) = 1.0
    elseif dim(this)==3    
        init(x,y,z) = 1.0
    end     

   psi1 = clone(this)
   if extrapolation_order >= 2
       psi2 = clone(this)
       psi3 = clone(this)
   end 
   psi_old = clone(this) 

   set!(this, init)

   copy!(psi_old, this)

   k = 0
   k_check = 20
   E_old = 1.0e6 #TODO: infinity

   for k=0:max_iters
       if mod(k, k_check)==0
           copy!(psi1, this)

           if extrapolation_order==1 then
               step_imaginary_time!(psi1, 0.5*dt, scheme, operator_sequence)
           else
               step_imaginary_time_extrapolated!(psi1, 0.5*dt, extrapolation_order, operator_sequence, psi2, psi3,)
           end
#ifndef _REAL_
#           to_recompute_groundstate(psi, extrapolation_order=2)al_space(psi1)
#           psi1%u = real(psi1%u, prec)
#endif 
println(norm(psi1))
           normalize!(psi1)
           E_mu1, E_dev1 = get_energy_expectation_deviation(psi1)
           E1 = E_mu1 - interaction_energy(psi1)
           err1 = E_dev1/E1
       end 

       if extrapolation_order==1 then
           step_imaginary_time!(this, dt, scheme, operator_sequence)
       else
           step_imaginary_time_extrapolated!(this, dt, extrapolation_order, operator_sequence, psi2, psi3,)
       end 

#ifndef _REAL_
#       call this%to_real_space()
#       this%u = real(this%u, prec)
#endif 
       normalize!(this)           
       E_mu, E_dev = get_energy_expectation_deviation(this)
       E = E_mu - interaction_energy(this)
       err = E_dev/E

       ddd = distance(this, psi_old)
       calc_time = time() - time0
       if mod(k,k_check)==0 then
           @printf("%5i%24.15e%24.15e%12.3e%12.3e%12.3e%10.2f%24.15e%12.3e%12.3e\n",
                    k, E, E_mu,  E_old-E, err, ddd, calc_time, E1, E_old-E1, err1)
           if err1<err then
               @printf("changed step size, old:%24.15e  new:%24.15e\n", dt, 0.5*dt)
                 
               dt = 0.5*dt
               copy!(this,psi1)
               #Note that if extrapolation_step>1, then psi1 has been overwritten and
               #cannot be used anymore
           end 
       else    
          @printf("%5i%24.15e%24.15e%12.3e%12.3e%12.3e%10.2f\n", 
                   k, E, E_mu, E_old-E,  err, ddd, calc_time)
       end 

       if (err<tol) 
          break
       end   

       E_old = E
       copy!(psi_old,this)
        
    end #for 
end




