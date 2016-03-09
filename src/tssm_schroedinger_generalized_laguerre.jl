# T, TSSM_HANDLE
println("including tssm_schroedinger_generalized_laguerre.jl for type ", T)
for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:SchroedingerGeneralizedLaguerre2D, :_schroedinger_gen_laguerre_2d, true, 2 ),             
                 (:SchroedingerGeneralizedLaguerreHermite3D, :_schroedinger_gen_laguerre_hermite_3d, true, 3 ),             
                 (:SchroedingerGeneralizedLaguerreReal2D, :_schroedinger_gen_laguerre_real_2d, false, 2 ),
                 (:SchroedingerGeneralizedLaguerreHermiteReal3D, :_schroedinger_gen_laguerre_hermite_real_3d, false, 3 ),
                )
println("    ", METHOD) 
if T == :Float128
    SUF = symbol(SUF, "_wf128")
end

if DIM==2
  @eval begin
    function ($METHOD)(T::Type{($T)}, ntheta::Integer, nfr::Integer,
                                      omega_r::Real, Omega::Real,
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_3D
        with_potential_t = potential_t!=none_4D
        V_c = cfunction(potential, ($T), (($T),($T)))
        V_t_c = cfunction(potential_t, ($T), (($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Void}, 
                   (Int32, Int32, ($T), ($T), 
                   ($T), ($T), Ptr{Void}, Bool, Ptr{Void}, Bool, ($T), Int32), 
                   ntheta, nfr, omega_r, Omega,  
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t,
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string(PRE,"finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
  end # eval
  if T == :Float64
  @eval begin
    function ($METHOD)(ntheta::Integer, nfr::Integer,
                       omega_r::Real, Omega::Real,
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), ntheta, nfr, omega_r, Omega,
                       hbar=hbar, mass=mass, potential=potential,
                       potential_t=potential_t,
                       cubic_coupling=cubic_coupling,
                       boundary_conditions=boundary_conditions)
    end      
  end # eval
  end # if

elseif DIM==3
   @eval begin
     function ($METHOD)(T::Type{($T)}, ntheta::Integer, nfr::Integer,
                                      omega_r::Real, Omega::Real,
                                      nz::Integer, omega_z::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_3D
        with_potential_t = potential_t!=none_4D
        V_c = cfunction(potential, ($T), (($T),($T),($T)))
        V_t_c = cfunction(potential_t, ($T), (($T),($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Void}, 
                   (Int32, Int32, ($T), ($T), Int32, ($T), 
                   ($T), ($T), Ptr{Void}, Bool, Ptr{Void}, Bool, ($T), Int32), 
                   ntheta, nfr, omega_r, Omega, nz, omega_z,  
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t,
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string(PRE,"finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
     end
   end # eval
   if T == :Float64
   @eval begin
    function ($METHOD)(ntheta::Integer, nfr::Integer,
                       omega_r::Real, Omega::Real,
                       nz::Integer, omega_z::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), ntheta, nfr, omega_r, Omega, nz, omega_z,
                       hbar=hbar, mass=mass, potential=potential,
                       potential_t=potential_t,
                       cubic_coupling=cubic_coupling,
                       boundary_conditions=boundary_conditions)
    end      
  end # eval
  end #if
 
end # if

@eval begin
    function get_ntheta(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ntheta",SUF))), Int32,
             (Ptr{Void}, ), m.m )
    end

    function get_nfr(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nfr",SUF))), Int32,
             (Ptr{Void}, ), m.m )
    end

    function get_nr(m::($METHOD){$T})
      ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nr",SUF))), Int32,
            (Ptr{Void}, ), m.m )    
    end

    function get_omega_r(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_omega_r",SUF))), ($T),
             (Ptr{Void}, ), m.m )
    end

    function get_Omega(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_Omega",SUF))), ($T),
             (Ptr{Void}, ), m.m )
    end
end # eval    

if DIM>=3 
    @eval begin
        function get_nz(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nz",SUF))), Int32,
                 (Ptr{Void}, ), m.m )
        end

        function get_omega_z(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_omega_z",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end
    end # eval    
end 

if DIM==2
    @eval begin

        function get_weights(m::($METHOD){$T})
           dim =Array(Int32, 1)
           np = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights_r",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}), m.m, dim)
           n = pointer_to_array(np, dim[1], false)     
           copy(n)
        end

        function get_transformation_matrices(m::($METHOD){$T}, unsafe_access::Bool=false)
            dims =Array(Int32, 3)
            Lp = ccall( Libdl.dlsym(tssm_handle, $(string(PRE,"get_L",SUF))), Ptr{$T},
                (Ptr{Void}, Ptr{Int32}), m.m, dims )
            L = pointer_to_array(Lp, (dims[1], dims[2], dims[3]), false)  
            if unsafe_access
                return L
            else
                return copy(L)
            end
        end
        

    end # eval    

elseif DIM==3
    @eval begin

        function get_weights(m::($METHOD){$T})
           dim =Array(Int32, 1)
           np1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights_r",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}), m.m, dim)
           n1 = pointer_to_array(n1p, dim[1], false)     
           np2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights_z",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}), m.m, dim)
           n2 = pointer_to_array(n2p, dim[1], false)     
           copy(n1), copy(n2)
        end

        function get_transformation_matrices(m::($METHOD){$T}, unsafe_access::Bool=false)
            dims =Array(Int32, 3)
            Lp = ccall( Libdl.dlsym(tssm_handle, $(string(PRE,"get_L",SUF))), Ptr{$T},
                (Ptr{Void}, Ptr{Int32}), m.m, dims )
            H = pointer_to_array(Hp, (dims[1], dims[2], dims[3]), false)  
            Hp = ccall( Libdl.dlsym(tssm_handle, $(string(PRE,"get_H_z",SUF))), Ptr{$T},
                (Ptr{Void}, Ptr{Int32}), m.m, dims )
            H = pointer_to_array(Hp, (dims[1], dims[2]), false)  
            if unsafe_access
                return L, H
            else
                return copy(L), copy(H)
            end
        end

    end # eval    
end #if

end # for

