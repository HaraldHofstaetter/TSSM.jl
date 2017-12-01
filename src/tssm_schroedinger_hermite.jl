# T , TSSM_HANDLE 
println("including tssm_schroedinger_hermite.jl for type ", T)
for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:SchroedingerHermite1D, :_schroedinger_hermite_1d, true, 1 ),             
                 (:SchroedingerHermite2D, :_schroedinger_hermite_2d, true, 2 ),             
                 (:SchroedingerHermite3D, :_schroedinger_hermite_3d, true, 3 ),             
                 (:SchroedingerHermiteReal1D, :_schroedinger_hermite_real_1d, false, 1 ), 
                 (:SchroedingerHermiteReal2D, :_schroedinger_hermite_real_2d, false, 2 ),
                 (:SchroedingerHermiteReal3D, :_schroedinger_hermite_real_3d, false, 3 ),
                )
println("    ", METHOD) 
if T == :Float128
    SUF = Symbol(SUF, "_wf128")
end

if DIM==1
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, omega_x::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       potential_t::Function=none_2D,
                       potential_t_derivative::Function=none_2D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_1D
        with_potential_t = potential_t!=none_2D
        with_potential_t_derivative = potential_t_derivative!=none_2D
        V_c = cfunction(potential, ($T), (($T),))
        V_t_c = cfunction(potential_t, ($T), (($T),($T)))
        V_t_derivative_c = cfunction(potential_t_derivative, ($T), (($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Void}, 
                   (Int32, ($T), ($T), ($T), Ptr{Void}, Bool, Ptr{Void}, Bool, Ptr{Void}, Bool, ($T)), 
                   nx, omega_x, 
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t, V_t_derivative_c, with_potential_t_derivative,
                   cubic_coupling) 
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string(PRE,"finalize",SUF))), Void, (Ptr{Void},), x.m) )
       m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, omega_x::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       potential_t::Function=none_2D,
                       potential_t_derivative::Function=none_2D,
                       cubic_coupling::Real=0.0)
        ($METHOD)(($T), nx, omega_x,
                       hbar=hbar, mass=mass, potential=potential, 
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling)
    end      
end # eval
end #if
elseif DIM==2
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, omega_x::Real,
                                      ny::Integer, omega_y::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       potential_t::Function=none_3D,
                       potential_t_derivative::Function=none_3D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_2D
        with_potential_t = potential_t!=none_3D
        with_potential_t_derivative = potential_t_derivative!=none_3D
        V_c = cfunction(potential, ($T), (($T),($T)))
        V_t_c = cfunction(potential_t, ($T), (($T),($T),($T)))
        V_t_derivative_c = cfunction(potential_t_derivative, ($T), (($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Void}, 
                   (Int32, ($T), Int32, ($T), 
                   ($T), ($T), Ptr{Void}, Bool, Ptr{Void}, Bool, Ptr{Void}, Bool, ($T)), 
                   nx, omega_x, ny, omega_y, 
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t, V_t_derivative_c, with_potential_t_derivative,
                   cubic_coupling) 
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string(PRE,"finalize",SUF))), Void, (Ptr{Void},), x.m) )
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, omega_x::Real, 
                       ny::Integer, omega_y::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       potential_t::Function=none_3D,
                       potential_t_derivative::Function=none_3D,
                       cubic_coupling::Real=0.0)
        ($METHOD)(($T), nx, omega_x,  ny, omega_y,
                       hbar=hbar, mass=mass, potential=potential, 
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling)
    end      
end # eval
end #if
elseif DIM==3
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, omega_x::Real,
                                      ny::Integer, omega_y::Real, 
                                      nz::Integer, omega_z::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       potential_t_derivative::Function=none_4D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_3D
        with_potential_t = potential_t!=none_4D
        with_potential_t_derivative = potential_t_derivative!=none_4D
        V_c = cfunction(potential, ($T), (($T),($T),($T)))
        V_t_c = cfunction(potential_t, ($T), (($T),($T),($T),($T)))
        V_t_derivative_c = cfunction(potential_t_derivative, ($T), (($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Void}, 
                   (Int32, ($T), Int32, ($T), Int32, ($T), 
                   ($T), ($T), Ptr{Void}, Bool, Ptr{Void}, Bool, Ptr{Void}, Bool, ($T)), 
                   nx, omega_x, ny, omega_y, nz, omega_z,  
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t, V_t_derivative_c, with_potential_t_derivative,
                   cubic_coupling) 
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string(PRE,"finalize",SUF))), Void, (Ptr{Void},), x.m) )
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, omega_x::Real, 
                       ny::Integer, omega_y::Real,
                       nz::Integer, omega_z::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       potential_t_derivative::Function=none_4D,
                       cubic_coupling::Real=0.0)
        ($METHOD)(($T), nx, omega_x, ny, omega_y, nz, omega_z,
                       hbar=hbar, mass=mass, potential=potential,
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling)
    end      
end # eval
end #if

end # if

@eval begin
    function get_nx(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nx",SUF))), Int32,
             (Ptr{Void}, ), m.m )
    end

    function get_omega_x(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_omega_x",SUF))), ($T),
             (Ptr{Void}, ), m.m )
    end
end # eval    

if DIM>=2        
    @eval begin
    
        function get_ny(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ny",SUF))), Int32,
                 (Ptr{Void}, ), m.m )
        end

        function get_omega_y(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_omega_y",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end
    end # eval    
end

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

if DIM==1
    @eval begin

        function get_weights(m::($METHOD){$T})
           dim =Array(Int32, 1)
           np = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
           n = pointer_to_array(np, dim[1], false)     
           copy(n)
        end

        function get_transformation_matrices(m::($METHOD){$T}, unsafe_access::Bool=false)
           dims =Array(Int32, 2)
           Hp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_H",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
           H = pointer_to_array(Hp1, (dims[1], dims[2]), false)  
           if unsafe_access
               return H
           else
               return copy(H)
           end
        end

    end  
elseif DIM==2
    @eval begin

        function get_weights(m::($METHOD){$T})
           dim =Array(Int32, 1)
           np1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
           n1 = pointer_to_array(n1p, dim[1], false)     
           np2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
           n2 = pointer_to_array(n2p, dim[1], false)     
           copy(n1), copy(n2)
        end

        function get_transformation_matrices(m::($METHOD){$T}, unsafe_access::Bool=false)
           dims =Array(Int32, 2)
           Hp1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_H",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
           H1 = pointer_to_array(Hp1, (dims[1], dims[2]), false)  
           Hp2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_H",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 2 )
           H2 = pointer_to_array(Hp2, (dims[1], dims[2]), false)  
           if unsafe_access
               return H1, H2
           else
               return copy(H1), copy(H2)
           end
        end
       
    end  
elseif DIM==2
    @eval begin

        function get_weights(m::($METHOD){$T})
           dim =Array(Int32, 1)
           np1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
           n1 = pointer_to_array(n1p, dim[1], false)     
           np2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
           n2 = pointer_to_array(n2p, dim[1], false)     
           np3 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_weights",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
           n3 = pointer_to_array(n3p, dim[1], false)     
           copy(n1), copy(n2), copy(n3)
        end

        function get_transformation_matrices(m::($METHOD){$T}, unsafe_access::Bool=false)
           dims =Array(Int32, 2)
           Hp1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_H",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
           H1 = pointer_to_array(Hp1, (dims[1], dims[2]), false)  
           Hp2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_H",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 2 )
           H2 = pointer_to_array(Hp2, (dims[1], dims[2]), false)  
           Hp3 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_H",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 3 )
           H3 = pointer_to_array(Hp3, (dims[1], dims[2]), false)  
           if unsafe_access
               return H1, H2, H3
           else
               return copy(H1), copy(H2), copy(H3)
           end
        end
        
    end #eval 
end


end # for

