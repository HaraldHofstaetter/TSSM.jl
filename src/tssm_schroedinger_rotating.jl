# T , TSSM_HANDLE 
println("including tssm_schroedinger_rotating.jl for type ", T)
for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:SchroedingerRotating2D, :_schroedinger_rotating_2d, true, 2 ),             
                 (:SchroedingerRotating3D, :_schroedinger_rotating_3d, true, 3 ),             
                 (:SchroedingerRotatingReal2D, :_schroedinger_rotating_real_2d, false, 2 ),             
                 (:SchroedingerRotatingReal3D, :_schroedinger_rotating_real_3d, false, 3 ),             
                )
println("    ", METHOD)     
if T == :Float128
    SUF = Symbol(SUF, "_wf128")
end
WF = Symbol(:Wf,METHOD)

if DIM==2
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                      ny::Integer, ymin::Real, ymax::Real; 
                       Omega::Real=0.0,
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       potential_t::Function=none_3D,
                       potential_t_derivative::Function=none_3D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_2D
        with_potential_t = potential_t!=none_3D
        with_potential_t_derivative = potential_t_derivative!=none_3D
        V_c = cfunction_check_return_type(potential, ($T), (($T),($T)))
        V_t_c = cfunction_check_return_type(potential_t, ($T), (($T),($T),($T)))
        V_t_derivative_c = cfunction_check_return_type(potential_t_derivative, ($T), (($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Void}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), ($T), Ptr{Void}, Bool, Ptr{Void}, Bool, Ptr{Void}, Bool, ($T)), 
                   nx, xmin, xmax, ny, ymin, ymax, 
                   Omega, hbar, mass, V_c, with_potential, 
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
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real, 
                       ny::Integer, ymin::Real, ymax::Real;
                       Omega::Real=0.0,
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       potential_t::Function=none_3D,
                       potential_t_derivative::Function=none_3D,
                       cubic_coupling::Real=0.0)
        ($METHOD)(($T), nx, xmin, xmax,  ny, ymin, ymax,
                       Omega=Omega, hbar=hbar, mass=mass, potential=potential, 
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling)
    end      
end # eval
end #if
elseif DIM==3
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                      ny::Integer, ymin::Real, ymax::Real, 
                                      nz::Integer, zmin::Real, zmax::Real; 
                       Omega::Real=0.0,
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       potential_t_derivative::Function=none_4D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_3D
        with_potential_t = potential_t!=none_4D
        with_potential_t_derivative = potential_t_derivative!=none_4D
        V_c = cfunction_check_return_type(potential, ($T), (($T),($T),($T)))
        V_t_c = cfunction_check_return_type(potential_t, ($T), (($T),($T),($T),($T)))
        V_t_derivative_c = cfunction_check_return_type(potential_t_derivative, ($T), (($T),($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Void}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), ($T), Ptr{Void}, Bool,  Ptr{Void}, Bool, Ptr{Void}, Bool, ($T)), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,  
                   Omega, hbar, mass, V_c, with_potential, 
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
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real, 
                       ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real;
                       Omega::Real=0.0,
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       potential_t_derivative::Function=none_4D,
                       cubic_coupling::Real=0.0)
        ($METHOD)(($T), nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,
                       Omega=Omega, hbar=hbar, mass=mass, potential=potential,
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling)
    end      
end # eval
end #if

end # if

@eval begin

    function get_Omega(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_Omega",SUF))), ($T),
             (Ptr{Void}, ), m.m )
    end

    function propagate_C!(psi::($WF){$T}, dt::Number)
        ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_C_wf",SUF))), Void,
                (Ptr{Void}, Complex{$T},), psi.p, dt)
    end

    function propagate_C_derivative!(this::($WF){$T}, other::($WF){$T},
                         dt::Number)
       if this.m ≠ other.m
           error("this and other must belong to the same method")
       end
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_C_derivative_wf",SUF))), Void,
                (Ptr{Void}, Ptr{Void}, Complex{$T}), 
                 this.p, other.p, dt)
    end

    function add_apply_C!(this::($WF){$T}, other::($WF){$T},
                         coefficient::Number=1.0)
       if this.m ≠ other.m
           error("this and other must belong to the same method")
       end
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"add_apply_C_wf",SUF))), Void,
                (Ptr{Void}, Ptr{Void}, Complex{$T}), 
                 this.p, other.p, coefficient)
    end

    function imaginary_time_propagate_C!(psi::($WF){$T}, dt::Number)
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"imaginary_time_propagate_C_wf",SUF))), Void,
                        (Ptr{Void}, Complex{($T)},), psi.p, dt)
    end
    
    function get_nx(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nx",SUF))), Int32,
             (Ptr{Void}, ), m.m )
    end

    function get_xmin(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_xmin",SUF))), ($T),
             (Ptr{Void}, ), m.m )
    end

    function get_xmax(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_xmax",SUF))), $T,
             (Ptr{Void}, ), m.m )
    end

    function is_real_space_x(psi::($WF){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_x_wf",SUF))), Int32,
           (Ptr{Void},), psi.p) == 1
    end

    function is_frequency_space_x(psi::($WF){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_x_wf",SUF))), Int32,
           (Ptr{Void},), psi.p) != 1
    end

    function to_real_space_x!(psi::($WF){$T})
         ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_real_space_x_wf",SUF))), Void,
                (Ptr{Void},), psi.p)
    end

    function to_frequency_space_x!(psi::($WF){$T})
         ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_frequency_space_x_wf",SUF))), Void,
                (Ptr{Void},), psi.p)
    end
end # eval    

if DIM>=2        
    @eval begin
    
        function get_ny(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ny",SUF))), Int32,
                 (Ptr{Void}, ), m.m )
        end

        function get_ymin(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ymin",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end

        function get_ymax(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ymax",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end
        function is_real_space_y(psi::($WF){$T})
        ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_y_wf",SUF))), Int32,
            (Ptr{Void},), psi.p) == 1
        end

        function is_frequency_space_y(psi::($WF){$T})
        ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_y_wf",SUF))), Int32,
            (Ptr{Void},), psi.p) != 1
        end

        function to_real_space_y!(psi::($WF){$T})
            ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_real_space_y_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end

        function to_frequency_space_y!(psi::($WF){$T})
            ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_frequency_space_y_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end
    end # eval    
end

if DIM>=3 
    @eval begin
        function get_nz(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nz",SUF))), Int32,
                 (Ptr{Void}, ), m.m )
        end

        function get_zmin(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_zmin",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end

        function get_zmax(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_zmax",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end

        function is_real_space_z(psi::($WF){$T})
        ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_z_wf",SUF))), Int32,
            (Ptr{Void},), psi.p) == 1
        end

        function is_frequency_space_z(psi::($WF){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_z_wf",SUF))), Int32,
            (Ptr{Void},), psi.p) != 1
        end
    
        function to_real_space_z!(psi::($WF){$T})
            ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_real_space_z_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end

        function to_frequency_space_z!(psi::($WF){$T})
            ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_frequency_space_z_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end
        
    end # eval    
end 

end # for

