# T , TSSM_HANDLE 
println("including tssm_schroedinger.jl for type ", T)
for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:Schroedinger1D, :_schroedinger_1d, true, 1 ),             
                 (:Schroedinger2D, :_schroedinger_2d, true, 2 ),             
                 (:Schroedinger3D, :_schroedinger_3d, true, 3 ),             
                 (:SchroedingerReal1D, :_schroedinger_real_1d, false, 1 ), 
                 (:SchroedingerReal2D, :_schroedinger_real_2d, false, 2 ),
                 (:SchroedingerReal3D, :_schroedinger_real_3d, false, 3 ),
                )
println("    ", METHOD)     
if T == :Float128
    SUF = Symbol(SUF, "_wf128")
end

if COMPLEX_METHOD

if DIM==1
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, xmin::Real, xmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       potential_t::Function=none_2D,
                       potential_t_derivative::Function=none_2D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_1D
        with_potential_t = potential_t!=none_2D
        with_potential_t_derivative = potential_t_derivative!=none_2D
        V_c = cfunction_check_return_type(potential, ($T), (($T),))
        V_t_c = cfunction_check_return_type(potential_t, ($T), (($T),($T)))
        V_t_derivative_c = cfunction_check_return_type(potential_t_derivative, ($T), (($T),($T)))        
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Nothing}, 
                   (Int32, ($T), ($T), ($T), ($T), Ptr{Nothing}, Bool, Ptr{Nothing}, Bool, Ptr{Nothing}, Bool, ($T), Int32), 
                   nx, xmin, xmax, 
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t, V_t_derivative_c, with_potential_t_derivative,
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
       m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       potential_t::Function=none_2D,
                       potential_t_derivative::Function=none_2D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax,
                       hbar=hbar, mass=mass, potential=potential, 
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling, boundary_conditions=boundary_conditions)
    end      
end # eval
end #if
elseif DIM==2
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                      ny::Integer, ymin::Real, ymax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       potential_t::Function=none_3D,
                       potential_t_derivative::Function=none_3D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_2D
        with_potential_t = potential_t!=none_3D
        with_potential_t_derivative = potential_t_derivative!=none_3D
        V_c = cfunction_check_return_type(potential, ($T), (($T),($T)))
        V_t_c = cfunction_check_return_type(potential_t, ($T), (($T),($T),($T)))
        V_t_derivative_c = cfunction_check_return_type(potential_t_derivative, ($T), (($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Nothing}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), Ptr{Nothing}, Bool, Ptr{Nothing}, Bool, Ptr{Nothing}, Bool, ($T), Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, 
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t, V_t_derivative_c, with_potential_t_derivative,
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real, 
                       ny::Integer, ymin::Real, ymax::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       potential_t::Function=none_3D,
                       potential_t_derivative::Function=none_3D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax,  ny, ymin, ymax,
                       hbar=hbar, mass=mass, potential=potential, 
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling, boundary_conditions=boundary_conditions)
    end      
end # eval
end #if
elseif DIM==3
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                      ny::Integer, ymin::Real, ymax::Real, 
                                      nz::Integer, zmin::Real, zmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       potential_t_derivative::Function=none_4D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_3D
        with_potential_t = potential_t!=none_4D
        with_potential_t_derivative = potential_t_derivative!=none_4D
        V_c = cfunction_check_return_type(potential, ($T), (($T),($T),($T)))
        V_t_c = cfunction_check_return_type(potential_t, ($T), (($T),($T),($T),($T)))
        V_t_derivative_c = cfunction_check_return_type(potential_t_derivative, ($T), (($T),($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Nothing}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), Ptr{Nothing}, Bool,  Ptr{Nothing}, Bool, Ptr{Nothing}, Bool, ($T), Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,  
                   hbar, mass, V_c, with_potential, 
                   V_t_c, with_potential_t, V_t_derivative_c, with_potential_t_derivative,
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real, 
                       ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       potential_t::Function=none_4D,
                       potential_t_derivative::Function=none_4D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,
                       hbar=hbar, mass=mass, potential=potential,
                       potential_t=potential_t, potential_t_derivative=potential_t_derivative,
                       cubic_coupling=cubic_coupling, boundary_conditions=boundary_conditions)
    end      
end # eval
end #if

end # if

else # !COMPLEX_METHOD

if DIM==1
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, xmin::Real, xmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_1D
        V_c = cfunction_check_return_type(potential, ($T), (($T),))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Nothing}, 
                   (Int32, ($T), ($T), ($T), ($T), Ptr{Nothing}, Bool, ($T), Int32), 
                   nx, xmin, xmax, 
                   hbar, mass, V_c, with_potential,
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
       m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax,
                       hbar=hbar, mass=mass, potential=potential,  
                       cubic_coupling=cubic_coupling, boundary_conditions=boundary_conditions)
    end      
end # eval
end #if
elseif DIM==2
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                      ny::Integer, ymin::Real, ymax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_2D
        V_c = cfunction_check_return_type(potential, ($T), (($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Nothing}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), Ptr{Nothing}, Bool, ($T), Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, 
                   hbar, mass, V_c, with_potential, 
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real, 
                       ny::Integer, ymin::Real, ymax::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax,  ny, ymin, ymax,
                       hbar=hbar, mass=mass, potential=potential,
                       cubic_coupling=cubic_coupling, boundary_conditions=boundary_conditions)
    end      
end # eval
end #if
elseif DIM==3
@eval begin
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                      ny::Integer, ymin::Real, ymax::Real, 
                                      nz::Integer, zmin::Real, zmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        with_potential = potential!=none_3D
        V_c = cfunction_check_return_type(potential, ($T), (($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), Ptr{Nothing}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), Ptr{Nothing}, Bool, ($T), Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,  
                   hbar, mass, V_c, with_potential,
                   cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real, 
                       ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,
                       hbar=hbar, mass=mass, potential=potential, 
                       cubic_coupling=cubic_coupling, boundary_conditions=boundary_conditions)
    end      
end # eval
end #if

end # if

end # if !COMPLEX_METHOD

@eval begin
    function get_nx(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nx",SUF))), Int32,
             (Ptr{Nothing}, ), m.m )
    end

    function get_xmin(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_xmin",SUF))), ($T),
             (Ptr{Nothing}, ), m.m )
    end

    function get_xmax(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_xmax",SUF))), $T,
             (Ptr{Nothing}, ), m.m )
    end
end # eval    

if DIM>=2        
    @eval begin
    
        function get_ny(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ny",SUF))), Int32,
                 (Ptr{Nothing}, ), m.m )
        end

        function get_ymin(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ymin",SUF))), ($T),
                 (Ptr{Nothing}, ), m.m )
        end

        function get_ymax(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_ymax",SUF))), ($T),
                 (Ptr{Nothing}, ), m.m )
        end
    end # eval    
end

if DIM>=3 
    @eval begin
        function get_nz(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nz",SUF))), Int32,
                 (Ptr{Nothing}, ), m.m )
        end

        function get_zmin(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_zmin",SUF))), ($T),
                 (Ptr{Nothing}, ), m.m )
        end

        function get_zmax(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_zmax",SUF))), ($T),
                 (Ptr{Nothing}, ), m.m )
        end
    end # eval    
end 


end # for

