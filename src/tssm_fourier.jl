# T , TSSM_HANDLE 

println("including tssm_fourier.jl for type ",T)
for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:Fourier1D, :_fourier_1d, true, 1 ),             
                 (:Fourier2D, :_fourier_2d, true, 2 ),             
                 (:Fourier3D, :_fourier_3d, true, 3 ),             
                 (:FourierReal1D, :_fourier_real_1d, false, 1 ), 
                 (:FourierReal2D, :_fourier_real_2d, false, 2 ),
                 (:FourierReal3D, :_fourier_real_3d, false, 3 ),
                )
println("    ", METHOD)     

if T == :Float128
    SUF = string(SUF, "_wf128")
end
if DIM==1
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, xmin::Real, xmax::Real; 
                       boundary_conditions::Integer=periodic)
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))),
            Ptr{Nothing}, (Int32, ($T), ($T), Int32), nx, xmin, xmax, boundary_conditions)
        m = ($METHOD){$T}(c)
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
       m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real; 
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax,  boundary_conditions=boundary_conditions)
    end      
end # eval
end #if
elseif DIM==2
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, xmin::Real, xmax::Real,
                                       ny::Integer, ymin::Real, ymax::Real; 
                       boundary_conditions::Integer=periodic)
        m = ($METHOD){$T}( ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), 
                        Ptr{Nothing}, (Int32, ($T), ($T), Int32, ($T), ($T), Int32), 
                        nx, xmin, xmax, ny, ymin, ymax, boundary_conditions))
        finalizer(x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string(PRE,"finalize",SUF))), Nothing, (Ptr{Nothing},), x.m), m)
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(nx::Integer, xmin::Real, xmax::Real, 
                       ny::Integer, ymin::Real, ymax::Real;
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax,  ny, ymin, ymax,
                  boundary_conditions=boundary_conditions)
    end      
end # eval
end #if
elseif DIM==3
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, xmin::Real, xmax::Real,
                                       ny::Integer, ymin::Real, ymax::Real, 
                                       nz::Integer, zmin::Real, zmax::Real; 
                       boundary_conditions::Integer=periodic)
        m = ($METHOD){$T}( ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new",SUF))), 
                        Ptr{Nothing}, (Int32, ($T), ($T), Int32, ($T), ($T), 
                        Int32, ($T), ($T), Int32), 
                        nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, 
                        boundary_conditions))
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
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,
                  boundary_conditions=boundary_conditions)
    end      
end # eval
end #if

end # if

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

