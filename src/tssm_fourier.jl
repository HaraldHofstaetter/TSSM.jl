# T , TSSM_HANDLE 


for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:Fourier1D, :_fourier_1d, true, 1 ),             
                 (:Fourier2D, :_fourier_2d, true, 2 ),             
                 (:Fourier3D, :_fourier_3d, true, 3 ),             
                 (:FourierReal1D, :_fourier_real_1d, false, 1 ), 
                 (:FourierReal2D, :_fourier_real_2d, false, 2 ),
                 (:FourierReal3D, :_fourier_real_3d, false, 3 ),
                )
if DIM==1
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, xmin::Real, xmax::Real; 
                       boundary_conditions::Integer=periodic)
        ccall( Libdl.dlsym(($TSSM_HANDLE), "c_initialize_tssm_fourier"), Void, ())
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new",SUF))),
            Ptr{Void}, (Int32, ($T), ($T), Int32), nx, xmin, xmax, boundary_conditions)
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string("c_finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
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
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                       ny::Integer, ymin::Real, ymax::Real; 
                       boundary_conditions::Integer=periodic)
        ccall( Libdl.dlsym(($TSSM_HANDLE), "c_initialize_tssm_fourier"), Void, ())
        m = ($METHOD){$T}( ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new",SUF))), 
                        Ptr{Void}, (Int32, ($T), ($T), Int32, ($T), ($T), Int32), 
                        nx, xmin, xmax, ny, ymin, ymax, boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string("c_finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
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
    function ($METHOD)(T::Type{($T)}, nx::Integer, xmin::Real, xmax::Real,
                                       ny::Integer, ymin::Real, ymax::Real, 
                                       nz::Integer, zmin::Real, zmax::Real; 
                       boundary_conditions::Integer=periodic)
        ccall( Libdl.dlsym(($TSSM_HANDLE), "c_initialize_tssm_fourier"), Void, ())
        m = ($METHOD){$T}( ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new",SUF))), 
                        Ptr{Void}, (Int32, ($T), ($T), Int32, ($T), ($T), 
                        Int32, ($T), ($T), Int32), 
                        nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, 
                        boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string("c_finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
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

end # for

