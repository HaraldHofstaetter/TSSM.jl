# T , TSSM_HANDLE 

## Method: Schroedinger ############################################################################

type Schroedinger1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Void}
end 

type Schroedinger2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Void}
end 

type Schroedinger3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Void}
end 

type SchroedingerReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Void}
end 

type SchroedingerReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Void}
end 

type SchroedingerReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Void}
end 

type WfSchroedinger1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Void}
    m::Schroedinger1D{T}
end

type WfSchroedinger2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Void}
    m::Schroedinger2D{T}
end

type WfSchroedinger3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Void}
    m::Schroedinger3D{T}
end

type WfSchroedingerReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Void}
    m::SchroedingerReal1D{T}
end

type WfSchroedingerReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Void}
    m::SchroedingerReal2D{T}
end

type WfSchroedingerReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Void}
    m::SchroedingerReal3D{T}
end

none_1D(x)=0.0
none_2D(x,y)=0.0
none_3D(x,y,z)=0.0


for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:Schroedinger1D, :_schroedinger_1d, true, 1 ),             
                 (:Schroedinger2D, :_schroedinger_2d, true, 2 ),             
                 (:Schroedinger3D, :_schroedinger_3d, true, 3 ),             
                 (:SchroedingerReal1D, :_schroedinger_real_1d, false, 1 ), 
                 (:SchroedingerReal2D, :_schroedinger_real_2d, false, 2 ),
                 (:SchroedingerReal3D, :_schroedinger_real_3d, false, 3 ),
                )
if DIM==1
@eval begin
    function ($METHOD)(T::Type{$T}, nx::Integer, xmin::Real, xmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ccall( Libdl.dlsym(($TSSM_HANDLE), "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_1D
        println("*** ", ($METHOD), " hbar=", hbar)
        println("*** ", ($METHOD), " mass=", mass)
        println("*** ", ($METHOD), " potential=", potential)
        println("*** ", ($METHOD), " cubic_coupling=", cubic_coupling)
        println("*** ", ($METHOD), " none_1D=", none_1D)
        V_c = cfunction(potential, ($T), (($T),))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new",SUF))), Ptr{Void}, 
                   (Int32, ($T), ($T), ($T), ($T), Ptr{Void}, Bool, ($T), Int32), 
                   nx, xmin, xmax, 
                   hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string("c_finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
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
                       hbar=hbar, mass=mass, potential=potential, cubic_coupling=cubic_coupling,
                       boundary_conditions=boundary_conditions)
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
        ccall( Libdl.dlsym(($TSSM_HANDLE), "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_2D
        V_c = cfunction(potential, ($T), (($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new",SUF))), Ptr{Void}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), Ptr{Void}, Bool, ($T), Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, 
                   hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                        $(string("c_finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
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
                       hbar=hbar, mass=mass, potential=potential, cubic_coupling=cubic_coupling,
                       boundary_conditions=boundary_conditions)
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
        ccall( Libdl.dlsym(($TSSM_HANDLE), "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_3D
        V_c = cfunction(potential, ($T), (($T),($T),($T)))
        c = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new",SUF))), Ptr{Void}, 
                   (Int32, ($T), ($T), Int32, ($T), ($T), Int32, ($T), ($T), 
                   ($T), ($T), Ptr{Void}, Bool, ($T), Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,  
                   hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) 
        m = ($METHOD){$T}(c)
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
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ($METHOD)(($T), nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,
                       hbar=hbar, mass=mass, potential=potential, cubic_coupling=cubic_coupling,
                       boundary_conditions=boundary_conditions)
    end      
end # eval
end #if

end # if

end # for

