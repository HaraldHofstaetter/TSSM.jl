# T , TSSM_HANDLE 

## Method: Fourier ############################################################################

type Fourier1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Void}
end 

type Fourier2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Void}
end 

type Fourier3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Void}
end 

type FourierReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Void}
end 

type FourierReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Void}
end 

type FourierReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Void}
end 

type WfFourier1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Void}
    m::Fourier1D{T}
end

type WfFourier2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Void}
    m::Fourier2D{T}
end

type WfFourier3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Void}
    m::Fourier3D{T}
end

type WfFourierReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Void}
    m::FourierReal1D{T}
end

type WfFourierReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Void}
    m::FourierReal2D{T}
end

type WfFourierReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Void}
    m::FourierReal3D{T}
end



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

