__precompile__()
module TSSM

push!(LOAD_PATH, dirname(@__FILE__))

import Base.copy!
import Base.scale!
import Base.norm

export WaveFunction, WaveFunction1D, WaveFunction2D, WaveFunction3D
export WaveFunctionComplex, WaveFunctionComplex1D, WaveFunctionComplex2D, WaveFunctionComplex3D
export WaveFunctionReal, WaveFunctionReal1D, WaveFunctionReal2D, WaveFunctionReal3D

export TimeSplittinSpectralMethod, TimeSplittinSpectralMethod1D, TimeSplittinSpectralMethod2D, TimeSplittinSpectralMethod3D
export TimeSplittinSpectralMethodComplex, TimeSplittinSpectralMethodComplex1D, TimeSplittinSpectralMethodComplex2D, TimeSplittinSpectralMethodComplex3D
export TimeSplittinSpectralMethodReal, TimeSplittinSpectralMethodReal1D, TimeSplittinSpectralMethodReal2D, TimeSplittinSpectralMethodReal3D

export dim

export periodic, dirichlet, neumann

export Fourier1D, WfFourier1D 
export Fourier2D, WfFourier2D 
export Fourier3D, WfFourier3D 
export FourierReal1D, WfFourierReal1D 
export FourierReal2D, WfFourierReal2D 
export FourierReal3D, WfFourierReal3D 

export Schroedinger1D, WfSchroedinger1D 
export Schroedinger2D, WfSchroedinger2D 
export Schroedinger3D, WfSchroedinger3D 
export SchroedingerReal1D, WfSchroedingerReal1D 
export SchroedingerReal2D, WfSchroedingerReal2D 
export SchroedingerReal3D, WfSchroedingerReal3D 

export SchroedingerHermite1D, WfSchroedingerHermite1D 
export SchroedingerHermite2D, WfSchroedingerHermite2D 
export SchroedingerHermite3D, WfSchroedingerHermite3D 
export SchroedingerHermiteReal1D, WfSchroedingerHermiteReal1D 
export SchroedingerHermiteReal2D, WfSchroedingerHermiteReal2D 
export SchroedingerHermiteReal3D, WfSchroedingerHermiteReal3D 

export GeneralizedLaguerre2D, WfGeneralizedLaguerre2D
export GeneralizedLaguerreHermite3D, WfGeneralizedLaguerreHermite3D
export GeneralizedLaguerreReal2D, WfGeneralizedLaguerreReal2D
export GeneralizedLaguerreHermiteReal3D, WfGeneralizedLaguerreHermiteReal3D

export gauss, radau, lobatto

export FourierBessel2D, WfFourierBessel2D 
export FourierBesselReal2D, WfFourierBesselReal2D
export BesselRotSym1D, WfBesselRotSym1D 
export BesselRotSymReal1D, WfBesselRotSymReal1D 

export wave_function, clone
export is_real_space, is_frequency_space, to_real_space!, to_frequency_space!
export propagate_A!, propagate_B!, add_apply_A!
export norm, norm_in_frequency_space, normalize!, distance, scale!, axpy! 
export inner_product, eigen_function!, evaluate
export save, load!, get_data, set!, copy!
export get_eigenvalues, get_nodes, get_weights
export get_nx, get_ny, get_nz
export get_xmin, get_xmax
export get_ymin, get_ymax
export get_zmin, get_zmax
export get_omega_x,  get_omega_y,  get_omega_z
export get_H
export get_L, get_Omega 

export get_hbar, get_mass, get_cubic_coupling
export set_potential!, get_potential, load_potential!, save_potential
export imaginary_time_propagate_A!, imaginary_time_propagate_B!, add_apply_B!
export kinetic_energy, potential_energy, interaction_energy, observable
export get_energy_expectation_deviation, get_realspace_observables
export selfconsistent_nonlinear_step!

## abstract types ############################################################################

abstract TimeSplittingSpectralMethod{T<:AbstractFloat}
abstract TimeSplittingSpectralMethodReal{T<:AbstractFloat} <: TimeSplittingSpectralMethod{T}
abstract TimeSplittingSpectralMethodReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal{T}
abstract TimeSplittingSpectralMethodReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal{T}
abstract TimeSplittingSpectralMethodReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal{T}
abstract TimeSplittingSpectralMethodComplex{T<:AbstractFloat} <: TimeSplittingSpectralMethod{T}
abstract TimeSplittingSpectralMethodComplex1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex{T}
abstract TimeSplittingSpectralMethodComplex2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex{T}
abstract TimeSplittingSpectralMethodComplex3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex{T}
TimeSplittingSpectralMethod1D = Union{TimeSplittingSpectralMethodReal1D, TimeSplittingSpectralMethodComplex1D}
TimeSplittingSpectralMethod2D = Union{TimeSplittingSpectralMethodReal2D, TimeSplittingSpectralMethodComplex2D}
TimeSplittingSpectralMethod3D = Union{TimeSplittingSpectralMethodReal3D, TimeSplittingSpectralMethodComplex3D}

abstract WaveFunction{T<:AbstractFloat}
abstract WaveFunctionReal{T<:AbstractFloat} <: WaveFunction{T}
abstract WaveFunctionReal1D{T<:AbstractFloat} <: WaveFunctionReal{T}
abstract WaveFunctionReal2D{T<:AbstractFloat} <: WaveFunctionReal{T}
abstract WaveFunctionReal3D{T<:AbstractFloat} <: WaveFunctionReal{T}
abstract WaveFunctionComplex{T<:AbstractFloat} <: WaveFunction{T}
abstract WaveFunctionComplex1D{T<:AbstractFloat} <: WaveFunctionComplex{T}
abstract WaveFunctionComplex2D{T<:AbstractFloat} <: WaveFunctionComplex{T}
abstract WaveFunctionComplex3D{T<:AbstractFloat} <: WaveFunctionComplex{T}
WaveFunction1D = Union{WaveFunctionReal1D, WaveFunctionComplex1D}
WaveFunction2D = Union{WaveFunctionReal2D, WaveFunctionComplex2D}
WaveFunction3D = Union{WaveFunctionReal3D, WaveFunctionComplex3D}

dim(psi::WaveFunction1D) = 1
dim(psi::WaveFunction2D) = 2
dim(psi::WaveFunction3D) = 3
dim(psi::TimeSplittingSpectralMethod1D) = 1
dim(psi::TimeSplittingSpectralMethod2D) = 2
dim(psi::TimeSplittingSpectralMethod3D) = 3


## Fourier types ############################################################################

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

## FourierBessel types ###########################################################################

type FourierBessel2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Void}
end 

type FourierBesselReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Void}
end 

type BesselRotSym1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Void}
end 

type BesselRotSymReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Void}
end 

type WfFourierBessel2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Void}
    m::FourierBessel2D{T}
end

type WfFourierBesselReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Void}
    m::FourierBesselReal2D{T}
end

type WfBesselRotSym1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Void}
    m::BesselRotSym1D{T}
end

type WfBesselRotSymReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Void}
    m::BesselRotSymReal1D{T}
end

## Schroedinger types ############################################################################

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

## SchroedingerHermite types ############################################################################

type SchroedingerHermite1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Void}
end 

type SchroedingerHermite2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Void}
end 

type SchroedingerHermite3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Void}
end 

type SchroedingerHermiteReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Void}
end 

type SchroedingerHermiteReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Void}
end 

type SchroedingerHermiteReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Void}
end 

type WfSchroedingerHermite1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Void}
    m::SchroedingerHermite1D{T}
end

type WfSchroedingerHermite2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Void}
    m::SchroedingerHermite2D{T}
end

type WfSchroedingerHermite3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Void}
    m::SchroedingerHermite3D{T}
end

type WfSchroedingerHermiteReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Void}
    m::SchroedingerHermiteReal1D{T}
end

type WfSchroedingerHermiteReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Void}
    m::SchroedingerHermiteReal2D{T}
end

type WfSchroedingerHermiteReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Void}
    m::SchroedingerHermiteReal3D{T}
end


const periodic = 0
const dirichlet = 1
const neumann = 2

const gauss = 1
const radau = 2
const lobatto = 3


none_1D(x)=zero(x)
none_2D(x,y)=zero(x)
none_3D(x,y,z)=zero(x)

const libtssm = joinpath(dirname(@__FILE__),  "..", "deps", "usr", "lib",
                     string("libtssm.", Libdl.dlext))
const libtssm_debug = joinpath(dirname(@__FILE__),  "..", "deps", "usr", "lib",
                     string("libtssm_debug.", Libdl.dlext))
const libtssmq = joinpath(dirname(@__FILE__),  "..", "deps", "usr", "lib",
                     string("libtssmq.", Libdl.dlext))

__use_Float128 = false
try
   h=Libdl.dlopen(libtssmq);
   using Quadmath
   __use_Float128 = true
   Libdl.dlclose(h);
end   
const use_Float128 = __use_Float128                     

function __init__()
    if searchindex(readall(`uname -a`), "juliabox")>0
        # In JuliaBox only 8 out of 16 cores are available.
        ENV["OMP_NUM_THREADS"] = "8"
    end

    global tssm_handle
    if "TSSM_DEBUG" in keys(ENV) && ENV["TSSM_DEBUG"]=="1"
        warn("using TSSM debug version")
        tssm_handle = Libdl.dlopen(libtssm_debug) 
    else
        tssm_handle = Libdl.dlopen(libtssm) 
    end    
    ccall( Libdl.dlsym(tssm_handle, "tssm_initialize"), Void, ())
    ccall( Libdl.dlsym(tssm_handle, "tssm_fourier_initialize"), Void, ())

    global tssmq_handle 
    if use_Float128
        tssmq_handle = Libdl.dlopen(libtssmq) 
        ccall( Libdl.dlsym(tssmq_handle, "tssmq_initialize"), Void, ())
        ccall( Libdl.dlsym(tssmq_handle, "tssmq_fourier_initialize"), Void, ())
    end    
    set_fftw_planning_rigor(FFTW.ESTIMATE)
end


function set_fftw_planning_rigor(flag::Integer=FFTW.ESTIMATE)
   if !(flag in [ FFTW.ESTIMATE, FFTW.PATIENT, FFTW.MEASURE])
       error("wrong planning rigor flag")
   end
   ccall( Libdl.dlsym(tssm_handle, "tssm_set_fftw_planning_rigor"), Void, (Int32,), flag )
    if use_Float128
       ccall( Libdl.dlsym(tssmq_handle, "tssmq_set_fftw_planning_rigor"), Void, (Int32,), flag )
   end
end


T = :Float64
TSSM_HANDLE = :tssm_handle
PRE = :tssm_
include("tssm_fourier.jl")
include("tssm_fourier_bessel.jl")
include("tssm_schroedinger.jl")
include("tssm_schroedinger_hermite.jl")
#include("tssm_generalized_laguerre.jl")
include("tssm_common.jl")
include("tssm_schroedinger_common.jl")


if use_Float128
    T = :Float128
    TSSM_HANDLE = :tssmq_handle
    PRE = :tssmq_
    include("tssm_fourier.jl")
    include("tssm_fourier_bessel.jl")
    include("tssm_schroedinger.jl")
    include("tssm_schroedinger_hermite.jl")
    #include("tssm_generalized_laguerre.jl")
    include("tssm_common.jl")
    include("tssm_schroedinger_common.jl")
end    



end # TSSM
