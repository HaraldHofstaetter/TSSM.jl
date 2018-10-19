__precompile__()
module TSSM

push!(LOAD_PATH, dirname(@__FILE__))

import Base.copy!

export WaveFunction, WaveFunction1D, WaveFunction2D, WaveFunction3D
export WaveFunctionComplex, WaveFunctionComplex1D, WaveFunctionComplex2D, WaveFunctionComplex3D
export WaveFunctionReal, WaveFunctionReal1D, WaveFunctionReal2D, WaveFunctionReal3D

export TimeSplittinSpectralMethod, TimeSplittinSpectralMethod1D, TimeSplittinSpectralMethod2D, TimeSplittinSpectralMethod3D
export TimeSplittinSpectralMethodComplex, TimeSplittinSpectralMethodComplex1D, TimeSplittinSpectralMethodComplex2D, TimeSplittinSpectralMethodComplex3D
export TimeSplittinSpectralMethodReal, TimeSplittinSpectralMethodReal1D, TimeSplittinSpectralMethodReal2D, TimeSplittinSpectralMethodReal3D

export dim, is_real

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

export SchroedingerRotating2D, WfSchroedingerRotating2D 
export SchroedingerRotating3D, WfSchroedingerRotating3D 
export SchroedingerRotatingReal2D, WfSchroedingerRotatingReal2D 
export SchroedingerRotatingReal3D, WfSchroedingerRotatingReal3D 

export SchroedingerHermite1D, WfSchroedingerHermite1D 
export SchroedingerHermite2D, WfSchroedingerHermite2D 
export SchroedingerHermite3D, WfSchroedingerHermite3D 
export SchroedingerHermiteReal1D, WfSchroedingerHermiteReal1D 
export SchroedingerHermiteReal2D, WfSchroedingerHermiteReal2D 
export SchroedingerHermiteReal3D, WfSchroedingerHermiteReal3D 

export SchroedingerGeneralizedLaguerre2D, WfSchroedingerGeneralizedLaguerre2D
export SchroedingerGeneralizedLaguerreHermite3D, WfSchroedingerGeneralizedLaguerreHermite3D
export SchroedingerGeneralizedLaguerreReal2D, WfSchroedingerGeneralizedLaguerreReal2D
export SchroedingerGeneralizedLaguerreHermiteReal3D, WfSchroedingerGeneralizedLaguerreHermiteReal3D

export gauss, radau, lobatto

export FourierBessel2D, WfFourierBessel2D 
export FourierBesselReal2D, WfFourierBesselReal2D
export BesselRotSym1D, WfBesselRotSym1D 
export BesselRotSymReal1D, WfBesselRotSymReal1D 

export wave_function, clone
export is_real_space, is_frequency_space, to_real_space!, to_frequency_space!
export is_real_space_x, is_frequency_space_x, to_real_space_x!, to_frequency_space_x!
export is_real_space_y, is_frequency_space_y, to_real_space_y!, to_frequency_space_y!
export is_real_space_z, is_frequency_space_z, to_real_space_z!, to_frequency_space_z!
export set_time!, get_time, propagate_time!
export set_propagate_time_together_with_A!, get_propagate_time_together_with_A
export propagate_A!, propagate_B!, propagate_C!, add_apply_A!, add_apply_B!, add_apply_C!
export add_phi_A!
export propagate_A_derivative!, propagate_B_derivative!, propagate_C_derivative!
export norm, norm_in_frequency_space, normalize!, distance, scale!, axpy! 
export inner_product, eigen_function!, evaluate
export save, load!, get_data, set!, copy!
export get_eigenvalues, get_nodes, get_weights, get_transformation_matrices
export get_nx, get_ny, get_nz
export get_xmin, get_xmax
export get_ymin, get_ymax
export get_zmin, get_zmax
export get_omega_x,  get_omega_y,  get_omega_z
export get_omega_r, get_Omega 
export get_nr, get_nfr, get_ntheta

export get_hbar, get_mass, get_cubic_coupling, set_cubic_coupling
export set_potential!, get_potential, load_potential!, save_potential
export set_potential_t!, get_potential_t
export set_potential_t_derivative!, get_potential_t_derivative
export imaginary_time_propagate_A!, imaginary_time_propagate_B!
export kinetic_energy, potential_energy, interaction_energy, observable
export get_energy_expectation_deviation, get_realspace_observables
export kinetic_matrix_element, potential_matrix_element
export selfconsistent_nonlinear_step!

## abstract types ############################################################################
abstract type TimeSplittingSpectralMethod{T<:AbstractFloat} end
abstract type TimeSplittingSpectralMethodReal{T<:AbstractFloat}<:TimeSplittingSpectralMethod{T} end
abstract type TimeSplittingSpectralMethodReal1D{T<:AbstractFloat}<:TimeSplittingSpectralMethodReal{T} end
abstract type TimeSplittingSpectralMethodReal2D{T<:AbstractFloat}<:TimeSplittingSpectralMethodReal{T} end
abstract type TimeSplittingSpectralMethodReal3D{T<:AbstractFloat}<:TimeSplittingSpectralMethodReal{T} end
abstract type TimeSplittingSpectralMethodComplex{T<:AbstractFloat}<:TimeSplittingSpectralMethod{T} end
abstract type TimeSplittingSpectralMethodComplex1D{T<:AbstractFloat}<:TimeSplittingSpectralMethodComplex{T} end
abstract type TimeSplittingSpectralMethodComplex2D{T<:AbstractFloat}<:TimeSplittingSpectralMethodComplex{T} end
abstract type TimeSplittingSpectralMethodComplex3D{T<:AbstractFloat}<:TimeSplittingSpectralMethodComplex{T} end
TimeSplittingSpectralMethod1D = Union{TimeSplittingSpectralMethodReal1D, TimeSplittingSpectralMethodComplex1D}
TimeSplittingSpectralMethod2D = Union{TimeSplittingSpectralMethodReal2D, TimeSplittingSpectralMethodComplex2D}
TimeSplittingSpectralMethod3D = Union{TimeSplittingSpectralMethodReal3D, TimeSplittingSpectralMethodComplex3D}

abstract type WaveFunction{T<:AbstractFloat} end
abstract type WaveFunctionReal{T<:AbstractFloat}<:WaveFunction{T} end
abstract type WaveFunctionReal1D{T<:AbstractFloat}<:WaveFunctionReal{T} end
abstract type WaveFunctionReal2D{T<:AbstractFloat}<:WaveFunctionReal{T} end
abstract type WaveFunctionReal3D{T<:AbstractFloat}<:WaveFunctionReal{T} end
abstract type WaveFunctionComplex{T<:AbstractFloat}<:WaveFunction{T} end
abstract type WaveFunctionComplex1D{T<:AbstractFloat}<:WaveFunctionComplex{T} end
abstract type WaveFunctionComplex2D{T<:AbstractFloat}<:WaveFunctionComplex{T} end
abstract type WaveFunctionComplex3D{T<:AbstractFloat}<:WaveFunctionComplex{T} end
WaveFunction1D = Union{WaveFunctionReal1D, WaveFunctionComplex1D}
WaveFunction2D = Union{WaveFunctionReal2D, WaveFunctionComplex2D}
WaveFunction3D = Union{WaveFunctionReal3D, WaveFunctionComplex3D}

dim(psi::WaveFunction1D) = 1
dim(psi::WaveFunction2D) = 2
dim(psi::WaveFunction3D) = 3
dim(psi::TimeSplittingSpectralMethod1D) = 1
dim(psi::TimeSplittingSpectralMethod2D) = 2
dim(psi::TimeSplittingSpectralMethod3D) = 3


is_real(psi::WaveFunctionReal) = true
is_real(psi::WaveFunctionComplex) = false 
is_real(psi::TimeSplittingSpectralMethodReal) = true
is_real(psi::TimeSplittingSpectralMethodComplex) = false


# Dummy propagators, to be overloaded
propagate_B!(psi::WaveFunction, dt::Number) = psi
propagate_C!(psi::WaveFunction, dt::Number) = psi


## Fourier types ############################################################################

mutable struct Fourier1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Nothing}
end 

mutable struct Fourier2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Nothing}
end 

mutable struct Fourier3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Nothing}
end 

mutable struct FourierReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Nothing}
end 

mutable struct FourierReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Nothing}
end 

mutable struct FourierReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Nothing}
end 

mutable struct WfFourier1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Nothing}
    m::Fourier1D{T}
end

mutable struct WfFourier2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Nothing}
    m::Fourier2D{T}
end

mutable struct WfFourier3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Nothing}
    m::Fourier3D{T}
end

mutable struct WfFourierReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Nothing}
    m::FourierReal1D{T}
end

mutable struct WfFourierReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Nothing}
    m::FourierReal2D{T}
end

mutable struct WfFourierReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Nothing}
    m::FourierReal3D{T}
end

## FourierBessel types ###########################################################################

mutable struct FourierBessel2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Nothing}
end 

mutable struct FourierBesselReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Nothing}
end 

mutable struct BesselRotSym1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Nothing}
end 

mutable struct BesselRotSymReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Nothing}
end 

mutable struct WfFourierBessel2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Nothing}
    m::FourierBessel2D{T}
end

mutable struct WfFourierBesselReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Nothing}
    m::FourierBesselReal2D{T}
end

mutable struct WfBesselRotSym1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Nothing}
    m::BesselRotSym1D{T}
end

mutable struct WfBesselRotSymReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Nothing}
    m::BesselRotSymReal1D{T}
end

## Schroedinger types ############################################################################

mutable struct Schroedinger1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Nothing}
end 

mutable struct Schroedinger2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Nothing}
end 

mutable struct Schroedinger3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Nothing}
end 

mutable struct WfSchroedinger1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Nothing}
    m::Schroedinger1D{T}
end

mutable struct WfSchroedinger2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Nothing}
    m::Schroedinger2D{T}
end

mutable struct WfSchroedinger3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Nothing}
    m::Schroedinger3D{T}
end

mutable struct WfSchroedingerReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Nothing}
    m::SchroedingerReal1D{T}
end

mutable struct WfSchroedingerReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Nothing}
    m::SchroedingerReal2D{T}
end

mutable struct WfSchroedingerReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Nothing}
    m::SchroedingerReal3D{T}
end

## SchroedingerRotating types ############################################################################

mutable struct SchroedingerRotating2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerRotating3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Nothing}
end 

mutable struct WfSchroedingerRotating2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Nothing}
    m::SchroedingerRotating2D{T}
end

mutable struct WfSchroedingerRotating3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Nothing}
    m::SchroedingerRotating3D{T}
end

mutable struct SchroedingerRotatingReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerRotatingReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Nothing}
end 

mutable struct WfSchroedingerRotatingReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Nothing}
    m::SchroedingerRotatingReal2D{T}
end

mutable struct WfSchroedingerRotatingReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Nothing}
    m::SchroedingerRotatingReal3D{T}
end


## SchroedingerHermite types ############################################################################

mutable struct SchroedingerHermite1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex1D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerHermite2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerHermite3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerHermiteReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal1D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerHermiteReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerHermiteReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Nothing}
end 

mutable struct WfSchroedingerHermite1D{T<:AbstractFloat} <: WaveFunctionComplex1D{T}
    p::Ptr{Nothing}
    m::SchroedingerHermite1D{T}
end

mutable struct WfSchroedingerHermite2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Nothing}
    m::SchroedingerHermite2D{T}
end

mutable struct WfSchroedingerHermite3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Nothing}
    m::SchroedingerHermite3D{T}
end

mutable struct WfSchroedingerHermiteReal1D{T<:AbstractFloat} <: WaveFunctionReal1D{T}
    p::Ptr{Nothing}
    m::SchroedingerHermiteReal1D{T}
end

mutable struct WfSchroedingerHermiteReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Nothing}
    m::SchroedingerHermiteReal2D{T}
end

mutable struct WfSchroedingerHermiteReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Nothing}
    m::SchroedingerHermiteReal3D{T}
end

## SchroedingerGeneralizedLaguerre(Hermite) types ########################################################


mutable struct SchroedingerGeneralizedLaguerre2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex2D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerGeneralizedLaguerreHermite3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex3D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerGeneralizedLaguerreReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal2D{T}
    m::Ptr{Nothing}
end 

mutable struct SchroedingerGeneralizedLaguerreHermiteReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal3D{T}
    m::Ptr{Nothing}
end 

mutable struct WfSchroedingerGeneralizedLaguerre2D{T<:AbstractFloat} <: WaveFunctionComplex2D{T}
    p::Ptr{Nothing}
    m::SchroedingerGeneralizedLaguerre2D{T}
end

mutable struct WfSchroedingerGeneralizedLaguerreHermite3D{T<:AbstractFloat} <: WaveFunctionComplex3D{T}
    p::Ptr{Nothing}
    m::SchroedingerGeneralizedLaguerreHermite3D{T}
end

mutable struct WfSchroedingerGeneralizedLaguerreReal2D{T<:AbstractFloat} <: WaveFunctionReal2D{T}
    p::Ptr{Nothing}
    m::SchroedingerGeneralizedLaguerreReal2D{T}
end

mutable struct WfSchroedingerGeneralizedLaguerreHermiteReal3D{T<:AbstractFloat} <: WaveFunctionReal3D{T}
    p::Ptr{Nothing}
    m::SchroedingerGeneralizedLaguerreHermiteReal3D{T}
end



const periodic = 0
const dirichlet = 1
const neumann = 2

const gauss = 1
const radau = 2
const lobatto = 3

const MEASURE         = UInt32(0)
const PATIENT         = UInt32(1 << 5)
const ESTIMATE        = UInt32(1 << 6)

none_1D(x)=zero(x)
none_2D(x,y)=zero(x)
none_3D(x,y,z)=zero(x)
none_4D(x,y,z,t)=zero(x)

using Libdl

const libtssm = joinpath(dirname(@__FILE__),  "..", "deps", "usr", "lib",
                     string("libtssm.", dlext))
const libtssm_debug = joinpath(dirname(@__FILE__),  "..", "deps", "usr", "lib",
                     string("libtssm_debug.", dlext))
const libtssmq = joinpath(dirname(@__FILE__),  "..", "deps", "usr", "lib",
                     string("libtssmq.", dlext))

# __use_Float128 = false
#try
#   h=dlopen(libtssmq);
#   using Quadmath
#   __use_Float128 = true
#   dlclose(h);
#finally   
#end   
#const use_Float128 = __use_Float128                     
const use_Float128 = false 

function __init__()
    global tssm_handle
    if "TSSM_DEBUG" in keys(ENV) && ENV["TSSM_DEBUG"]=="1"
        warn("using TSSM debug version")
        tssm_handle = dlopen(libtssm_debug) 
    else
        tssm_handle = dlopen(libtssm) 
    end    
    ccall( dlsym(tssm_handle, "tssm_initialize"), Nothing, ())
    ccall( dlsym(tssm_handle, "tssm_fourier_initialize"), Nothing, ())

    global tssmq_handle 
    if use_Float128
        tssmq_handle = dlopen(libtssmq) 
        ccall( dlsym(tssmq_handle, "tssmq_initialize"), Nothing, ())
        ccall( dlsym(tssmq_handle, "tssmq_fourier_initialize"), Nothing, ())
    end    
    set_fftw_planning_rigor(ESTIMATE)
end

function set_fftw_planning_rigor(flag::Integer=ESTIMATE)
   if !(flag in [ ESTIMATE, PATIENT, MEASURE])
       error("wrong planning rigor flag")
   end
   ccall( dlsym(tssm_handle, "tssm_set_fftw_planning_rigor"), Nothing, (Int32,), flag )
    if use_Float128
       ccall( dlsym(tssmq_handle, "tssmq_set_fftw_planning_rigor"), Nothing, (Int32,), flag )
   end
end

function cfunction_check_return_type(f, r, a)
    rt = Base.return_types(f, a)
    @assert(length(rt)==1 && rt[1]==r, "wrong return type of function")
    cfunction(f, r, a)
end


T = :Float64
TSSM_HANDLE = :tssm_handle
PRE = :tssm_
include("tssm_fourier.jl")
include("tssm_fourier_bessel.jl")
include("tssm_schroedinger.jl")
include("tssm_schroedinger_rotating.jl")
include("tssm_schroedinger_hermite.jl")
include("tssm_schroedinger_generalized_laguerre.jl")
include("tssm_common.jl")
include("tssm_schroedinger_common.jl")


if use_Float128
    T = :Float128
    TSSM_HANDLE = :tssmq_handle
    PRE = :tssmq_
    include("tssm_fourier.jl")
    include("tssm_fourier_bessel.jl")
    include("tssm_schroedinger.jl")
    include("tssm_schroedinger_rotating.jl")
    include("tssm_schroedinger_hermite.jl")
    include("tssm_schroedinger_generalized_laguerre.jl")
    include("tssm_common.jl")
    include("tssm_schroedinger_common.jl")
end    



end # TSSM
