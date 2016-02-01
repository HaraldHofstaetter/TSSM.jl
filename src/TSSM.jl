module TSSM

import Base.Libdl: dlsym, dlopen 

import Base.copy!
import Base.scale!
import Base.norm

export WaveFunction, WaveFunction1D, WaveFunction2D, WaveFunction3D, TimeSplittingSpectralMethod
export dim
#export tssm_handle

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

#export gauss, radau, lobatto

export FourierBessel2D, WfFourierBessel2D 
export FourierBesselReal2D, WfFourierBesselReal2D

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

function __init__()
    global tssm_handle
    if searchindex(readall(`uname -a`), "juliabox")>0
        # In JuliaBox only 8 out of 16 cores are available.
        ENV["OMP_NUM_THREADS"] = "8"
    end
    libtssm = string(joinpath(dirname(@__FILE__),  "..", "deps", "usr", "lib", "libtssm."), 
                   (@windows? :"dll" : ( @osx? "dylib" : :"so" )) )
    tssm_handle = Libdl.dlopen(libtssm);
    ccall( Libdl.dlsym(tssm_handle, "c_initialize_tssm"), Void, ())
    set_fftw_planning_rigor(FFTW.ESTIMATE)
end

function set_fftw_planning_rigor(flag::Integer=FFTW.ESTIMATE)
   if !(flag in [ FFTW.ESTIMATE, FFTW.PATIENT, FFTW.MEASURE])
       error("wrong planning rigor flag")
   end
   ccall( dlsym(tssm_handle, "c_set_fftw_planning_rigor"), Void, (Int32,), flag )
end

abstract TimeSplittingSpectralMethod

abstract WaveFunction
abstract WaveFunction1D <: WaveFunction
abstract WaveFunction2D <: WaveFunction
abstract WaveFunction3D <: WaveFunction

dim(psi::WaveFunction1D) = 1
dim(psi::WaveFunction2D) = 2
dim(psi::WaveFunction3D) = 3

## Boundary conditions ############################################################################

const periodic = 0
const dirichlet = 1
const neumann = 2

none_1D(x)=0.0
none_2D(x,y)=0.0
none_3D(x,y,z)=0.0

include("tssm_fourier.jl")
include("tssm_schroedinger.jl")
include("tssm_schroedinger_hermite.jl")
include("tssm_generalized_laguerre.jl")
include("tssm_fourier_bessel.jl")

function clone(psi::WaveFunction)
    wave_function(psi.m)
end

function eigen_function!(psi::WaveFunction1D, k)
    to_frequency_space!(psi)
    u=get_data(psi,true)
    u[:]=0.0
    u[k]=1
    psi
end

function eigen_function!(psi::WaveFunction2D, k, l)
    to_frequency_space!(psi)
    u=get_data(psi,true)
    u[:,:]=0.0
    u[k,l]=1
    psi
end

function eigen_function!(psi::WaveFunction3D, k, l, m)
    to_frequency_space!(psi)
    u=get_data(psi,true)
    u[:,:,:]=0.0
    u[k,l,m]=1
    psi
end

end # TSSM
