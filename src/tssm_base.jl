module tssm_base

export TimeSplittingSpectralMethod, TimeSplittingSpectralMethodReal, TimeSplittingSpectralMethodComplex
export TimeSplittingSpectralMethodReal1D, TimeSplittingSpectralMethodReal2D, TimeSplittingSpectralMethodReal3D
export TimeSplittingSpectralMethodComplex1D, TimeSplittingSpectralMethodComplex2D, TimeSplittingSpectralMethodComplex3D
export WaveFunction, WaveFuncionReal, WavefunctionComplex
export WaveFunctionReal1D, WaveFunctionReal2D, WaveFunctionReal3D
export WaveFunctionComplex1D, WaveFunctionComplex2D, WaveFunctionComplex3D

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

export dim

export periodic, dirichlet, neumann


abstract TimeSplittingSpectralMethod{T<:AbstractFloat}
abstract TimeSplittingSpectralMethodReal{T<:AbstractFloat} <: TimeSplittingSpectralMethod{T}
abstract TimeSplittingSpectralMethodReal1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal{T}
abstract TimeSplittingSpectralMethodReal2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal{T}
abstract TimeSplittingSpectralMethodReal3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodReal{T}
abstract TimeSplittingSpectralMethodComplex{T<:AbstractFloat} <: TimeSplittingSpectralMethod{T}
abstract TimeSplittingSpectralMethodComplex1D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex{T}
abstract TimeSplittingSpectralMethodComplex2D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex{T}
abstract TimeSplittingSpectralMethodComplex3D{T<:AbstractFloat} <: TimeSplittingSpectralMethodComplex{T}

abstract WaveFunction{T<:AbstractFloat}
abstract WaveFunctionReal{T<:AbstractFloat} <: WaveFunction{T}
abstract WaveFunctionReal1D{T<:AbstractFloat} <: WaveFunctionReal{T}
abstract WaveFunctionReal2D{T<:AbstractFloat} <: WaveFunctionReal{T}
abstract WaveFunctionReal3D{T<:AbstractFloat} <: WaveFunctionReal{T}
abstract WaveFunctionComplex{T<:AbstractFloat} <: WaveFunction{T}
abstract WaveFunctionComplex1D{T<:AbstractFloat} <: WaveFunctionComplex{T}
abstract WaveFunctionComplex2D{T<:AbstractFloat} <: WaveFunctionComplex{T}
abstract WaveFunctionComplex3D{T<:AbstractFloat} <: WaveFunctionComplex{T}



## Boundary conditions ############################################################################

const periodic = 0
const dirichlet = 1
const neumann = 2

none_1D(x)=0.0
none_2D(x,y)=0.0
none_3D(x,y,z)=0.0

dim(psi::WaveFunctionReal1D) = 1
dim(psi::WaveFunctionReal2D) = 2
dim(psi::WaveFunctionReal3D) = 3
dim(psi::WaveFunctionComplex1D) = 1
dim(psi::WaveFunctionComplex2D) = 2
dim(psi::WaveFunctionComplex3D) = 3
dim(psi::TimeSplittingSpectralMethodReal1D) = 1
dim(psi::TimeSplittingSpectralMethodReal2D) = 2
dim(psi::TimeSplittingSpectralMethodReal3D) = 3
dim(psi::TimeSplittingSpectralMethodComplex1D) = 1
dim(psi::TimeSplittingSpectralMethodComplex2D) = 2
dim(psi::TimeSplittingSpectralMethodComplex3D) = 3



end # module tssm_base
