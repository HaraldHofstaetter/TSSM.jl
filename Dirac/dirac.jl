mutable struct Dirac1D <:  TSSM.TimeSplittingSpectralMethodComplex1D{Float64}
    m::Schroedinger1D
    epsilon::Float64
    delta::Float64
    nu::Float64
    d::Array{Float64, 1}
    Qd::Array{Float64, 1}
    Qo::Array{Float64, 1}
    function Dirac1D(nx::Integer, xmin::Real, xmax::Real; epsilon::Real=1, delta::Real=1, nu::Real=1,
                     potential::Function=TSSM.none_1D) #, potential_t::Function=TSSM.none_2D)
        m = Schroedinger1D(nx, xmin, xmax, potential=potential)
        #CHECK: maybe the order of elements in mu is different!!!
        mu = [2*pi*l/(xmax-xmin) for l =vcat(0:(div(nx,2)-1), (-div(nx,2):-1))] # nx should be even        
        #mu = [2*pi*l/(xmax-xmin) for l =-div(nx,2):div(nx,2)-1] # nx should be even        
        eta = [sqrt(nu^2+delta^2*epsilon^2*x^2) for x=mu] 
        d = [x/(delta*epsilon^2) for x=eta]
        f = [1/sqrt(2*x*(x+nu)) for x=eta]
        Qd = f.*[x+nu for x=eta]
        Qo = f.*[delta*epsilon*x for x=mu]
        new(m, epsilon, delta, nu, d, Qd, Qo)
    end
end

mutable struct WfDirac1D <: TSSM.WaveFunctionComplex1D{Float64}
    m::Dirac1D
    psi1::WfSchroedinger1D
    psi2::WfSchroedinger1D
    function WfDirac1D(m::Dirac1D)
        psi1 = WfSchroedinger1D(m.m)
        psi2 = WfSchroedinger1D(m.m)
        new(m, psi1, psi2)
    end
end


function TSSM.wave_function(m::Dirac1D )
    WfDirac1D(m) 
end

set_potential!(m::Dirac1D, V::Function) = TSSM.set_potential!(m.m, V)
set_potential_t!(m::Dirac1D, V::Function) = TSSM.set_potential_t!(m.m, V)




function TSSM.to_real_space!(psi::WfDirac1D)
    to_real_space!(psi.psi1)
    to_real_space!(psi.psi2)
end

function TSSM.to_frequency_space!(psi::WfDirac1D)
    to_frequency_space!(psi.psi1)
    to_frequency_space!(psi.psi2)
end

function set!(psi::WfDirac1D, f1::Function, f2::Function)
    TSSM.set!(psi.psi1, f1)
    TSSM.set!(psi.psi2, f2)
end

function distance(psi1::WfDirac1D, psi2::WfDirac1D)  
    m = psi1.m
    if m ≠ psi2.m
        error("psi1 and psi2 must belong to the same method")
    end
    sqrt(TSSM.distance(psi1.psi1, psi2.psi1)^2+TSSM.distance(psi1.psi1, psi2.psi1)^2)
end

function TSSM.set_time!(psi::WfDirac1D, t::Number)
    set_time!(psi.psi1, t)
    set_time!(psi.psi2, t)
end

TSSM.get_time(psi::WfDirac1D) = get_time(psi.psi1)
TSSM.set_propagate_time_together_with_A!(m::Dirac1D, flag::Bool) = set_propagate_time_together_with_A!(m.m, flag)
TSSM.get_propagate_time_together_with_A(m::Dirac1D) = get_propagate_time_together_with_A(m.m)

function TSSM.copy!(psi1::WfDirac1D, psi2::WfDirac1D)
    TSSM.copy!(psi1.psi1, psi2.psi1)
    TSSM.copy!(psi1.psi2, psi2.psi2)
end


function TSSM.propagate_B!(psi::WfDirac1D, dt::Number)
    propagate_B!(psi.psi1, dt)
    propagate_B!(psi.psi2, dt)
end

function TSSM.propagate_A!(psi::WfDirac1D, dt::Number)
    to_frequency_space!(psi)
    u1 = get_data(psi.psi1, true)
    u2 = get_data(psi.psi2, true)
    m = psi.m
    for l=1:length(u1)
        ee = exp(-1im*dt*m.d[l])
        h1 = ee*(m.Qd[l]*u1[l]+m.Qo[l]*u2[l])
        h2 =   (-m.Qo[l]*u1[l]+m.Qd[l]*u2[l])/ee
        u1[l] = (m.Qd[l]*h1-m.Qo[l]*h2)
        u2[l] = (m.Qo[l]*h1+m.Qd[l]*h2)
    end
end
