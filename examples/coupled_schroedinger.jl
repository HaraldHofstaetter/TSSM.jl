type CoupledSchroedinger1D <: TSSM.TimeSplittingSpectralMethodComplex1D
    m1::Schroedinger1D
    m2::Schroedinger1D
    delta
    e
    function CoupledSchroedinger1D(nx::Integer, xmin::Real, xmax::Real,
                        delta::Real, e::Real)
       m1 = Schroedinger1D(nx, xmin, xmax)
       m2 = Schroedinger1D(nx, xmin, xmax)
       ev1 = get_eigenvalues(m1, true)
       ev2 = get_eigenvalues(m2, true)
       f = 2*2*pi/(xmax-xmin)
       ev1[1:end] = ev1 - f*delta*[0:nx/2-1;-nx/2:-1]
       ev2[1:end] = ev2 + f*delta*[0:nx/2-1;-nx/2:-1]
       new(m1, m2, delta, e) 
    end                    
end

type WfCoupledSchroedinger1D <: TSSM.WaveFunctionComplex1D
   psi1::WfSchroedinger1D
   psi2::WfSchroedinger1D
   m::CoupledSchroedinger1D
   function WfCoupledSchroedinger1D( m::CoupledSchroedinger1D )
      new( WfSchroedinger1D(m.m1), WfSchroedinger1D(m.m2), m )
   end
end   

function TSSM.wave_function(m::CoupledSchroedinger1D )
    WfCoupledSchroedinger1D(m) 
end

function TSSM.propagate_A!(psi::WfCoupledSchroedinger1D, dt::Number)
    propagate_A!(psi.psi1, dt)
    propagate_A!(psi.psi2, dt)
end

function TSSM.copy!(target::WfCoupledSchroedinger1D, source::WfCoupledSchroedinger1D)
    copy!(target.psi1, source.psi1)
    copy!(target.psi2, source.psi2)
end

function TSSM.scale!(psi::WfCoupledSchroedinger1D, factor::Number)
    scale!(psi.psi1, factor)
    scale!(psi.psi2, factor)
end

function TSSM.axpy!(psi1::WfCoupledSchroedinger1D, psi2::WfCoupledSchroedinger1D, factor::Number)
    axpy!(psi1.psi1, psi2.psi1, factor)
    axpy!(psi1.psi2, psi2.psi2, factor)
end

type Function2
    f1::Function
    f2::Function
end

function TSSM.set!(psi::WfCoupledSchroedinger1D, f::Function2)
    set!(psi.psi1, f.f1)
    set!(psi.psi2, f.f2)
end

function TSSM.set!(psi::WfCoupledSchroedinger1D, f::Function2, t::Real)
    set!(psi.psi1, f.f1, t)
    set!(psi.psi2, f.f2, t)
end


function TSSM.norm(psi::WfCoupledSchroedinger1D )
   sqrt(norm(psi.psi1)^2+norm(psi.psi2)^2)
end

function TSSM.distance(psi1::WfCoupledSchroedinger1D, psi2::WfCoupledSchroedinger1D)
   sqrt(distance(psi1.psi1, psi2.psi1)^2+distance(psi1.psi2,psi2.psi2)^2)
end

function TSSM.propagate_B!(psi::WfCoupledSchroedinger1D, dt::Number)
   u1 = get_data(psi.psi1, true)
   u2 = get_data(psi.psi2, true)
   to_real_space!(psi.psi1)
   to_real_space!(psi.psi2)
   u1[1:end] = exp(1im*dt*(abs(u1).^2 + psi.m.e*abs(u2).^2)).*u1
   u2[1:end] = exp(1im*dt*(psi.m.e*abs(u1).^2 + abs(u2).^2)).*u2
end


