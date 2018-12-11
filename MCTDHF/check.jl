function get_full_coeffs(psi::WfMCTDHF1D)
    m = psi.m
    A = zeros(Complex{Float64}, [m.N for j=1:m.f]...)
    perms = [(p, sign(p)) for p in Combinatorics.Permutations(1:m.f,m.f)]
    for j=1:m.lena
        J = m.slater_indices[j]
        for (p,s) in perms
            J1 = [J[p[i]] for i=1:m.f]
            setindex!(A,s*psi.a[j], J1...)
        end
    end
    A
end

function density_matrix_from_full_coeffs(psi::WfMCTDHF1D)
    m = psi.m
    A = get_full_coeffs(psi)
    rho = zeros(Complex{Float64},m.N,m.N)
    for p=1:m.N
        for q=1:m.N
            h = 0im
            for J = MultiFor([m.N for j=2:m.f])
                h += conj(getindex(A, p, J...))*getindex(A, q, J...)
            end
            rho[p,q] = h
        end
    end
    rho
end

function density2_tensor_from_full_coeffs(psi::WfMCTDHF1D)
    m = psi.m
    A = get_full_coeffs(psi)
    rho = zeros(Complex{Float64},m.N,m.N,m.N,m.N)
    for p=1:m.N
        for q=1:m.N
            for r=1:m.N
                for s=1:m.N
                    h = 0im
                    if m.f>=3
                        for J = MultiFor([m.N for j=3:m.f])
                            h += conj(getindex(A, p, r, J...))*getindex(A, q, s, J...)
                        end
                    else
                        h += conj(A[p,r])*A[q,s]
                    end
                    rho[p,q,r,s] = h
                end
            end
        end
    end
    rho
end


mutable struct Schroedinger2Electrons <: TSSM.TimeSplittingSpectralMethodComplex2D{Float64}
    m::Schroedinger2D
    function Schroedinger2Electrons(nx::Integer, xmin::Real, xmax::Real; 
                                    potential::Function=TSSM.none_2D,
                                    potential_t::Function=TSSM.none_3D)
        new(Schroedinger2D(nx, xmin, xmax, nx, xmin, xmax, potential=potential, potential_t=potential_t))
    end
end

TSSM.set_potential!(m::Schroedinger2Electrons, V::Function) = set_potential!(m.m, V)

mutable struct WfSchroedinger2Electrons <: TSSM.WaveFunctionComplex2D{Float64}
    m::Schroedinger2Electrons
    singlet::WfSchroedinger2D
    triplet_up::WfSchroedinger2D
    triplet_down::WfSchroedinger2D
    triplet_symm::WfSchroedinger2D
    function WfSchroedinger2Electrons(m::Schroedinger2Electrons)
        new(m, WfSchroedinger2D(m.m),
               WfSchroedinger2D(m.m),
               WfSchroedinger2D(m.m),
               WfSchroedinger2D(m.m))
    end
end

function TSSM.wave_function(m::Schroedinger2Electrons)
    WfSchroedinger2Electrons(m) 
end



function convert_to_full!(psi2::WfSchroedinger2Electrons, psi::WfMCTDHF1D)
    m = psi.m
    n = get_nx(m.m)
    @assert m.f==2
    u_singlet = get_data(psi2.singlet, true)
    u_triplet_up = get_data(psi2.triplet_up, true)
    u_triplet_down = get_data(psi2.triplet_down, true)
    u_triplet_symm = get_data(psi2.triplet_symm, true)
    u_singlet[:,:] .= 0.0
    u_triplet_up[:,:] .= 0.0
    u_triplet_down[:,:] .= 0.0
    u_triplet_symm[:,:] .= 0.0
    f = 1/sqrt(2)
    to_real_space!(psi)
    for j = 1:m.lena # eval slater determinants
        J = m.slater_indices[j]
        s = sign(J) 
        v1 = get_data(psi.o[J[1]].phi, true)
        v2 = get_data(psi.o[J[2]].phi, true)
        s1 = psi.o[J[1]].spin 
        s2 = psi.o[J[2]].spin
        if s1==s2
            if s1==+1
                for i1=1:n
                    for i2=1:n
                        u_triplet_up[i1,i2] += (s*f*psi.a[j])*(v1[i1]*v2[i2] - v1[i2]*v2[i1])
                    end
                end
            else
                for i1=1:n
                    for i2=1:n
                        u_triplet_down[i1,i2] += (s*f*psi.a[j])*(v1[i1]*v2[i2] - v1[i2]*v2[i1])
                    end
                end
            end
        else
            if s1==+1
                for i1=1:n
                    for i2=1:n
                        u_singlet[i1,i2] += (0.5*s*psi.a[j])*(v1[i1]*v2[i2] + v1[i2]*v2[i1])
                        u_triplet_symm[i1,i2] += (0.5*s*psi.a[j])*(v1[i1]*v2[i2] - v1[i2]*v2[i1])
                    end
                end
            else
                for i1=1:n
                    for i2=1:n
                        u_singlet[i1,i2] -= (0.5*s*psi.a[j])*(v1[i1]*v2[i2] + v1[i2]*v2[i1])
                        u_triplet_symm[i1,i2] += (0.5*s*psi.a[j])*(v1[i1]*v2[i2] - v1[i2]*v2[i1])
                    end
                end
            end
        end
    end    
    psi2
end

TSSM.norm(psi::WfSchroedinger2Electrons) = sqrt(
    TSSM.norm(psi.singlet)^2 + TSSM.norm(psi.triplet_symm)^2 +
    TSSM.norm(psi.triplet_up)^2 + TSSM.norm(psi.triplet_down)^2 )

TSSM.potential_energy(psi::WfSchroedinger2Electrons) = (
    potential_energy(psi.singlet) + potential_energy(psi.triplet_symm) +
    potential_energy(psi.triplet_up) + potential_energy(psi.triplet_down) )

TSSM.kinetic_energy(psi::WfSchroedinger2Electrons) = (
    kinetic_energy(psi.singlet) + kinetic_energy(psi.triplet_symm) +
    kinetic_energy(psi.triplet_up) + kinetic_energy(psi.triplet_down) )



function expand_slater_determinants!(psi3::WfSchroedinger3D, psi::WfMCTDHF1D)
    m = psi.m
    @assert m.f==3
    m3 = psi3.m
    u3 = get_data(psi3, true)
    u3[:,:,:] = 0.0
    f = 1/sqrt(6)
    to_real_space!(psi)
    for j = 1:m.lena # eval slater determinants
        J = m.slater_indices[j]
        s = sign(J)
        v1 = get_data(psi.o[J[1]].phi, true)
        v2 = get_data(psi.o[J[2]].phi, true)
        v3 = get_data(psi.o[J[3]].phi, true)
        #TODO: spin-dependent signs !!!        
        for i1=1:size(u3,1)
            for i2=1:size(u3,2)
                for i3=1:size(u3,3)
                    u3[i1,i2,i3] += (s*f*psi.a[j])*(v1[i1]*v2[i2]*v3[i3] + v1[i2]*v2[i3]*v3[i1] + v1[i3]*v2[i1]*v3[i2]
                                                   -v1[i3]*v2[i2]*v3[i1] - v1[i1]*v2[i3]*v3[i2] - v1[i2]*v2[i1]*v3[i3])
                end
            end
        end
    end    
    psi3
end

