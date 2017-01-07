using Combinatorics

function Base.sign(p::Array{Int,1})
    n = length(p)
    k = 0    
    for i=1:n
        for j=i+1:n
            if p[i]>p[j]
                k += 1
            end
        end
    end
   (-1)^k
end


immutable MultiFor
    k::Array{Int,1}
end

Base.start(MF::MultiFor) = Int[]

Base.done(MF::MultiFor, k::Array{Int,1}) = MF.k==k

function Base.next(MF::MultiFor, k::Array{Int,1}) 
    if k==Int[]
        k = ones(Int, length(MF.k))
        return(copy(k), k)
    end
    for i=1:length(k)
        if k[i]<MF.k[i]
            k[i] += 1
            for j = 1:i-1                 
                k[j] = 1       
            end
            return (copy(k), k)
        end
    end            
end


function find_xxx(v::Vector{Int}, lena::Int, slater_exchange, slater_indices; must_be_included::Int=0)
    res = Tuple{Int,Int,Float64}[]
    for i=1:lena
        for j=1:lena
            if must_be_included==0 || ((must_be_included in slater_indices[i]) || (must_be_included in slater_indices[i]))
                if v==slater_exchange[i,j][1]
                    push!(res, (i,j, slater_exchange[i,j][2]))
                end
            end
        end
    end
    res
end

function add_xxx!(d::Dict{Tuple{Int,Int}, Float64}, v::Tuple{Int,Int,Float64})
    key = (v[1],v[2])
    f = get!(d, key, 0.0) + v[3]
    if f==0
        delete!(d, key)
    else
        d[key] = f
    end
    d
end

function add_xxx!(d::Dict{Tuple{Int,Int}, Float64}, v::Vector{Tuple{Int,Int,Float64}})
    for x in v
        add_xxx!(d, x)
    end
    d
end

function init_mctdhf_combinatorics(f::Int, N::Int)
    slater_indices = collect(Combinatorics.Combinations(1:N, f))
    
    density_rules = [Tuple{Int,Int,Int}[] for j=1:N, l=1:N]
    for J = MultiFor([N for j=2:f])
        for j=1:N
            for l=j:N
                J1 = vcat(j,J)
                J2 = vcat(l,J)
                j1 = findfirst(slater_indices, sort(J1))
                if j1>0
                    j2 = findfirst(slater_indices, sort(J2))
                    if j2>0
                        s1 = sign(J1)
                        s2 = sign(J2)
                        t = (j1, j2, s1*s2)
                        push!(density_rules[j,l], t)
                        if j!=l
                            push!(density_rules[l,j], t)
                        end
                    end
                end
            end
        end
    end    
    
    density2_rules = [Tuple{Int,Int,Int}[] for j=1:N, l=1:N, p=1:N, q=1:N]
    if f==2
        for j=1:N
            for l=j:N
                for p=1:N
                    for q=1:N
                        J1 = vcat(j, p)
                        J2 = vcat(l, q)
                        j1 = findfirst(slater_indices, sort(J1))
                        if j1>0
                            j2 = findfirst(slater_indices, sort(J2))
                            if j2>0
                                s1 = sign(J1)
                                s2 = sign(J2)
                                t = (j1, j2, s1*s2)
                                push!(density2_rules[j,l,p,q], t)
                                if j!=l
                                    push!(density2_rules[l,j,q,p], t)
                                end
                            end
                        end
                    end
                end
            end    
        end        
    else
    for J = MultiFor([N for j=3:f])
        for j=1:N
            for l=j:N
                for p=1:N
                    for q=1:N
                        J1 = vcat(j, p, J)
                        J2 = vcat(l, q, J)
                        j1 = findfirst(slater_indices, sort(J1))
                        if j1>0
                            j2 = findfirst(slater_indices, sort(J2))
                            if j2>0
                                s1 = sign(J1)
                                s2 = sign(J2)
                                t = (j1, j2, s1*s2)
                                push!(density2_rules[j,l,p,q], t)
                                if j!=l
                                    push!(density2_rules[l,j,q,p], t)
                                end
                            end
                        end
                    end
                end
            end    
        end
    end
    end
    
    lena = binomial(N,f)
    slater_exchange = [(Int[], 1) for j=1:lena, l=1:lena]
    # entry = (exchange_index_pairs, sign)
    for j=1:lena
        for l=1:lena
            u = symdiff(slater_indices[j], slater_indices[l])
            if length(u)==2
                i1 = findfirst(slater_indices[j], u[1])
                i2 = findfirst(slater_indices[l], u[2])
                s = (-1)^(i1+i2)
                slater_exchange[j,l] = (u, s)
            elseif length(u)==4
                u = [u[1], u[3], u[2], u[4]]
                i1 = findfirst(slater_indices[j], u[1])
                i2 = findfirst(slater_indices[l], u[2])
                i3 = findfirst(slater_indices[j], u[3])
                i4 = findfirst(slater_indices[l], u[4])
                s = (-1)^(i1+i2+i3+i4)
                slater_exchange[j,l] = (u, s)
            end
        end
    end 
    
    res = [Dict{Tuple{Int,Int},Float64}([]) for p=1:N, q=1:N]
    for p=1:N
        for q=1:N
            if p==q 
                for j=1:lena
                    if p in slater_indices[j]
                       add_xxx!(res[p,q], (j,j, 1.0))
                    end
                end
            end
            v = find_xxx([q,p], lena, slater_exchange, slater_indices)
            add_xxx!(res[p,q], v)
        end
    end
    slater1_rules = [[(key[1],key[2],val) for (key,val) in res[p,q]] for p=1:N, q=1:N]
    
    res = [Dict{Tuple{Int,Int},Float64}([]) for p=1:N, q=1:N, r=1:N, s=1:N]
    for p=1:N
        for q=1:N
            for r=1:N
                for s=1:N
                    if p==q && r==s && p!=r
                        for j=1:lena
                            if p in slater_indices[j] && r in slater_indices[j]                       
                                add_xxx!(res[p,q,r,s], (j,j, 0.5))
                            end
                        end
                    end
                    if p==s && r==q && p!=r
                        for j=1:lena
                            if p in slater_indices[j] && r in slater_indices[j]                       
                                add_xxx!(res[p,q,r,s], (j,j, -0.5)) 
                            end
                        end
                    end
                    if p==q && r!=s 
                        v = find_xxx([s,r], lena, slater_exchange, slater_indices, must_be_included=p)
                        add_xxx!(res[p,q,r,s], v)
                    end
                    if q==r && s!=p 
                        v0 = find_xxx([s,p], lena, slater_exchange, slater_indices, must_be_included=q)
                        v = [(j,k,-sigma) for (j,k, sigma) in v0]
                        add_xxx!(res[p,q,r,s], v)     
                    end     
                    v = find_xxx([s,r,q,p], lena, slater_exchange, slater_indices)
                    add_xxx!(res[p,q,r,s], v)
                    v0 = find_xxx([s,p,q,r], lena, slater_exchange, slater_indices)
                    v = [(j,k,-sigma) for (j,k, sigma) in v0]
                    add_xxx!(res[p,q,r,s], v)
                end
            end
        end
    end
    slater2_rules = [[(key[1],key[2],val) for (key,val) in res[p,q,r,s]] for p=1:N, q=1:N, r=1:N, s=1:N]
    
    orthogonalization_rules = [Tuple{Int,Int,Int}[] for p=1:N, q=1:N]
    for p=1:N
        for j=1:lena
            k = findfirst(slater_indices[j], p)
            if k>0
                for q=1:p-1
                    J = copy(slater_indices[j])
                    J[k] = q
                    J1 = sort(J)
                    l = findfirst(slater_indices, J1)
                    if l>0
                        s = sign(J)*sign(J1)
                        push!(orthogonalization_rules[p,q], (j,l, s))
                    end
                end
            end
        end
    end    
    
    slater_indices, density_rules, density2_rules, slater_exchange, slater1_rules, slater2_rules, orthogonalization_rules
end


function init_Vee(x::Vector{Float64})
    n = length(x)
    1./sqrt(( kron(x',ones(n)) - kron(x,ones(1,n)) ).^2 + 1)    
end


type MCTDHF1D <: TSSM.TimeSplittingSpectralMethodComplex1D
    m::Schroedinger1D
    f::Int # number of electrons
    N::Int # number of orbitals
    lena::Int # number of (independent) coefficients
    spins::Array{Int, 1}
    
    slater_indices
    density_rules 
    density2_rules
    slater_exchange
    slater1_rules
    slater2_rules
    orthogonalization_rules

    Vee
    density_matrix
    density2_tensor
    k1
    k2

    function MCTDHF1D(f::Integer, N::Integer, 
                      nx::Integer, xmin::Real, xmax::Real; spins::Array{Int,1}=ones(Int, N))
        m = Schroedinger1D(nx, xmin, xmax)
        lena = binomial(N,f)
        slater_indices, density_rules, density2_rules, slater_exchange, slater1_rules, slater2_rules, orthogonalization_rules = 
              init_mctdhf_combinatorics(f, N)
        Vee = init_Vee(get_nodes(m))
        new(m, f, N, lena, spins, slater_indices, density_rules, 
            density2_rules, slater_exchange, slater1_rules, slater2_rules, orthogonalization_rules, Vee,
        zeros(Complex{Float64},N,N),  zeros(Complex{Float64},N,N,N,N), nothing, nothing )
    end
end

type Orbital
    phi::WfSchroedinger1D
    spin::Int 
end

type WfMCTDHF1D <: TSSM.WaveFunctionComplex1D
    a::Array{Complex{Float64}, 1}
    o::Array{Orbital, 1}
    m::MCTDHF1D
    function WfMCTDHF1D(m::MCTDHF1D)
        a = zeros(Complex{Float64}, m.lena)
        o = [Orbital(WfSchroedinger1D(m.m), m.spins[j]) for j=1:m.N]
        new(a, o, m)
    end
end

function TSSM.wave_function(m::MCTDHF1D )
    WfMCTDHF1D(m) 
end

Base.norm(o::Orbital) = Base.norm(o.phi)
inner_product(o1::Orbital, o2::Orbital) = (o1.spin==o2.spin ? TSSM.inner_product(o1.phi, o2.phi) : 0.0im)
potential_matrix_element(o1::Orbital, o2::Orbital) = (o1.spin==o2.spin ? TSSM.potential_matrix_element(o1.phi, o2.phi) : 0.0im)
kinetic_matrix_element(o1::Orbital, o2::Orbital) = (o1.spin==o2.spin ? TSSM.kinetic_matrix_element(o1.phi, o2.phi) : 0.0im)

function Base.scale!(o::Orbital, f::Number) 
    TSSM.scale!(o.phi, f)
end

function axpy!(o1::Orbital, o2::Orbital, f::Number)
    TSSM.axpy!(o1.phi, o2.phi, f)
end


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
                h += conj(getindex(A, p, J...)*getindex(A, q, J...))
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
                            h += conj(getindex(A, p, r, J...)*getindex(A, q, s, J...))
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


function gen_density_matrix(psi::WfMCTDHF1D)
    N = psi.m.N
    rho = psi.m.density_matrix
    rho[:,:] = 0.0
    for j=1:N
        for l=j:N
            for (u, v, s) in psi.m.density_rules[j,l]
                rho[j,l] += s * conj(psi.a[u]) * psi.a[v]
                if l!=j
                   rho[l,j] = conj(rho[j,l])
                end
            end    
        end
    end
    nothing
end


function gen_density2_tensor(psi::WfMCTDHF1D; mult_inverse_density_matrix::Bool=true)
    N = psi.m.N
    rho = psi.m.density2_tensor
    rho[:,:,:,:] = 0.0
    for j=1:N
        for l=j:N
            for p=1:N
                for q=1:N
                    for (u, v, s) in psi.m.density2_rules[j,l,p,q]
                        rho[j,l,p,q] += s * conj(psi.a[u]) * psi.a[v]
                        if l!=j
                           rho[l,j,q,p] = conj(rho[j,l,p,q])
                        end
                    end
                end
            end    
        end
    end
    if mult_inverse_density_matrix
        for p=1:N
            for q=1:N
                rho[:,:,p,q] = psi.m.density_matrix \ rho[:,:,p,q]
            end
        end
    end
    nothing
end


norm(psi::WfMCTDHF1D) = Base.norm(psi.a)
#Note, only correct if orbitals are ortonormal

function TSSM.to_real_space!(psi::WfMCTDHF1D)
    for j=1:psi.m.N
        to_real_space!(psi.o[j].phi)
    end
end


function TSSM.to_frequency_space!(psi::WfMCTDHF1D)
    for j=1:psi.m.N
        to_frequency_space!(psi.o[j].phi)
    end
end


function set_zero!(psi::WfMCTDHF1D)
    for j=1:psi.m.N 
        get_data(psi.o[j].phi, true)[:] = 0.0
    end
    psi.a[:] = 0.0
end



function gen_rhs1!(rhs::WfMCTDHF1D, psi::WfMCTDHF1D)
    m = rhs.m
    if m ≠ psi.m
        error("rhs and psi must belong to the same method")
    end
    for q=1:m.N
        add_apply_B!(psi.o[q].phi, rhs.o[q].phi, 1im)
        for p=1:m.N
            if psi.o[p].spin==psi.o[q].spin
                h = inner_product(psi.o[p], rhs.o[q])
                for (j,l,f) in m.slater1_rules[q,p]
                    rhs.a[j] += h*f*psi.a[l] 
                end
            end
        end
    end
end

function project_out_orbitals!(rhs:: WfMCTDHF1D, psi::WfMCTDHF1D)
    m = psi.m
    c = zeros(Complex{Float64},m.N)
    for p = 1:m.N
        for q = 1:m.N
            c[q] = inner_product(psi.o[q], rhs.o[p])
        end
        for q = 1:m.N
            axpy!(rhs.o[p], psi.o[q], -c[q])
        end
    end            
end



function gen_rhs2!(rhs::WfMCTDHF1D, psi::WfMCTDHF1D)
    m = rhs.m
    if m ≠ psi.m
        error("rhs and psi must belong to the same method")
    end
    n = get_nx(m.m)
    dx = (get_xmax(m.m)-get_xmin(m.m))/n
    u_pq = zeros(Complex{Float64}, n)
    u_pqs = zeros(Complex{Float64}, n)
    to_real_space!(psi)
    to_real_space!(rhs)
    for p=1:m.N
        for q=1:m.N
            if psi.o[p].spin==psi.o[q].spin
                u_pq[:] = m.Vee * (conj(get_data(psi.o[p].phi, true)).*get_data(psi.o[q].phi, true))
                for s=1:m.N
                    u_pqs[:] = u_pq .* get_data(psi.o[s].phi, true)
                    for r=1:m.N
                        if psi.o[r].spin==psi.o[s].spin
                            u = get_data(rhs.o[r].phi, true)
                            u[:] += (m.density2_tensor[r,s,p,q]*(m.f-1)*dx) * u_pqs                
                            h = dot(get_data(psi.o[r].phi, true), u_pqs) * dx^2
                            for (j,l,f) in m.slater2_rules[q,p,s,r]
                                rhs.a[j] += h*f*psi.a[l] 
                            end
                        end
                    end
                end
            end
        end
    end
end


function gen_rhs!(rhs::WfMCTDHF1D, psi::WfMCTDHF1D)
    m = rhs.m
    if m ≠ psi.m
        error("rhs and psi must belong to the same method")
    end
    #gen_density_matrix(psi)
    #gen_density2_tensor(psi)
    set_zero!(rhs)
    gen_rhs1!(rhs, psi)
    #gen_rhs2!(rhs, psi)
    project_out_orbitals!(rhs, psi)
end


type Schroedinger2Electrons <: TSSM.TimeSplittingSpectralMethodComplex2D
    m::Schroedinger2D
    function Schroedinger2Electrons(nx::Integer, xmin::Real, xmax::Real)
        new(Schroedinger2D(nx, xmin, xmax, nx, xmin, xmax))
    end
end

type WfSchroedinger2Electrons <: TSSM.WaveFunctionComplex2D
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
    u_singlet[:,:] = 0.0
    u_triplet_up[:,:] = 0.0
    u_triplet_down[:,:] = 0.0
    u_triplet_symm[:,:] = 0.0
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


function orthonormalize_orbitals!(psi::WfMCTDHF1D)
    m = psi.m
    g = zeros(Complex{Float64}, m.N)
    for p = 1:m.N
        a1 = zeros(Complex{Float64}, m.lena)
        for q=1:p-1
            g[q] = inner_product(psi.o[q], psi.o[p])
        end
        for q=1:p-1
            if psi.o[p].spin==psi.o[q].spin
                axpy!(psi.o[p], psi.o[q], -g[q])
                for (j, l, s) in m.orthogonalization_rules[p,q]
                    a1[l] += s*g[q]*psi.a[j]
                end
            end    
        end
        f = Base.norm(psi.o[p])
        scale!(psi.o[p],1/f)
        for j=1:m.lena
            if p in m.slater_indices[j]
                a1[j] += f*psi.a[j]
            else
                a1[j] += psi.a[j]
            end
        end
        psi.a[:] = a1
    end
    psi
end


function potential_energy_1(psi::WfMCTDHF1D)
    m = psi.m
    V = 0
    for p=1:m.N        
        for q=1:m.N
            if psi.o[p].spin==psi.o[q].spin
                h = potential_matrix_element(psi.o[p], psi.o[q])
                for (j,l,f) in m.slater1_rules[q,p]
                    V += h*f*conj(psi.a[j])*psi.a[l]
                    #V += h*f*psi.a[j]*conj(psi.a[l])
                end
            end
        end
    end
    real(V)
end


function potential_energy_1_A(psi::WfMCTDHF1D)
    m = psi.m
    V = 0
    for p=1:m.N    
        h = potential_matrix_element(psi.o[p], psi.o[p])
        for (j,l,f) in m.slater1_rules[p,p]
           #V += real(h*f*conj(psi.a[j])*psi.a[l])
            V += real(h*f*psi.a[j]*conj(psi.a[l]))
        end
        for q=1:p-1
            if psi.o[p].spin==psi.o[q].spin
                h = potential_matrix_element(psi.o[p], psi.o[q])
                for (j,l,f) in m.slater1_rules[p,q]
                    #V += 2*real(h*f*conj(psi.a[j])*psi.a[l])
                    V += 2*real(h*f*psi.a[j]*conj(psi.a[l]))
                end
            end
        end
    end
    real(V)
end



function TSSM.kinetic_energy(psi::WfMCTDHF1D)
    m = psi.m
    T = 0
    for p=1:m.N        
        for q=1:m.N
            if psi.o[p].spin==psi.o[q].spin
                h = kinetic_matrix_element(psi.o[p], psi.o[q])
                for (j,l,f) in m.slater1_rules[q,p]
                    T += h*f*conj(psi.a[j])*psi.a[l]
                    #T += h*f*psi.a[j]*conj(psi.a[l])
                end
            end
        end
    end
    real(T)    
end


function kinetic_energy_A(psi::WfMCTDHF1D)
    m = psi.m
    V = 0
    for p=1:m.N    
        h = kinetic_matrix_element(psi.o[p], psi.o[p])
        for (j,l,f) in m.slater1_rules[p,p]
           #V += real(h*f*conj(psi.a[j])*psi.a[l])
            V += real(h*f*psi.a[j]*conj(psi.a[l]))
        end
        for q=1:p-1
            if psi.o[p].spin==psi.o[q].spin
                h = kinetic_matrix_element(psi.o[p], psi.o[q])
                for (j,l,f) in m.slater1_rules[p,q]
                    #V += 2*real(h*f*conj(psi.a[j])*psi.a[l])
                    V += 2*real(h*f*psi.a[j]*conj(psi.a[l]))
                end
            end
        end
    end
    V
end


function potential_energy_2(psi::WfMCTDHF1D)
    m = psi.m
    V = 0
    n = get_nx(m.m)
    dx = (get_xmax(m.m)-get_xmin(m.m))/n
    u_pq = zeros(Complex{Float64}, n)
    u_pqs = zeros(Complex{Float64}, n)
    to_real_space!(psi)
    for p=1:m.N
        for q=1:m.N
            if psi.o[p].spin==psi.o[q].spin
                u_pq[:] = m.Vee * (conj(get_data(psi.o[p].phi, true)).*get_data(psi.o[q].phi, true))
                for s=1:m.N
                    u_pqs[:] = u_pq .* get_data(psi.o[s].phi, true)
                    for r=1:m.N
                        if psi.o[r].spin==psi.o[s].spin
                            h = dot(get_data(psi.o[r].phi, true), u_pqs)
                            for (j,l,f) in m.slater2_rules[q,p,s,r]
                                V += h*f*conj(psi.a[j])*psi.a[l]
                                #V += h*f*psi.a[j]*conj(psi.a[l])
                            end
                        end
                    end
                end
            end
        end
    end
    V *= dx^2
    real(V)
end


TSSM.potential_energy(psi::WfMCTDHF1D) = potential_energy_1(psi) + potential_energy_2(psi)


function TSSM.imaginary_time_propagate_A!(psi::WfMCTDHF1D, dt::Real)
    for j=1:psi.m.N
        imaginary_time_propagate_A!(psi.o[j].phi, dt)
    end
end

function TSSM.imaginary_time_propagate_B!(psi::WfMCTDHF1D, dt::Real)
    for j=1:psi.m.N
        imaginary_time_propagate_B!(psi.o[j].phi, dt)
    end
end

function TSSM.scale!(psi::WfMCTDHF1D, f::Number)
    for j=1:psi.m.N
        scale!(psi.o[j], f)
    end
    psi.a[:] *= f
end

function axpy!(psi1::WfMCTDHF1D, psi2::WfMCTDHF1D, f::Number)
    for j=1:psi.m.N
        axpy!(psi1.o[j], psi2.o[j], f)
    end
    psi1.a[:] += f*psi2.a[:]
end

function RK2_step!(psi::WfMCTDHF1D, dt::Number)
    m = psi.m
    gen_rhs!(m.k1, psi)
    scale!(m.k1, -0.5im*dt)
    axpy!(m.k1, psi, 1.0)
    #orthonormalize_orbitals!(m.k1)
    gen_rhs!(m.k2, m.k1)
    axpy!(psi, m.k2, -1im*dt)
end


function groundstate!(psi::WfMCTDHF1D, dt::Real, n::Int)
    m = psi.m
    m.k1 = wave_function(m)
    m.k2 = wave_function(m)
    
    orthonormalize_orbitals!(psi)
    for k=1:n
        imaginary_time_propagate_A!(psi, 0.5*dt)
        orthonormalize_orbitals!(psi)
        RK2_step!(psi, -1im*dt)
        #imaginary_time_propagate_B!(psi, dt)
        orthonormalize_orbitals!(psi)
        imaginary_time_propagate_A!(psi, 0.5*dt)
        orthonormalize_orbitals!(psi)
        norm_psi = norm(psi)
        psi.a[:] *= 1/norm_psi
        
        E_pot = potential_energy(psi)
        E_kin = kinetic_energy(psi)
        E = E_pot + E_kin
        println("step =", k,"  E_pot =", E_pot, "  E_kin=", E_kin, "  E=", E)                
    end
    
    m.k1 = nothing
    m.k2 = nothing
end
    
