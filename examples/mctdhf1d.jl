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


function find_xxx(v::Vector{Int}, lena::Int, slater_exchange)
    res = Tuple{Int,Int,Float64}[]
    for i=1:lena
        for j=1:lena
            if v==slater_exchange[i,j][1]
                push!(res, (i,j, slater_exchange[i,j][2]))
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
            v = find_xxx([q,p], lena, slater_exchange)
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
                        v = find_xxx([s,r], lena, slater_exchange)
                        add_xxx!(res[p,q,r,s], v)
                    end
                    if q==r && s!=p 
                        v0 = find_xxx([s,p], lena, slater_exchange)
                        v = [(j,k,-sigma) for (j,k, sigma) in v0]
                        add_xxx!(res[p,q,r,s], v)     
                    end     
                    v = find_xxx([s,r,q,p], lena, slater_exchange)
                    add_xxx!(res[p,q,r,s], v)
                    v0 = find_xxx([s,p,q,r], lena, slater_exchange)
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
                      nx::Integer, xmin::Real, xmax::Real)
        m = Schroedinger1D(nx, xmin, xmax)
        lena = binomial(N,f)
        slater_indices, density_rules, density2_rules, slater_exchange, slater1_rules, slater2_rules, orthogonalization_rules = 
              init_mctdhf_combinatorics(f, N)
        Vee = init_Vee(get_nodes(m))
        new(m, f, N, lena, slater_indices, density_rules, 
            density2_rules, slater_exchange, slater1_rules, slater2_rules, orthogonalization_rules, Vee,
        zeros(Complex{Float64},N,N),  zeros(Complex{Float64},N,N,N,N), nothing, nothing )
    end
end

type WfMCTDHF1D <: TSSM.WaveFunctionComplex1D
    a::Array{Complex{Float64}, 1}
    phi::Array{WfSchroedinger1D, 1}
    m::MCTDHF1D
    function WfMCTDHF1D(m::MCTDHF1D)
        a = zeros(Complex{Float64}, m.lena)
        phi = [WfSchroedinger1D(m.m) for j=1:m.N]
        new(a, phi, m)
    end
end

function TSSM.wave_function(m::MCTDHF1D )
    WfMCTDHF1D(m) 
end


function gen_density_matrix(psi::WfMCTDHF1D)
    N = psi.m.N
    rho = psi.m.density_matrix
    rho[:,:] = 0.0
    for j=1:N
        for l=1:N
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
#Note, only correct if psi.phi's are ortonormal

function TSSM.to_real_space!(psi::WfMCTDHF1D)
    for j=1:psi.m.N
        to_real_space!(psi.phi[j])
    end
end


function TSSM.to_frequency_space!(psi::WfMCTDHF1D)
    for j=1:psi.m.N
        to_frequency_space!(psi.phi[j])
    end
end


function set_zero!(psi::WfMCTDHF1D)
    for j=1:psi.m.N 
        get_data(psi.phi[j], true)[:] = 0.0
    end
    psi.a[:] = 0.0
end



function gen_rhs1!(rhs::WfMCTDHF1D, psi::WfMCTDHF1D)
    m = rhs.m
    if m ≠ psi.m
        error("rhs and psi must belong to the same method")
    end
    #n = get_nx(m.m)
    #u_save = zeros(Complex{Float64}, n)    
    for q=1:m.N
        #to_real_space!(rhs.phi[q])
        #u_save[:] = get_data(rhs.phi[q], false)
        #u = get_data(rhs.phi[q], true)  
        #u[:] = 0.0
        #if A
        #    add_apply_A!(psi.phi[q], rhs.phi[q], 1im)
        #end
        #if B
            add_apply_B!(psi.phi[q], rhs.phi[q], 1im)
        #end    
        #to_real_space!(rhs.phi[q])
        for p=1:m.N
            h = inner_product(psi.phi[p], rhs.phi[q])
            for (j,l,f) in m.slater1_rules[q,p]
                rhs.a[j] += h*f*psi.a[l] 
            end
        end
        #u[:] += u_save        
    end
end

function project_out_orbitals!(rhs:: WfMCTDHF1D, psi::WfMCTDHF1D)
    m = psi.m
    c = zeros(Complex{Float64},m.N)
    for p = 1:m.N
        for q = 1:m.N
            c[q] = inner_product(psi.phi[q], rhs.phi[p])
        end
        for q = 1:m.N
            axpy!(rhs.phi[p], psi.phi[q], -c[q])
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
            u_pq[:] = m.Vee * (conj(get_data(psi.phi[p], true)).*get_data(psi.phi[q], true))
            for s=1:m.N
                u_pqs[:] = u_pq .* get_data(psi.phi[s], true)
                for r=1:m.N
                    u = get_data(rhs.phi[r], true)
                    u[:] += (m.density2_tensor[r,s,p,q]*(m.f-1)*dx) * u_pqs                
                    h = dot(get_data(psi.phi[r], true), u_pqs) * dx^2
                    for (j,l,f) in m.slater2_rules[q,p,s,r]
                        rhs.a[j] += h*f*psi.a[l] 
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


function expand_slater_determinants!(psi2::WfSchroedinger2D, psi::WfMCTDHF1D)
    m = psi.m
    @assert m.f==2
    m2 = psi2.m
    u2 = get_data(psi2, true)
    u2[:,:] = 0.0
    f = 1/sqrt(2)
    to_real_space!(psi)
    for j = 1:m.lena # eval slater determinants
        J = m.slater_indices[j]
        s = sign(J)
        v1 = get_data(psi.phi[J[1]], true)
        v2 = get_data(psi.phi[J[2]], true)
        #u2[:,:] += (s*f*psi.a[j])*(v1*v2.' - v2*v1.')
        for i1=1:size(u2,1)
            for i2=1:size(u2,2)
                u2[i1,i2] += (s*f*psi.a[j])*(v1[i1]*v2[i2]-v1[i2]*v2[i1])
            end
        end
    end    
    psi2
end


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
        v1 = get_data(psi.phi[J[1]], true)
        v2 = get_data(psi.phi[J[2]], true)
        v3 = get_data(psi.phi[J[3]], true)
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
            g[q] = inner_product(psi.phi[q], psi.phi[p])
        end
        for q=1:p-1
            axpy!(psi.phi[p], psi.phi[q], -g[q])
            for (j, l, s) in m.orthogonalization_rules[p,q]
                a1[l] += s*g[q]*psi.a[j]
            end
        end
        f = TSSM.norm(psi.phi[p])
        scale!(psi.phi[p],1/f)
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
            h = potential_matrix_element(psi.phi[p], psi.phi[q])
            for (j,l,f) in m.slater1_rules[q,p]
                V += h*f*conj(psi.a[j])*psi.a[l]
                #V += h*f*psi.a[j]*conj(psi.a[l])
            end
        end
    end
    real(V)
end


function potential_energy_1_A(psi::WfMCTDHF1D)
    m = psi.m
    V = 0
    for p=1:m.N    
        h = potential_matrix_element(psi.phi[p], psi.phi[p])
        for (j,l,f) in m.slater1_rules[p,p]
           #V += real(h*f*conj(psi.a[j])*psi.a[l])
            V += real(h*f*psi.a[j]*conj(psi.a[l]))
        end
        for q=1:p-1
            h = potential_matrix_element(psi.phi[p], psi.phi[q])
            for (j,l,f) in m.slater1_rules[p,q]
                #V += 2*real(h*f*conj(psi.a[j])*psi.a[l])
                V += 2*real(h*f*psi.a[j]*conj(psi.a[l]))
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
            h = kinetic_matrix_element(psi.phi[p], psi.phi[q])
            for (j,l,f) in m.slater1_rules[q,p]
                T += h*f*conj(psi.a[j])*psi.a[l]
                #T += h*f*psi.a[j]*conj(psi.a[l])
            end
        end
    end
    real(T)    
end


function kinetic_energy_A(psi::WfMCTDHF1D)
    m = psi.m
    V = 0
    for p=1:m.N    
        h = kinetic_matrix_element(psi.phi[p], psi.phi[p])
        for (j,l,f) in m.slater1_rules[p,p]
           #V += real(h*f*conj(psi.a[j])*psi.a[l])
            V += real(h*f*psi.a[j]*conj(psi.a[l]))
        end
        for q=1:p-1
            h = kinetic_matrix_element(psi.phi[p], psi.phi[q])
            for (j,l,f) in m.slater1_rules[p,q]
                #V += 2*real(h*f*conj(psi.a[j])*psi.a[l])
                V += 2*real(h*f*psi.a[j]*conj(psi.a[l]))
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
            u_pq[:] = m.Vee * (conj(get_data(psi.phi[p], true)).*get_data(psi.phi[q], true))
            for s=1:m.N
                u_pqs[:] = u_pq .* get_data(psi.phi[s], true)
                for r=1:m.N
                    h = dot(get_data(psi.phi[r], true), u_pqs)
                    for (j,l,f) in m.slater2_rules[q,p,s,r]
                        V += h*f*conj(psi.a[j])*psi.a[l]
                        #V += h*f*psi.a[j]*conj(psi.a[l])
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
        imaginary_time_propagate_A!(psi.phi[j], dt)
    end
end

function TSSM.imaginary_time_propagate_B!(psi::WfMCTDHF1D, dt::Real)
    for j=1:psi.m.N
        imaginary_time_propagate_B!(psi.phi[j], dt)
    end
end

function TSSM.scale!(psi::WfMCTDHF1D, f::Number)
    for j=1:psi.m.N
        scale!(psi.phi[j], f)
    end
    psi.a[:] *= f
end

function TSSM.axpy!(psi1::WfMCTDHF1D, psi2::WfMCTDHF1D, f::Number)
    for j=1:psi.m.N
        axpy!(psi1.phi[j], psi2.phi[j], f)
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
    