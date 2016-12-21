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
                   add_xxx!(res[p,q], (j,j, 1.0))
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
                            add_xxx!(res[p,q,r,s], (j,j, 0.5))
                        end
                    end
                    if p==s && r==q && p!=r
                        for j=1:lena
                            add_xxx!(res[p,q,r,s], (j,j, -0.5)) 
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
    
    slater_indices, density_rules, density2_rules, slater_exchange, slater1_rules, slater2_rules
end


function init_Vee(x::Vector{Float64})
    n = length(x)
    1./(sqrt(kron(x',ones(n)) - kron(x,ones(1,n)))+1)    
end


type MCTDHF1D <: TSSM.TimeSplittingSpectralMethodComplex1D
    m::Schroedinger1D
    f::Int # number of electrond
    N::Int # number of orbitals
    lena::Int # number of (independent) coefficients
    
    slater_indices
    density_rules 
    density2_rules
    slater_exchange
    slater1_rules
    slater2_rules

    Vee
    density_matrix
    density2_tensor

    function MCTDHF1D(f::Integer, N::Integer, 
                      nx::Integer, xmin::Real, xmax::Real)
        m = Schroedinger1D(nx, xmin, xmax)
        lena = binomial(N,f)
        slater_indices, density_rules, density2_rules, slater_exchange, slater1_rules, slater2_rules = 
              init_mctdhf_combinatorics(f, N)
        Vee = init_Vee(get_nodes(m))
        new(m, f, N, lena, slater_indices, density_rules, 
            density2_rules, slater_exchange, slater1_rules, slater2_rules, Vee,
            zeros(Complex{Float64},N,N),  zeros(Complex{Float64},N,N,N,N) )
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
            for q=q:N
                rho[:,:,p,q] = psi.m.density_matrix \ rho[:,:,p,q]
            end
        end
    end
    nothing
end


function norm(psi::WfMCTDHF1D)
    m = psi.m
    res = 0
    for p=1:m.N
        for q=1:m.N
            h = inner_product(psi.phi[p], psi.phi[q])
            for (j,l,f) in m.slater1_rules[p,q]
                res += h*f*conj(psi.a[j])*psi.a[l]
            end
        end
    end
    res    
end


function to_real_space!(psi::WfMCTDHF1D)
    for j=1:psi.m.N
        to_real_space!(psi.phi[j])
    end
end


function to_fourier_space!(psi::WfMCTDHF1D)
    for j=1:psi.m.N
        to_fourier_space!(psi.phi[j])
    end
end


function set_zero!(psi::WfMCTDHF1D)
    for j=1:psi.m.N 
        get_data(psi.phi[j], unsafe_access=true)[:] = 0.0
    end
    psi.a[:] = 0.0
end



function gen_rhs1(rhs::WfMCTDHF1D, psi::WfMCTDHF1D; A::Bool=true, b::Bool=true)
    m = rhs.m
    if m ≠ psi.m
        error("rhs and psi must belong to the same method")
    end
    n = size(m.Vee, 1)
    u_save = zeros(n)    
    for q=1:m.N
        to_real_space!(rhs.phi[q])
        u_save[:] = get_data(rhs.phi[q], unsafe_access=false)
        u = get_data(rhs.phi[q], unsafe_access=true)  
        u[:] = 0.0
        if A
            add_apply_A!(rhs.phi[q], psi.phi[q])
        end
        if B
            add_apply_B!(rhs.phi[q], psi.phi[q])
        end    
        to_real_space!(rhs.phi[q])
        for p=1:m.N
            h = inner_product(psi.phi[p], rhs.phi[q])
            for (j,l,f) in m.slater1_rules[p,q]
                rhs.a[j] += h*f*psi.a[l] 
            end
        end
        u += u_save        
    end
end



function gen_rhs2!(rhs::WfMCTDHF1D, psi::WfMCTDHF1D)
    m = rhs.m
    if m ≠ psi.m
        error("rhs and psi must belong to the same method")
    end
    n = size(m.Vee, 1)
    u_pq = zeros(n)
    u_pqs = zeros(n)
    to_real_space!(psi)
    to_real_space!(rhs)
    for p=1:m.N
        for q=1:m.N
            u_pq[:] = m.Vee * (conj(get_data(psi.phi[p], unsafe_access=true)).*get_data(psi.phi[q], unsafe_access=true))
            for s=1:m.N
                u_pqs[:] = u_pq .* get_data(psi.phi[s], unsafe_access=true)
                for r=1:m.N
                    u = get_data(rhs.phi[r], unsafe_access=true)
                    u += (m.density2_tensor[s,r,p,q]*(m.f-1)) * u_pqs                
                    h = get_data(psi.phi[r], unsafe_access=true)' * u_pqs # inner product...
                    # maybe factor 1/f! or something similar necessary...
                    for (j,l,f) in m.slater2_rules[p,q,r,s]
                        rhs.a[j] += h*f*psi.a[l] 
                    end
                end
            end
        end
    end
    #TODO: projection (1-P) ...
end

