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
    slater_rules = [(Int[], 1) for j=1:lena, l=1:lena]
    # entry = (exchange_index_pairs, sign)
    for j=1:lena
        for l=1:lena
            u = symdiff(slater_indices[j], slater_indices[l])
            if length(u)==2
                i1 = findfirst(slater_indices[j], u[1])
                i2 = findfirst(slater_indices[l], u[2])
                s = (-1)^(i1+i2)
                slater_rules[j,l] = (u, s)
            elseif length(u)==4
                u = [u[1], u[3], u[2], u[4]]
                i1 = findfirst(slater_indices[j], u[1])
                i2 = findfirst(slater_indices[l], u[2])
                i3 = findfirst(slater_indices[j], u[3])
                i4 = findfirst(slater_indices[l], u[4])
                s = (-1)^(i1+i2+i3+i4)
                slater_rules[j,l] = (u, s)
            end
        end
    end
    
    slater_indices, density_rules, density2_rules, slater_rules
end


type MCTDHF1D <: TSSM.TimeSplittingSpectralMethodComplex1D
    m::Schroedinger1D
    f::Int # number of electrond
    N::Int # number of orbitals
    lena::Int # number of (independent) coefficients
    slater_indices
    density_rules 
    density2_rules 
    slater_rules

    function MCTDHF1D(f::Integer, N::Integer, 
                      nx::Integer, xmin::Real, xmax::Real)
        m = Schroedinger1D(nx, xmin, xmax)
        lena = binomial(N,f)
        slater_indices, density_rules, density2_rules, slater_rules = 
              init_mctdhf_combinatorics(f, N)
        new(m, f, N, lena, slater_indices, density_rules, 
                           density2_rules, slater_rules)
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


function density_matrix(psi::WfMCTDHF1D)
    N = psi.m.N
    rho = zeros(N, N)
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
    rho
end


function density2_tensor(psi::WfMCTDHF1D)
    N = psi.m.N
    rho = zeros(N, N, N, N)
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
    rho
end

