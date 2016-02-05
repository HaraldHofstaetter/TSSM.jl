# T , TSSM_HANDLE 

println("including tssm_fourier_bessel.jl for type ",T)
for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:FourierBessel2D, :_fourier_bessel_2d, true, 2 ),             
                 (:FourierBesselReal2D, :_fourier_bessel_real_2d, false, 2 ),
                )
println("    ", METHOD)    

@eval begin
    function ($METHOD)(T::Type{$T}, ntheta::Integer, nr::Integer, nfr::Integer;
                       rmax::Real=1.0, #one($T),
                       boundary_conditions::Integer = dirichlet,
                       quadrature_rule::Integer = 
                       (boundary_conditions==neumann ? radau : lobatto))
        ccall( Libdl.dlsym(($TSSM_HANDLE), "c_initialize_tssm_fourier"), Void, ())
        m = ($METHOD){$T}( ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new",SUF))), 
                       Ptr{Void}, (Int32, Int32, Int32, ($T), Int32, Int32), 
                       ntheta, nr, nfr, rmax, boundary_conditions, quadrature_rule))
        finalizer(m, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), 
                       $(string("c_finalize",SUF))), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
end # eval
if T == :Float64
@eval begin
    function ($METHOD)(ntheta::Integer, nr::Integer, nfr::Integer;
                       rmax::Real=1.0, #one($T),
                       boundary_conditions::Integer = dirichlet,
                       quadrature_rule::Integer = 
                       (boundary_conditions==neumann ? radau : lobatto))
        ($METHOD)(($T), ntheta, nr, nfr, rmax=rmax,
                  boundary_conditions=boundary_conditions,
                  quadrature_rule=quadrature_rule)
    end      
end # eval
end # if

@eval begin
    function get_ntheta(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_ntheta",SUF))), Int32,
             (Ptr{Void}, ), m.m )
    end

    function get_nr(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nr",SUF))), Int32,
             (Ptr{Void}, ), m.m )
    end

    function get_nfr(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nfr",SUF))), Int32,
             (Ptr{Void}, ), m.m )
    end

    function get_rmax(m::($METHOD){$T})
       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_rmax",SUF))), ($T),
             (Ptr{Void}, ), m.m )
    end
    
    function get_weights(m::($METHOD){$T})
       dim =Array(Int32, 1)
       np = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_weights",SUF))), Ptr{$T},
             (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
       n = pointer_to_array(np, dim[1], false)     
       copy(n)
    end

    function get_L(m::($METHOD){$T}, unsafe_access::Bool=false)
       dims =Array(Int32, 3)
       Lp = ccall( dlsym(tssm_handle, $(string("c_get_l",SUF))), Ptr{$T},
             (Ptr{Void}, Ptr{Int32}), m.m, dims )
       L = pointer_to_array(Lp, (dims[1], dims[2], dims[3]), false)  
       if unsafe_access
           return L
       else
           return copy(L)
       end
    end

end # eval

end # for

