
type GeneralizedLaguerre2D <: TSSM
    m::Ptr{Void}
    function GeneralizedLaguerre2D(M::Integer, K::Integer, gamma_r::Real, Omega::Real) 
        meth = new( ccall( dlsym(tssm_handle, "c_new_gen_laguerre_2d"), Ptr{Void}, 
                   (Int32, Int32, Float64, Float64), M, K, gamma_r, Omega) )
        finalizer(meth, x -> ccall( dlsym(tssm_handle, "c_finalize_gen_laguerre_2d"), 
                  Void, (Ptr{Ptr{Void}},), &x.m) )
        meth           
    end
end

type GeneralizedLaguerreHermite3D <: TSSM
    m::Ptr{Void}
    function GeneralizedLaguerreHermite3D(M::Integer, K::Integer, gamma_r::Real, Omega::Real, nz::Integer, gamma_z::Real) 
        meth = new( ccall( dlsym(tssm_handle, "c_new_gen_laguerre_hermite_3d"), Ptr{Void}, 
                   (Int32, Int32, Float64, Float64, Int32, Float64), M, K, gamma_r, Omega, nz, gamma_z) )
        finalizer(meth, x -> ccall( dlsym(tssm_handle, "c_finalize_gen_laguerre_hermite_3d"), 
                  Void, (Ptr{Ptr{Void}},), &x.m) )
        meth           
    end
end

type GeneralizedLaguerreReal2D <: TSSM
    m::Ptr{Void}
    function GeneralizedLaguerreReal2D(M::Integer, K::Integer, gamma_r::Real, Omega::Real) 
        meth = new( ccall( dlsym(tssm_handle, "c_new_gen_laguerre_real_2d"), Ptr{Void}, 
                   (Int32, Int32, Float64, Float64), M, K, gamma_r, Omega) )
        finalizer(meth, x -> ccall( dlsym(tssm_handle, "c_finalize_gen_laguerre_real_2d"), 
                  Void, (Ptr{Ptr{Void}},), &x.m) )
        meth           
    end
end

type GeneralizedLaguerreHermiteReal3D <: TSSM
    m::Ptr{Void}
    function GeneralizedLaguerreHermiteReal3D(M::Integer, K::Integer, gamma_r::Real, Omega::Real, nz::Integer, gamma_z::Real) 
        meth = new( ccall( dlsym(tssm_handle, "c_new_gen_laguerre_hermite_real_3d"), Ptr{Void}, 
                   (Int32, Int32, Float64, Float64, Int32, Float64), M, K, gamma_r, Omega, nz, gamma_z) )
        finalizer(meth, x -> ccall( dlsym(tssm_handle, "c_finalize_gen_laguerre_hermite_real_3d"), 
                  Void, (Ptr{Ptr{Void}},), &x.m) )
        meth           
    end
end

## method: get_eigenvalues #######################################################################

function get_eigenvalues(m::GeneralizedLaguerre2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::GeneralizedLaguerreHermite3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2, ev3
   else
       return copy(ev1), copy(ev2), copy(ev3)
   end
end

function get_eigenvalues(m::GeneralizedLaguerreReal2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::GeneralizedLaguerreHermiteReal3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2, ev3
   else
       return copy(ev1), copy(ev2), copy(ev3)
   end
end

## method: get_nodes #######################################################################

function get_nodes(m::GeneralizedLaguerre2D)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   return copy(ev1), copy(ev2)
end

function get_nodes(m::GeneralizedLaguerreHermite3D)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp2, dim[1], false)     
   return copy(ev1), copy(ev2), copy(ev3)
end

function get_nodes(m::GeneralizedLaguerreReal2D)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   return copy(ev1), copy(ev2)
end

function get_nodes(m::GeneralizedLaguerreHermiteReal3D)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_nodes_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp2, dim[1], false)     
   return copy(ev1), copy(ev2), copy(ev3)
end

## method: get_weights #######################################################################

function get_weights(m::GeneralizedLaguerre2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_weights_gen_laguerre_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   if unsafe_access
       return ev1
   else
       return copy(ev1)
   end
end

function get_weights(m::GeneralizedLaguerreHermite3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_weights_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_weights_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_weights(m::GeneralizedLaguerreReal2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_weights_gen_laguerre_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   if unsafe_access
       return ev1
   else
       return copy(ev1)
   end
end

function get_weights(m::GeneralizedLaguerreHermiteReal3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_weights_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_weights_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

## method: get_L #######################################################################

function get_L(m::GeneralizedLaguerre2D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Lp = ccall( dlsym(tssm_handle, "c_get_l_gen_laguerre_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   L = pointer_to_array(Lp, (dims[1], dims[2], dims[3]), false)  
   if unsafe_access
       return L
   else
       return copy(L)
   end
end

function get_L(m::GeneralizedLaguerreHermite3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Lp = ccall( dlsym(tssm_handle, "c_get_l_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   L = pointer_to_array(Lp, (dims[1], dims[2], dims[3]), false)
   if unsafe_access
       return L
   else
       return copy(L)
   end
end

function get_L(m::GeneralizedLaguerreReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Lp = ccall( dlsym(tssm_handle, "c_get_l_gen_laguerre_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   L = pointer_to_array(Hp, (dims[1], dims[2], dims[3]), false)  
   if unsafe_access
       return L
   else
       return copy(L)
   end
end

function get_L(m::GeneralizedLaguerreHermiteReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Lp = ccall( dlsym(tssm_handle, "c_get_l_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   L = pointer_to_array(Hp, (dims[1], dims[2], dims[3]), false)
   if unsafe_access
       return L
   else
       return copy(L)
   end
end

## method: get_H #######################################################################

function get_H(m::GeneralizedLaguerreHermite3D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp = ccall( dlsym(tssm_handle, "c_get_h_gen_laguerre_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   H = pointer_to_array(Hp, (dims[1], dims[2]), false)  
   if unsafe_access
       return H
   else
       return copy(H)
   end
end

function get_H(m::GeneralizedLaguerreHermiteReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp = ccall( dlsym(tssm_handle, "c_get_h_gen_laguerre_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   H = pointer_to_array(Hp, (dims[1], dims[2]), false)  
   if unsafe_access
       return H
   else
       return copy(H)
   end
end


########################################################################################################
########################################################################################################


## Wave function types ############################################################################

type WfGeneralizedLaguerre2D <: WaveFunction2D
    p::Ptr{Void}
    m::GeneralizedLaguerre2D
    function WfGeneralizedLaguerre2D( m::GeneralizedLaguerre2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_gen_laguerre_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_gen_laguerre_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfGeneralizedLaguerreHermite3D <: WaveFunction3D
    p::Ptr{Void}
    m::GeneralizedLaguerreHermite3D
    function WfGeneralizedLaguerreHermite3D( m::GeneralizedLaguerreHermite3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_gen_laguerre_hermite_3d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_gen_laguerre_hermite_3d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfGeneralizedLaguerreReal2D <: WaveFunction2D
    p::Ptr{Void}
    m::GeneralizedLaguerreReal2D
    function WfGeneralizedLaguerreReal2D( m::GeneralizedLaguerreReal2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_gen_laguerre_real_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_gen_laguerre_real_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfGeneralizedLaguerreHermiteReal3D <: WaveFunction3D
    p::Ptr{Void}
    m::GeneralizedLaguerreHermiteReal3D
    function WfGeneralizedLaguerreHermiteReal3D( m::GeneralizedLaguerreHermiteReal3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_gen_laguerre_hermite_real_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_gen_laguerre_hermite_real_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

## method: wave_function ##########################################################################

function wave_function(m::GeneralizedLaguerre2D )
    WfGeneralizedLaguerre2D(m) 
end

function wave_function(m::GeneralizedLaguerreHermite3D )
    WfGeneralizedLaguerreHermite3D(m) 
end

function wave_function(m::GeneralizedLaguerreReal2D )
    WfGeneralizedLaguerreReal2D(m) 
end

function wave_function(m::GeneralizedLaguerreHermiteReal3D )
    WfGeneralizedLaguerreHermiteReal3D(m) 
end

## method: is_real_space###########################################################################

function is_real_space(psi::WfGeneralizedLaguerre2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfGeneralizedLaguerreHermite3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_hermite_3d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfGeneralizedLaguerreReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_real_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfGeneralizedLaguerreHermiteReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_hermite_real_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

## method: is_frequency_space###########################################################################

function is_frequency_space(psi::WfGeneralizedLaguerre2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfGeneralizedLaguerreHermite3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_hermite_3d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfGeneralizedLaguerreReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_real_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfGeneralizedLaguerreHermiteReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_gen_laguerre_hermite_real_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end


## method: to_real_space! ##########################################################################

function to_real_space!(psi::WfGeneralizedLaguerre2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_gen_laguerre_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfGeneralizedLaguerreHermite3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_gen_laguerre_hermite_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfGeneralizedLaguerreReal2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_gen_laguerre_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfGeneralizedLaguerreHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_gen_laguerre_hermite_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: to_frequency_space! ######################################################################

function to_frequency_space!(psi::WfGeneralizedLaguerre2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_gen_laguerre_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfGeneralizedLaguerreHermite3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_gen_laguerre_hermite_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfGeneralizedLaguerreReal2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_gen_laguerre_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfGeneralizedLaguerreHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_gen_laguerre_hermite_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: propagate_A! ############################################################################

function propagate_A!(psi::WfGeneralizedLaguerre2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_gen_laguerre_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfGeneralizedLaguerreHermite3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_gen_laguerre_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfGeneralizedLaguerreReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_gen_laguerre_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_A!(psi::WfGeneralizedLaguerreHermiteReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_gen_laguerre_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_B! ############################################################################

function propagate_B!(psi::WfGeneralizedLaguerre2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_gen_laguerre_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfGeneralizedLaguerreHermite3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_gen_laguerre_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfGeneralizedLaguerreReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_gen_laguerre_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_B!(psi::WfGeneralizedLaguerreHermiteReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_gen_laguerre_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_C! ############################################################################

function propagate_C!(psi::WfGeneralizedLaguerre2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_gen_laguerre_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_C!(psi::WfGeneralizedLaguerreHermite3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_gen_laguerre_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_C!(psi::WfGeneralizedLaguerreReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_gen_laguerre_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_C!(psi::WfGeneralizedLaguerreHermiteReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_gen_laguerre_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end
## method: add_apply_A! ############################################################################

function add_apply_A!(this::WfGeneralizedLaguerre2D, other::WfGeneralizedLaguerre2D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_gen_laguerre_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfGeneralizedLaguerreHermite3D, other::WfGeneralizedLaguerreHermite3D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_gen_laguerre_hermite_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfGeneralizedLaguerreReal2D, other::WfGeneralizedLaguerreReal2D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_gen_laguerre_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfGeneralizedLaguerreHermiteReal3D, other::WfGeneralizedLaguerreHermiteReal3D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_gen_laguerre_hermite_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

## method: save ###################################################################################

function save(psi::WfGeneralizedLaguerre2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_gen_laguerre_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfGeneralizedLaguerreHermite3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_gen_laguerre_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfGeneralizedLaguerreReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_gen_laguerre_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfGeneralizedLaguerreHermiteReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_gen_laguerre_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: load! ###################################################################################

function load!(psi::WfGeneralizedLaguerre2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_gen_laguerre_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfGeneralizedLaguerreHermite3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_gen_laguerre_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfGeneralizedLaguerreReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_gen_laguerre_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfGeneralizedLaguerreHermiteReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_gen_laguerre_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: norm ###################################################################################

function norm(psi::WfGeneralizedLaguerre2D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_gen_laguerre_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfGeneralizedLaguerreHermite3D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_gen_laguerre_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfGeneralizedLaguerreReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_gen_laguerre_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfGeneralizedLaguerreHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_gen_laguerre_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: inner_product ##########################################################################

function inner_product(psi1::WfGeneralizedLaguerre2D, psi2::WfGeneralizedLaguerre2D)
   ccall( dlsym(tssm_handle, "c_inner_product_wf_gen_laguerre_2d"), Complex128,
         (Ptr{Void}, Ptr{Void} ), psi1.p, psi2.p )
end

function inner_product(psi1::WfGeneralizedLaguerreReal2D, psi2::WfGeneralizedLaguerreReal2D)
   ccall( dlsym(tssm_handle, "c_inner_product_wf_gen_laguerre_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi1.p, psi2.p )
end

function inner_product(psi1::WfGeneralizedLaguerreHermite3D, psi2::WfGeneralizedLaguerreHermite3D)
   ccall( dlsym(tssm_handle, "c_inner_product_wf_gen_laguerre_hermite_3d"), Complex128,
         (Ptr{Void}, Ptr{Void} ), psi1.p, psi2.p )
end

function inner_product(psi1::WfGeneralizedLaguerreHermiteReal3D, psi2::WfGeneralizedLaguerreHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_inner_product_wf_gen_laguerre_hermite_real_3d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi1.p, psi2.p )
end



## method: norm_in_frequency_space ################################################################

function norm_in_frequency_space(psi::WfGeneralizedLaguerre2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_gen_laguerre_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfGeneralizedLaguerreHermite3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_gen_laguerre_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfGeneralizedLaguerreReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_gen_laguerre_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfGeneralizedLaguerreHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_gen_laguerre_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end


## method: distance ##############################################################################

function distance(psi1::WfGeneralizedLaguerre2D, psi2::WfGeneralizedLaguerre2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_gen_laguerre_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfGeneralizedLaguerreHermite3D, psi2::WfGeneralizedLaguerreHermite3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_gen_laguerre_hermite_3d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfGeneralizedLaguerreReal2D, psi2::WfGeneralizedLaguerreReal2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_gen_laguerre_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfGeneralizedLaguerreHermiteReal3D, psi2::WfGeneralizedLaguerreHermiteReal3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_gen_laguerre_hermite_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

## method: normalize! ##############################################################################

function normalize!(psi::WfGeneralizedLaguerre2D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_gen_laguerre_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfGeneralizedLaguerreHermite3D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_gen_laguerre_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfGeneralizedLaguerreReal2D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_gen_laguerre_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfGeneralizedLaguerreHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_gen_laguerre_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: scale! ##################################################################################

function scale!(psi::WfGeneralizedLaguerre2D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wf_gen_laguerre_2d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfGeneralizedLaguerreHermite3D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wf_gen_laguerre_hermite_3d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfGeneralizedLaguerreReal2D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wf_gen_laguerre_real_2d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

function scale!(psi::WfGeneralizedLaguerreHermiteReal3D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wf_gen_laguerre_hermite_real_2d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

## method: axpy! ###################################################################################

function axpy!(this::WfGeneralizedLaguerre2D, other::WfGeneralizedLaguerre2D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_gen_laguerre_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfGeneralizedLaguerreHermite3D, other::WfGeneralizedLaguerreHermite3D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_gen_laguerre_hermite_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfGeneralizedLaguerreReal2D, other::WfGeneralizedLaguerreReal2D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_gen_laguerre_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

function axpy!(this::WfGeneralizedLaguerreHermiteReal3D, other::WfGeneralizedLaguerreHermiteReal3D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_gen_laguerre_hermite_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

## method: get_data ###########################################################################

function get_data(psi::WfGeneralizedLaguerre2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_gen_laguerre_2d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfGeneralizedLaguerreHermite3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_gen_laguerre_hermite_3d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfGeneralizedLaguerreReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   if is_real_space(psi)
      up = ccall( dlsym(tssm_handle, "c_get_data_wf_gen_laguerre_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
      data = pointer_to_array(up, (dims[1], dims[2]), false)     
   else
      up = ccall( dlsym(tssm_handle, "c_get_data_wf_gen_laguerre_real_2d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
      data = pointer_to_array(up, (dims[1], div(dims[2],2)), false)     
   end
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfGeneralizedLaguerreHermiteReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   if is_real_space(psi)
      up = ccall( dlsym(tssm_handle, "c_get_data_wf_gen_laguerre_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
      data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
   else
      up = ccall( dlsym(tssm_handle, "c_get_data_wf_gen_laguerre_hermite_real_2d"),Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
      data = pointer_to_array(up, (dims[1], dims[2], div(dims[3],2)), false)     
   end
   if unsafe_access
      return data
   else
      return copy(data)
   end
end
#TODO: acces_data_frequency_space for *real* spectral methods

## method: set! ###################################################################################

function set!(psi::WfGeneralizedLaguerre2D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wf_gen_laguerre_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wf_gen_laguerre_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end      
end

function set!(psi::WfGeneralizedLaguerreHermite3D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wf_gen_laguerre_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wf_gen_laguerre_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end      
end

function set!(psi::WfGeneralizedLaguerreReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wf_gen_laguerre_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

function set!(psi::WfGeneralizedLaguerreHermiteReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wf_gen_laguerre_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

## method: set! with t argument ###########################################################################

function set!(psi::WfGeneralizedLaguerre2D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wf_gen_laguerre_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wf_gen_laguerre_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfGeneralizedLaguerreHermite3D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wf_gen_laguerre_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wf_gen_laguerre_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfGeneralizedLaguerreReal2D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wf_gen_laguerre_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

function set!(psi::WfGeneralizedLaguerreHermiteReal3D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wf_gen_laguerre_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

## method: copy! ##################################################################################

function copy!(target::WfGeneralizedLaguerre2D, source::WfGeneralizedLaguerre2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_gen_laguerre_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfGeneralizedLaguerreHermite3D, source::WfGeneralizedLaguerreHermite3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_gen_laguerre_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfGeneralizedLaguerreReal2D, source::WfGeneralizedLaguerreReal2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_gen_laguerre_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfGeneralizedLaguerreHermiteReal3D, source::WfGeneralizedLaguerreHermiteReal3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_gen_laguerre_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end


