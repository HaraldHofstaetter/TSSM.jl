# to be includes by tssm.jl

#quadrature formulas:

const gauss = 1
const radau = 2
const lobatto = 3


type FourierBessel2D <: TSSM
    m::Ptr{Void}
    function FourierBessel2D(M::Integer, nr::Integer, nfr::Integer; 
                             r_max::Real = 1.0,
                             boundary_conditions::Integer = dirichlet,
                             quadrature_rule::Integer = (boundary_conditions==neumann ? radau : lobatto)
                            ) 
        meth = new( ccall( dlsym(tssm_handle, "c_new_fourier_bessel_2d"), Ptr{Void}, 
                   (Int32, Int32, Int32, Float64, Int32, Int32), M, nr, nfr,  r_max, boundary_conditions, quadrature_rule) )
        finalizer(meth, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_bessel_2d"), 
                  Void, (Ptr{Ptr{Void}},), &x.m) )
        meth           
    end
end

type FourierBesselReal2D <: TSSM
    m::Ptr{Void}
    function FourierBesselReal2D(M::Integer, nr::Integer, nfr::Integer;
                                 r_max::Real = 1.0,
                                 boundary_conditions::Integer = dirichlet,
                                 quadrature_rule::Integer = (boundary_conditions==neumann ? radau : lobatto)
                                )
        meth = new( ccall( dlsym(tssm_handle, "c_new_fourier_bessel_real_2d"), Ptr{Void}, 
                   (Int32, Int32, Int32, Float64, Int32, Int32), M, nr, nfr, r_max, boundary_conditions, quadrature_rule) )
        finalizer(meth, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_bessel_real_2d"), 
                  Void, (Ptr{Ptr{Void}},), &x.m) )
        meth           
    end
end

## method: get_eigenvalues #######################################################################

function get_eigenvalues(m::FourierBessel2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_bessel_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   ev = pointer_to_array(evp, (dims[1], dims[2]), false)     
   if unsafe_access
       return ev
   else
       return copy(ev)
   end
end

function get_eigenvalues(m::FourierBesselReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_bessel_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   ev = pointer_to_array(evp, (dims[1], dims[2]), false)     
   if unsafe_access
       return ev
   else
       return copy(ev)
   end
end


## method: get_nodes #######################################################################

function get_nodes(m::FourierBessel2D)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_bessel_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_bessel_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   return copy(ev1), copy(ev2)
end

function get_nodes(m::FourierBesselReal2D)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_bessel_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_bessel_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   return copy(ev1), copy(ev2)
end

## method: get_weights #######################################################################

function get_weights(m::FourierBessel2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_weights_fourier_bessel_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   if unsafe_access
       return ev1
   else
       return copy(ev1)
   end
end

function get_weights(m::FourierBesselReal2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_weights_fourier_bessel_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   if unsafe_access
       return ev1
   else
       return copy(ev1)
   end
end

## method: get_L #######################################################################

function get_L(m::FourierBessel2D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Lp = ccall( dlsym(tssm_handle, "c_get_l_fourier_bessel_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   L = pointer_to_array(Lp, (dims[1], dims[2], dims[3]), false)  
   if unsafe_access
       return L
   else
       return copy(L)
   end
end

function get_L(m::FourierBesselReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Lp = ccall( dlsym(tssm_handle, "c_get_l_fourier_bessel_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   L = pointer_to_array(Lp, (dims[1], dims[2], dims[3]), false)  
   if unsafe_access
       return L
   else
       return copy(L)
   end
end


########################################################################################################
########################################################################################################


## Wave function types ############################################################################

type WfFourierBessel2D <: WaveFunction2D
    p::Ptr{Void}
    m::FourierBessel2D
    function WfFourierBessel2D( m::FourierBessel2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_bessel_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_bessel_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfFourierBesselReal2D <: WaveFunction2D
    p::Ptr{Void}
    m::FourierBesselReal2D
    function WfFourierBesselReal2D( m::FourierBesselReal2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_bessel_real_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_bessel_real_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

## method: wave_function ##########################################################################

function wave_function(m::FourierBessel2D )
    WfFourierBessel2D(m) 
end

function wave_function(m::FourierBesselReal2D )
    WfFourierBesselReal2D(m) 
end

## method: is_real_space###########################################################################

function is_real_space(psi::WfFourierBessel2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_bessel_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfFourierBesselReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_bessel_real_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end


## method: is_frequency_space###########################################################################

function is_frequency_space(psi::WfFourierBessel2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_bessel_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfFourierBesselReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_bessel_real_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end


## method: to_real_space! ##########################################################################

function to_real_space!(psi::WfFourierBessel2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_bessel_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfFourierBesselReal2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_bessel_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: to_frequency_space! ######################################################################

function to_frequency_space!(psi::WfFourierBessel2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_bessel_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfFourierBesselReal2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_bessel_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: propagate_A! ############################################################################

function propagate_A!(psi::WfFourierBessel2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_bessel_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfFourierBesselReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_bessel_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_B! ############################################################################

function propagate_B!(psi::WfFourierBessel2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_bessel_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfFourierBesselReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_bessel_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_C! ############################################################################

function propagate_C!(psi::WfFourierBessel2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_bessel_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_C!(psi::WfFourierBesselReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_bessel_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: add_apply_A! ############################################################################

function add_apply_A!(this::WfFourierBessel2D, other::WfFourierBessel2D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_bessel_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfFourierBesselReal2D, other::WfFourierBesselReal2D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_bessel_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

## method: save ###################################################################################

function save(psi::WfFourierBessel2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_bessel_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfFourierBesselReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_bessel_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: load! ###################################################################################

function load!(psi::WfFourierBessel2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_bessel_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfFourierBesselReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_bessel_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: eval ###################################################################################

function eval(psi::WfFourierBessel2D, x::Real, y::Real)
   ccall( dlsym(tssm_handle, "c_eval_wf_fourier_bessel_2d"), Complex128,
         (Ptr{Void}, Float64, Float64), psi.p, x, y )
end

function eval(psi::WfFourierBesselReal2D, x::Real, y::Real)
   ccall( dlsym(tssm_handle, "c_eval_wf_fourier_bessel_real_2d"), Float64,
         (Ptr{Void}, Float64, Float64), psi.p, x, y )
end


## method: norm ###################################################################################

function norm(psi::WfFourierBessel2D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_bessel_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfFourierBesselReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_bessel_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: inner_product ##########################################################################

function inner_product(psi1::WfFourierBessel2D, psi2::WfFourierBessel2D)
   ccall( dlsym(tssm_handle, "c_inner_product_wf_fourier_bessel_2d"), Complex128,
         (Ptr{Void}, Ptr{Void} ), psi1.p, psi2.p )
end

function inner_product(psi1::WfFourierBesselReal2D, psi2::WfFourierBesselReal2D)
   ccall( dlsym(tssm_handle, "c_inner_product_wf_fourier_bessel_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi1.p, psi2.p )
end

## method: norm_in_frequency_space ################################################################

function norm_in_frequency_space(psi::WfFourierBessel2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_fourier_bessel_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfFourierBesselReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_fourier_bessel_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end


## method: distance ##############################################################################

function distance(psi1::WfFourierBessel2D, psi2::WfFourierBessel2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_bessel_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfFourierBesselReal2D, psi2::WfFourierBesselReal2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_bessel_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

## method: normalize! ##############################################################################

function normalize!(psi::WfFourierBessel2D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_bessel_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfFourierBesselReal2D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_bessel_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: scale! ##################################################################################

function scale!(psi::WfFourierBessel2D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_bessel_2d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfFourierBesselReal2D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_bessel_real_2d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

## method: axpy! ###################################################################################

function axpy!(this::WfFourierBessel2D, other::WfFourierBessel2D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_bessel_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfFourierBesselReal2D, other::WfFourierBesselReal2D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_bessel_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

## method: get_data ###########################################################################

function get_data(psi::WfFourierBessel2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_bessel_2d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfFourierBesselReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   if is_real_space(psi)
      up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_bessel_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
      data = pointer_to_array(up, (dims[1], dims[2]), false)     
   else
      up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_bessel_real_2d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
      data = pointer_to_array(up, (dims[1], div(dims[2],2)), false)     
   end
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

#TODO: acces_data_frequency_space for *real* spectral methods

## method: set! ###################################################################################

function set!(psi::WfFourierBessel2D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wf_fourier_bessel_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wf_fourier_bessel_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end      
end

function set!(psi::WfFourierBesselReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wf_fourier_bessel_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

## method: set! with t argument ###########################################################################

function set!(psi::WfFourierBessel2D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wf_fourier_bessel_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_bessel_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfFourierBesselReal2D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_bessel_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

## method: copy! ##################################################################################

function copy!(target::WfFourierBessel2D, source::WfFourierBessel2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_bessel_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfFourierBesselReal2D, source::WfFourierBesselReal2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_bessel_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end



