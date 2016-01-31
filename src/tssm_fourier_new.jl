module tssm_fourier_new

import float128: Float128, Complex256

using Base.Libdl.dlopen
using Base.Libdl.dlsym
using Base.Libdl.dlext

import Base.copy!
import Base.scale!
import Base.norm

export WaveFunction, WaveFunction1D, WaveFunction2D, WaveFunction3D, TSSM
export dim

export periodic, dirichlet, neumann

export Fourier1D, WfFourier1D 
export Fourier2D, WfFourier2D 
export Fourier3D, WfFourier3D 
export FourierReal1D, WfFourierReal1D 
export FourierReal2D, WfFourierReal2D 
export FourierReal3D, WfFourierReal3D 

function __init__()
    global tssm_handle = 0
    global tssmq_handle = 0
end

function lazy_init_tssm(T::DataType)
    if T==Float64
        if tssm_handle==0
            #tssm_handle = Libdl.dlopen( (@windows? :"libtssm.dll" : ( @osx? "libtssm.dylib" : :"libtssm.so" )) );
            tssm_handle = Libdl.dlopen( string("libtssm.", Libdl.dlext) );
            ccall( Libdl.dlsym(tssm_handle, "c_initialize_tssm"), Void, ())
        end
        return tssm_handle
    elseif T==Float128
        if tssmq_handle==0
            #tssmq_handle = Libdl.dlopen( (@windows? :"libtssmq.dll" : ( @osx? "libtssmq.dylib" : :"libtssmq.so" )) );
            tssmq_handle = Libdl.dlopen( string("libtssmq.", Libdl.dlext) );
            ccall( Libdl.dlsym(tssmq_handle, "c_initialize_tssm"), Void, ())
        end
        return tssmq_handle
    else
        error(string("tssm not implemented for ",T))
    end
end

abstract TSSM

abstract WaveFunction
abstract WaveFunction1D <: WaveFunction
abstract WaveFunction2D <: WaveFunction
abstract WaveFunction3D <: WaveFunction

dim(psi::WaveFunction1D) = 1
dim(psi::WaveFunction2D) = 2
dim(psi::WaveFunction3D) = 3

## Boundary conditions ############################################################################

const periodic = 0
const dirichlet = 1
const neumann = 2

none_1D(x)=0.0
none_2D(x,y)=0.0
none_3D(x,y,z)=0.0

function clone(psi::WaveFunction)
    wave_function(psi.m)
end



########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# Note: the follwing types cannot be declared as 'immutable' because
# otherwise finalizer won't work.

type Fourier1D{T} <: TSSM
    m::Ptr{Void}
    function Fourier1D(nx::Integer, xmin::Real, xmax::Real; 
                        boundary_conditions::Integer=periodic)
        h = lazy_init_tssm(T)
        ccall( Libdl.dlsym(h, "c_initialize_tssm_fourier"), Void, ())
        #m = new( ccall( Libdl.dlsym(h, "c_new_fourier_1d"), Ptr{Void}, 
        m = new( ccall( Libdl.dlsym(h, "c_new_fourier_1d"), Ptr{Void}, 
                   (Int32, T, T, Int32), 
                   nx, xmin, xmax, boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(h, "c_finalize_fourier_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
end

type Fourier2D{T} <: TSSM
    m::Ptr{Void}
    function Fourier2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                        boundary_conditions::Integer=periodic)
        h = lazy_init_tssm(T)
        ccall( Libdl.dlsym(h, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( Libdl.dlsym(h, "c_new_fourier_2d"), Ptr{Void}, 
                   (Int32, T, T, Int32, T, T, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(h, "c_finalize_fourier_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
end

type Fourier3D{T} <: TSSM
    m::Ptr{Void}
    function Fourier3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real; boundary_conditions::Integer=periodic)
        h = lazy_init_tssm(T)
        ccall( Libdl.dlsym(h, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( Libdl.dlsym(h, "c_new_fourier_3d"), Ptr{Void}, 
                   (Int32, T, T, Int32, T, T, Int32, T, T, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(h, "c_finalize_fourier_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
end

type FourierReal1D{T} <: TSSM
    m::Ptr{Void}
    function FourierReal1D(nx::Integer, xmin::Real, xmax::Real; 
                        boundary_conditions::Integer=periodic)
        h = lazy_init_tssm(T)
        ccall( Libdl.dlsym(h, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( Libdl.dlsym(h, "c_new_fourier_real_1d"), Ptr{Void}, 
                   (Int32, T, T, Int32), 
                   nx, xmin, xmax, boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(h, "c_finalize_fourier_real_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
end

type FourierReal2D{T} <: TSSM
    m::Ptr{Void}
    function FourierReal2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                        boundary_conditions::Integer=periodic)
        h = lazy_init_tssm(T)
        ccall( Libdl.dlsym(h, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( Libdl.dlsym(h, "c_new_fourier_real_2d"), Ptr{Void}, 
                   (Int32, T, T, Int32, T, T, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(h, "c_finalize_fourier_real_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
end

type FourierReal3D{T} <: TSSM
    m::Ptr{Void}
    function FourierReal3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real; boundary_conditions::Integer=periodic)
        h = lazy_init_tssm(T)
        ccall( Libdl.dlsym(h, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( Libdl.dlsym(h, "c_new_fourier_real_3d"), Ptr{Void}, 
                   (Int32, T, T, Int32, T, T, Int32, T, T, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, boundary_conditions))
        finalizer(m, x -> ccall( Libdl.dlsym(h, "c_finalize_fourier_real_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m
    end
end


#default type is Float64:
function Fourier1D(nx::Integer, xmin::Real, xmax::Real; 
                   boundary_conditions::Integer=periodic)
    Fourier1D{Float64}(nx, xmin, xmax,  boundary_conditions=boundary_conditions)
end

function Fourier2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                    boundary_conditions::Integer=periodic)
    Fourier2D{Float64}(nx, xmin, xmax, ny, ymin, ymax, boundary_conditions=boundary_conditions)
end

function Fourier3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                   nz::Integer, zmin::Real, zmax::Real; boundary_conditions::Integer=periodic)
    Fourier3D{Float64}(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,  boundary_conditions=boundary_conditions)
end

function FourierReal1D(nx::Integer, xmin::Real, xmax::Real; 
                    boundary_conditions::Integer=periodic)
    Fourier1DReal{Float64}(nx, xmin, xmax,  boundary_conditions=boundary_conditions)
end

function FourierReal2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                    boundary_conditions::Integer=periodic)
    Fourier2DReal{Float64}(nx, xmin, xmax, ny, ymin, ymax, boundary_conditions=boundary_conditions)
end

function FourierReal3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                   nz::Integer, zmin::Real, zmax::Real; boundary_conditions::Integer=periodic)
    Fourier3DReal{Float64}(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax,  boundary_conditions=boundary_conditions)
end



for (METHOD, WF, PARENT_WF, SUF, COMPLEX_METHOD, DIM) in (
                 (:Fourier1D, :WfFourier1D, :WaveFunction1D, :_fourier_1d, true, 1 ),             
                 (:Fourier2D, :WfFourier2D, :WaveFunction2D, :_fourier_2d, true, 2 ),             
                 (:Fourier3D, :WfFourier3D, :WaveFunction3D, :_fourier_3d, true, 3 ),             
                 (:FourierReal1D, :WfFourierReal1D, :WaveFunction1D, :_fourier_real_1d, false, 1 ), 
                 (:FourierReal2D, :WfFourierReal2D, :WaveFunction2D, :_fourier_real_2d, false, 2 ),
                 (:FourierReal3D, :WfFourierReal3D, :WaveFunction3D, :_fourier_real_3d, false, 3 ),
                )

    @eval begin

        type ($WF){T}  <: ($PARENT_WF)
            p::Ptr{Void}
            m::$METHOD{T}
            function ($WF){T}( m::($METHOD){T} )
                h = lazy_init_tssm(T)
                wf = new( ccall( Libdl.dlsym(h, $(string("c_new_wf",SUF))), Ptr{Void}, (Ptr{Void},), m.m) , m)   
                finalizer(wf, x -> ccall( Libdl.dlsym(h, $(string("c_finalize_wf",SUF))), Void, (Ptr{Ptr{Void}},), &x.p) )
                wf
            end
        end    
    end #eval    
   
    for (T, TSSM_HANDLE, LIBTSSM) in (
          (:Float64, :tssm_handle, "libtssm"),
          (:Float128, :tssmq_handle, "libtssmq"), 
          )
    @eval begin


        function wave_function(m::($METHOD){$T})
            ($WF){$T}(m) 
        end

        function is_real_space(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_is_real_space_wf",SUF))), Int32,
               (Ptr{Void},), psi.p) == 1
        end

        function is_frequency_space(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_is_real_space_wf",SUF))), Int32,
               (Ptr{Void},), psi.p) != 1
        end

        function to_real_space!(psi::($WF){$T})
             ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_to_real_space_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end

        function to_frequency_space_space!(psi::($WF){$T})
             ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_to_frequency_space_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end

        
        function save(psi::($WF){$T}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_save_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
        end         

        function load!(psi::($WF){$T}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_load_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
        end    

        function copy!(target::($WF){$T}, source::($WF){$T})
           if target.m ≠ source.m
               error("source and target must belong to the same method")
           end
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_copy_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{Void}), target.p, source.p )
        end
        
        function norm(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_norm2_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function norm_in_frequency_space(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_norm_in_frequency_space2_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function distance(psi1::($WF){$T}, psi2::($WF){$T})
           if psi1.m ≠ psi2.m
               error("psi1 and psi2 must belong to the same method")
           end
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_distance_wf",SUF))), ($T),
                 (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
        end

        function normalize!(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_normalize_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end


        function get_eigenvalues(m::($METHOD){$T}, unsafe_access::Bool=false)
           dim =Array(Int32, 1)
           evp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_eigenvalues",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
           ev = pointer_to_array(evp, dim[1], false)   
           if unsafe_access
               return ev
           else    
               return copy(ev)
           end
        end

    end #eval

    if DIM==1
        @eval begin

            function get_nodes(m::($METHOD){$T})
               println("get_nodes T=", string($T))
               dim =Array(Int32, 1)
               #np = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nodes",SUF))), Ptr{$T},
               np = ccall( ( ($(string("c_get_nodes",SUF))), ($LIBTSSM)), Ptr{($T)},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
               n = pointer_to_array(np, dim[1], false)     
               copy(n)
            end

        end  
    elseif DIM==2
        @eval begin

            function get_nodes(m::($METHOD){$T})
               dim =Array(Int32, 1)
               np1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
               n1 = pointer_to_array(n1p, dim[1], false)     
               np2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
               n2 = pointer_to_array(n2p, dim[1], false)     
               copy(n1), copy(n2)
            end

        end  
    elseif DIM==2
        @eval begin

            function get_nodes(m::($METHOD){$T})
               dim =Array(Int32, 1)
               np1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
               n1 = pointer_to_array(n1p, dim[1], false)     
               np2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
               n2 = pointer_to_array(n2p, dim[1], false)     
               np3 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
               n3 = pointer_to_array(n3p, dim[1], false)     
               copy(n1), copy(n2), copy(n3)
            end

          end  
     end


    @eval begin
        function get_nx(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nx",SUF))), Int32,
                 (Ptr{Void}, ), m.m )
        end

        function get_xmin(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_xmin",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end

        function get_xmax(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_xmax",SUF))), $T,
                 (Ptr{Void}, ), m.m )
        end
    end # eval    

    if DIM>=2        
        @eval begin
        
            function get_ny(m::($METHOD){$T})
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_ny",SUF))), Int32,
                     (Ptr{Void}, ), m.m )
            end

            function get_ymin(m::($METHOD){$T})
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_ymin",SUF))), ($T),
                     (Ptr{Void}, ), m.m )
            end

            function get_ymax(m::($METHOD){$T})
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_ymax",SUF))), ($T),
                     (Ptr{Void}, ), m.m )
            end
        end # eval    
    end

    if DIM>=3 
        @eval begin
            function get_nz(m::($METHOD){$T})
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nz",SUF))), Int32,
                     (Ptr{Void}, ), m.m )
            end

            function get_zmin(m::($METHOD){$T})
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_zmin",SUF))), ($T),
                     (Ptr{Void}, ), m.m )
            end

            function get_zmax(m::($METHOD){$T})
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_zmax",SUF))), ($T),
                     (Ptr{Void}, ), m.m )
            end
        end # eval    
    end 


    if COMPLEX_METHOD
        @eval begin

            function propagate_A!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_a_wf",SUF))), Void,
                        (Ptr{Void}, Complex{$T},), psi.p, dt)
            end
    
            function propagate_B!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_b_wf",SUF))), Void,
                        (Ptr{Void}, Complex{$T},), psi.p, dt)
            end
    
            function propagate_C!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_c_wf",SUF))), Void,
                        (Ptr{Void}, Complex{$T},), psi.p, dt)
            end
    
            function add_apply_A!(this::($WF){$T}, other::($WF){$T},
                                 coefficient::Number=1.0)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_add_apply_a_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, Complex{$T}), 
                         this.p, other.p, coefficient)
            end
    
            function scale!(psi::($WF){$T}, factor::Number)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_scale_wf",SUF))), Void,
                     (Ptr{Void}, Complex{$T} ), psi.p, factor )
            end

            function axpy!(this::($WF){$T}, other::($WF){$T},
                                 factor::Number)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_axpy_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, Complex{$T}), 
                         this.p, other.p, factor)
            end

            function set!(psi::($WF){$T}, f::Function)
               try
                   f_c = cfunction(f, T, (T,))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
               catch
                   f_c = cfunction(f, Complex{T}, (T,))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
               end      
            end
    
            function set!(psi::($WF){$T}, f::Function, t::Real)
               try
                   f_c = cfunction(f, T, (T, T,))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
               catch
                   f_c = cfunction(f, Complex{T}, (T, T,))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
               end      
            end
        end # eval    
            

        if DIM==1
            @eval begin

                function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 1)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_data_wf",SUF))), Ptr{Complex{$T}},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, dims[1], false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

            end # eval    
        elseif DIM==2
            @eval begin

               function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 2)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_data_wf",SUF))), Ptr{Complex{$T}},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

            end # eval    
        elseif DIM==3
            @eval begin

               function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 3)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_data_wf",SUF))), Ptr{Complex{$T}},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end
            
            end # eval    
        end

    else    
        @eval begin
            function propagate_A!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_a_wf",SUF))), Void,
                        (Ptr{Void}, $T,), psi.p, dt)
            end
    
            function propagate_B!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_b_wf",SUF))), Void,
                        (Ptr{Void}, $T,), psi.p, dt)
            end
    
            function propagate_C!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_c_wf",SUF))), Void,
                        (Ptr{Void}, $T,), psi.p, dt)
            end

            function add_apply_A!(this::($WF){$T}, other::($WF){$T},
                                 coefficient::Number=1.0)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_add_apply_a_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, ($T)), 
                         this.p, other.p, coefficient)
            end

            #function scale!(psi::($WF){$T}, factor::Number)
            #   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_scale_wf",SUF))), Void,
            #         (Ptr{Void}, ($T) ), psi.p, factor )
            #end

            function scale!(psi::($WF){$T}, factor::Number)
               ccall( ( ($(string(:c_scale_wf,SUF))), $LIBTSSM ) , Void,
                     (Ptr{Void}, ($T) ), psi.p, factor )
            end
    
            function axpy!(this::($WF){$T}, other::($WF){$T},
                                 factor::Number)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_axpy_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, ($T)), 
                         this.p, other.p, factor)
            end

            function set!(psi::($WF){$T}, f::Function)
               f_c = cfunction(f, T, (T,))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}), psi.p, f_c )
            end
    
            function set!(psi::($WF){$T}, f::Function, t::Real)
               f_c = cfunction(f, T, (T, T,))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
            end

        end # eval    

        if DIM==1
            @eval begin

                function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 1)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_data_wf",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, dims[1], false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

            end # eval    
        elseif DIM==2
            @eval begin

               function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 2)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_data_wf",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

            end # eval    
        elseif DIM==3
            @eval begin

               function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 3)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_data_wf",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end
            
            end # eval    
        end
    end        

    end #for T






    #for (T, TSSM_HANDLE) in (
    #      (:Float64, :tssm_handle),
    #     #(:Float128, :tssmq_handle), 
    #      )
    #    @eval begin
    #        function wave_function(m::($METHOD){$T})
    #            ($WF){$T}(m) 
    #        end
    #     end
    # end
 
end 

end # module

