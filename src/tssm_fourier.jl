
type Fourier1D <: TSSM
    m::Ptr{Void}
    function Fourier1D(nx::Integer, xmin::Real, xmax::Real; 
                        boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( dlsym(tssm_handle, "c_new_fourier_1d"), Ptr{Void}, 
                   (Int32, Float64, Float64, Int32), 
                   nx, xmin, xmax, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type Fourier2D <: TSSM
    m::Ptr{Void}
    function Fourier2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                        boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( dlsym(tssm_handle, "c_new_fourier_2d"), Ptr{Void}, 
                   (Int32, Float64, Float64,Int32, Float64, Float64, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type Fourier3D <: TSSM
    m::Ptr{Void}
    function Fourier3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real; boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( dlsym(tssm_handle, "c_new_fourier_3d"), Ptr{Void}, 
                   (Int32, Float64, Float64,Int32, Float64, Float64,Int32, Float64, Float64, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type FourierReal1D <: TSSM
    m::Ptr{Void}
    function FourierReal1D(nx::Integer, xmin::Real, xmax::Real; 
                        boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( dlsym(tssm_handle, "c_new_fourier_real_1d"), Ptr{Void}, 
                   (Int32, Float64, Float64, Int32), 
                   nx, xmin, xmax, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_real_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type FourierReal2D <: TSSM
    m::Ptr{Void}
    function FourierReal2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                        boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( dlsym(tssm_handle, "c_new_fourier_real_2d"), Ptr{Void}, 
                   (Int32, Float64, Float64,Int32, Float64, Float64, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_real_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type FourierReal3D <: TSSM
    m::Ptr{Void}
    function FourierReal3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real; boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        m = new( ccall( dlsym(tssm_handle, "c_new_fourier_real_3d"), Ptr{Void}, 
                   (Int32, Float64, Float64,Int32, Float64, Float64,Int32, Float64, Float64, Int32), 
                   nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_fourier_real_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

## Wave function types ############################################################################

type WfFourier1D <: WaveFunction1D
    p::Ptr{Void}
    m::Fourier1D
    function WfFourier1D( m::Fourier1D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_1d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_1d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfFourier2D <: WaveFunction2D
    p::Ptr{Void}
    m::Fourier2D
    function WfFourier2D( m::Fourier2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfFourier3D <: WaveFunction3D
    p::Ptr{Void}
    m::Fourier3D
    function WfFourier3D( m::Fourier3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_3d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_3d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfFourierReal1D <: WaveFunction1D
    p::Ptr{Void}
    m::FourierReal1D
    function WfFourierReal1D( m::FourierReal1D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_real_1d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_real_1d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfFourierReal2D <: WaveFunction2D
    p::Ptr{Void}
    m::FourierReal2D
    function WfFourierReal2D( m::FourierReal2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_real_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_real_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfFourierReal3D <: WaveFunction3D
    p::Ptr{Void}
    m::FourierReal3D
    function WfFourierReal3D( m::FourierReal3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_fourier_real_3d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wf_fourier_real_3d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

## method: wave_function ##########################################################################

function wave_function(m::Fourier1D )
    WfFourier1D(m) 
end

function wave_function(m::Fourier2D )
    WfFourier2D(m) 
end

function wave_function(m::Fourier3D )
    WfFourier3D(m) 
end

function wave_function(m::FourierReal1D )
    WfFourierReal1D(m) 
end

function wave_function(m::FourierReal2D )
    WfFourierReal2D(m) 
end

function wave_function(m::FourierReal3D )
    WfFourierReal3D(m) 
end

## method: is_real_space###########################################################################

function is_real_space(psi::WfFourier1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_1d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfFourier2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfFourier3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_3d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfFourierReal1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_real_1d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfFourierReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_real_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfFourierReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_real_3d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

## method: is_frequency_space###########################################################################

function is_frequency_space(psi::WfFourier1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_1d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfFourier2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfFourier3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_3d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfFourierReal1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_real_1d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfFourierReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_real_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfFourierReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wf_fourier_real_3d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end


## method: to_real_space! ##########################################################################

function to_real_space!(psi::WfFourier1D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfFourier2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfFourier3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfFourierReal1D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_real_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfFourierReal2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfFourierReal3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wf_fourier_real_3d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: to_frequency_space! ######################################################################

function to_frequency_space!(psi::WfFourier1D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfFourier2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfFourier3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfFourierReal1D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_real_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfFourierReal2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfFourierReal3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wf_fourier_real_3d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: propagate_A! ############################################################################

function propagate_A!(psi::WfFourier1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfFourier2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfFourier3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfFourierReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_real_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_A!(psi::WfFourierReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_A!(psi::WfFourierReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wf_fourier_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_B! ############################################################################

function propagate_B!(psi::WfFourier1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfFourier2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfFourier3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfFourierReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_B!(psi::WfFourierReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_B!(psi::WfFourierReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wf_fourier_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_C! ############################################################################

function propagate_C!(psi::WfFourier1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_C!(psi::WfFourier2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_C!(psi::WfFourier3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_C!(psi::WfFourierReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_C!(psi::WfFourierReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_C!(psi::WfFourierReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_c_wf_fourier_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end
## method: add_apply_A! ############################################################################

function add_apply_A!(this::WfFourier1D, other::WfFourier1D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfFourier2D, other::WfFourier2D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfFourier3D, other::WfFourier3D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfFourierReal1D, other::WfFourierReal1D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfFourierReal2D, other::WfFourierReal2D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfFourierReal3D, other::WfFourierReal3D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wf_fourier_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

## method: save ###################################################################################

function save(psi::WfFourier1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfFourier2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfFourier3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfFourierReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfFourierReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfFourierReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wf_fourier_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: load! ###################################################################################

function load!(psi::WfFourier1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfFourier2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfFourier3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfFourierReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfFourierReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfFourierReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wf_fourier_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: norm ###################################################################################

function norm(psi::WfFourier1D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfFourier2D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfFourier3D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfFourierReal1D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfFourierReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfFourierReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_wf_fourier_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: norm_in_frequency_space ###############################################################

function norm_in_frequency_space(psi::WfFourier1D)
   ccall( dlsym(tssm_handle, "c_norm_in_frequency_space2_wf_fourier_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfFourier2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_fourier_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfFourier3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_fourier_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfFourierReal1D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_fourier_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfFourierReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_fourier_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfFourierReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wf_fourier_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end


## method: distance ##############################################################################

function distance(psi1::WfFourier1D, psi2::WfFourier1D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_1d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfFourier2D, psi2::WfFourier2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfFourier3D, psi2::WfFourier3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_3d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfFourierReal1D, psi2::WfFourierReal1D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_real_1d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfFourierReal2D, psi2::WfFourierReal2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfFourierReal3D, psi2::WfFourierReal3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wf_fourier_real_3d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

## method: normalize! ##############################################################################

function normalize!(psi::WfFourier1D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfFourier2D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfFourier3D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfFourierReal1D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfFourierReal2D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfFourierReal3D)
   ccall( dlsym(tssm_handle, "c_normalize_wf_fourier_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: scale! ##################################################################################

function scale!(psi::WfFourier1D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_1d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfFourier2D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_2d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfFourier3D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_3d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfFourierReal1D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_real_1d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

function scale!(psi::WfFourierReal2D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_real_2d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

function scale!(psi::WfFourierReal3D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wf_fourier_real_3d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

## method: axpy! ###################################################################################

function axpy!(this::WfFourier1D, other::WfFourier1D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfFourier2D, other::WfFourier2D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfFourier3D, other::WfFourier3D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfFourierReal1D, other::WfFourierReal1D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

function axpy!(this::WfFourierReal2D, other::WfFourierReal2D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

function axpy!(this::WfFourierReal3D, other::WfFourierReal3D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wf_fourier_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

## method: get_data ###########################################################################

function get_data(psi::WfFourier1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_1d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, dims[1], false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfFourier2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_2d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfFourier3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_3d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfFourierReal1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, dims[1], false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfFourierReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfFourierReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   up = ccall( dlsym(tssm_handle, "c_get_data_wf_fourier_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end
#TODO: acces_data_frequency_space for *real* spectral methods

## method: get_eigenvalues #######################################################################

function get_eigenvalues(m::Fourier1D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev = pointer_to_array(evp, dim[1], false)   
   if unsafe_access
       return ev
   else    
       return copy(ev)
   end
end

function get_eigenvalues(m::Fourier2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::Fourier3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp3, dim[1], false)     
   if unsafe_access
       return ev1, ev2, ev3
   else
       return copy(ev1), copy(ev2), copy(ev3)
   end
end

function get_eigenvalues(m::FourierReal1D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev = pointer_to_array(evp, dim[1], false)     
   if unsafe_access
       return ev
   else
       return copy(ev)
   end
end

function get_eigenvalues(m::FourierReal2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::FourierReal3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_fourier_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp3, dim[1], false)     
   if unsafe_access
       return ev1, ev2, ev3
   else
       return copy(ev1), copy(ev2), copy(ev3)
   end
end

## method: get_nodes #############################################################################

function get_nodes(m::Fourier1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_nodes(m::Fourier2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_nodes(m::Fourier3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

function get_nodes(m::FourierReal1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_nodes(m::FourierReal2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_nodes(m::FourierReal3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_nodes_fourier_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

## method: set! ###################################################################################

function set!(psi::WfFourier1D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble,))
       ccall( dlsym(tssm_handle, "c_rset_wf_fourier_1d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble,))
       ccall( dlsym(tssm_handle, "c_set_wf_fourier_1d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end      
end

function set!(psi::WfFourier2D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wf_fourier_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wf_fourier_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end      
end

function set!(psi::WfFourier3D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wf_fourier_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wf_fourier_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end      
end

function set!(psi::WfFourierReal1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_wf_fourier_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

function set!(psi::WfFourierReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wf_fourier_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

function set!(psi::WfFourierReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wf_fourier_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

## method: set! with t argument ###########################################################################

function set!(psi::WfFourier1D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble,))
       ccall( dlsym(tssm_handle, "c_rset_t_wf_fourier_1d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble,))
       ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_1d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfFourier2D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wf_fourier_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfFourier3D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wf_fourier_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfFourierReal1D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

function set!(psi::WfFourierReal2D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

function set!(psi::WfFourierReal3D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wf_fourier_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

## method: copy! ##################################################################################

function copy!(target::WfFourier1D, source::WfFourier1D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_1d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfFourier2D, source::WfFourier2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfFourier3D, source::WfFourier3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfFourierReal1D, source::WfFourierReal1D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfFourierReal2D, source::WfFourierReal2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfFourierReal3D, source::WfFourierReal3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wf_fourier_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

## method: get_nx #############################################################################

function get_nx(m::Fourier1D)
   ccall( dlsym(tssm_handle, "c_get_nx_fourier_1d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::Fourier2D)
   ccall( dlsym(tssm_handle, "c_get_nx_fourier_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_nx_fourier_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::FourierReal1D)
   ccall( dlsym(tssm_handle, "c_get_nx_fourier_real_1d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::FourierReal2D)
   ccall( dlsym(tssm_handle, "c_get_nx_fourier_real_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_nx_fourier_real_3d"), Int32,
         (Ptr{Void}, ), m.m )
end


## method: get_xmin #############################################################################

function get_xmin(m::Fourier1D)
   ccall( dlsym(tssm_handle, "c_get_xmin_fourier_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::Fourier2D)
   ccall( dlsym(tssm_handle, "c_get_xmin_fourier_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_xmin_fourier_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::FourierReal1D)
   ccall( dlsym(tssm_handle, "c_get_xmin_fourier_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::FourierReal2D)
   ccall( dlsym(tssm_handle, "c_get_xmin_fourier_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_xmin_fourier_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_xmax #############################################################################

function get_xmax(m::Fourier1D)
   ccall( dlsym(tssm_handle, "c_get_xmax_fourier_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::Fourier2D)
   ccall( dlsym(tssm_handle, "c_get_xmax_fourier_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_xmax_fourier_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::FourierReal1D)
   ccall( dlsym(tssm_handle, "c_get_xmax_fourier_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::FourierReal2D)
   ccall( dlsym(tssm_handle, "c_get_xmax_fourier_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_xmax_fourier_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_ny #############################################################################

function get_ny(m::Fourier2D)
   ccall( dlsym(tssm_handle, "c_get_ny_fourier_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_ny(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_ny_fourier_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_ny(m::FourierReal2D)
   ccall( dlsym(tssm_handle, "c_get_ny_fourier_real_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_ny(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_ny_fourier_real_3d"), Int32,
         (Ptr{Void}, ), m.m )
end


## method: get_ymin #############################################################################

function get_ymin(m::Fourier2D)
   ccall( dlsym(tssm_handle, "c_get_ymin_fourier_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymin(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_ymin_fourier_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymin(m::FourierReal2D)
   ccall( dlsym(tssm_handle, "c_get_ymin_fourier_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymin(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_ymin_fourier_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_ymax #############################################################################

function get_ymax(m::Fourier2D)
   ccall( dlsym(tssm_handle, "c_get_ymax_fourier_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymax(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_ymax_fourier_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymax(m::FourierReal2D)
   ccall( dlsym(tssm_handle, "c_get_ymax_fourier_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymax(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_ymax_fourier_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_nz #############################################################################

function get_nz(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_nz_fourier_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nz(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_nz_fourier_real_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

## method: get_zmin #############################################################################
function get_zmin(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_zmin_fourier_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_zmin(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_zmin_fourier_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_zmax #############################################################################

function get_zmax(m::Fourier3D)
   ccall( dlsym(tssm_handle, "c_get_zmax_fourier_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_zmax(m::FourierReal3D)
   ccall( dlsym(tssm_handle, "c_get_zmax_fourier_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end


# end
