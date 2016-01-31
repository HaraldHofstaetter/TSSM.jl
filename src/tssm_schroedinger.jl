
type Schroedinger1D <: TSSM
    m::Ptr{Void}
    function Schroedinger1D(nx::Integer, xmin::Real, xmax::Real; 
                        hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                        cubic_coupling::Real=0.0,
                        boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_1D
        V_c = cfunction(potential, Cdouble, (Cdouble,))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_1d"), Ptr{Void}, 
                   (Int32, Float64, Float64, Float64, Float64, Ptr{Void}, Bool, Float64, Int32), 
                   nx, xmin, xmax, hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m          
    end
end

type Schroedinger2D <: TSSM
    m::Ptr{Void}
    function Schroedinger2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                        hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                        cubic_coupling::Real=0.0,
                        boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_2D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_2d"), Ptr{Void}, 
               (Int32, Float64, Float64, Int32, Float64, Float64, Float64, Float64, Ptr{Void}, Bool, Float64, Int32), 
               nx, xmin, xmax, ny, ymin, ymax, hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type Schroedinger3D <: TSSM
    m::Ptr{Void}
    function Schroedinger3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_3D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_3d"), Ptr{Void}, 
               (Int32, Float64, Float64, Int32, Float64, Float64, Int32, Float64, Float64, Float64, Float64, Ptr{Void}, Bool, Float64, Int32), 
               nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type SchroedingerReal1D <: TSSM
    m::Ptr{Void}
    function SchroedingerReal1D(nx::Integer, xmin::Real, xmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_1D
        V_c = cfunction(potential, Cdouble, (Cdouble,))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_real_1d"), Ptr{Void}, 
               (Int32, Float64, Float64, Float64, Float64, Ptr{Void}, Bool, Float64, Int32), 
               nx, xmin, xmax, hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_real_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type SchroedingerReal2D <: TSSM
    m::Ptr{Void}
    function SchroedingerReal2D(nx::Integer, xmin::Real, xmax::Real, ny::Integer, ymin::Real, ymax::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_2D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_real_2d"), Ptr{Void}, 
               (Int32, Float64, Float64, Int32, Float64, Float64, Float64, Float64, Ptr{Void}, Bool, Float64, Int32), 
               nx, xmin, xmax, ny, ymin, ymax, hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_real_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

type SchroedingerReal3D <: TSSM
    m::Ptr{Void}
    function SchroedingerReal3D(nx::Integer, xmin::Real, xmax::Real,  ny::Integer, ymin::Real, ymax::Real,
                       nz::Integer, zmin::Real, zmax::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       cubic_coupling::Real=0.0,
                       boundary_conditions::Integer=periodic)
        ccall( dlsym(tssm.tssm_handle, "c_initialize_tssm_fourier"), Void, ())
        with_potential = potential!=none_3D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_real_3d"), Ptr{Void}, 
              (Int32, Float64, Float64, Int32, Float64, Float64, Int32, Float64, Float64, Float64, Float64, Ptr{Void}, Bool, Float64, Int32), 
               nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, hbar, mass, V_c, with_potential, cubic_coupling, boundary_conditions) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_real_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m           
    end
end

## Wave function types ############################################################################

type WfSchroedinger1D <: WaveFunction1D
    p::Ptr{Void}
    m::Schroedinger1D
    function WfSchroedinger1D( m::Schroedinger1D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_1d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_1d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedinger2D <: WaveFunction2D
    p::Ptr{Void}
    m::Schroedinger2D
    function WfSchroedinger2D( m::Schroedinger2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedinger3D <: WaveFunction3D
    p::Ptr{Void}
    m::Schroedinger3D
    function WfSchroedinger3D( m::Schroedinger3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_3d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_3d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerReal1D <: WaveFunction1D
    p::Ptr{Void}
    m::SchroedingerReal1D
    function WfSchroedingerReal1D( m::SchroedingerReal1D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_real_1d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_real_1d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerReal2D <: WaveFunction2D
    p::Ptr{Void}
    m::SchroedingerReal2D
    function WfSchroedingerReal2D( m::SchroedingerReal2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_real_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_real_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerReal3D <: WaveFunction3D
    p::Ptr{Void}
    m::SchroedingerReal3D
    function WfSchroedingerReal3D( m::SchroedingerReal3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_real_3d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_real_3d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    


## method: wave_function ##########################################################################

function wave_function(m::Schroedinger1D )
    WfSchroedinger1D(m) 
end

function wave_function(m::Schroedinger2D )
    WfSchroedinger2D(m) 
end

function wave_function(m::Schroedinger3D )
    WfSchroedinger3D(m) 
end

function wave_function(m::SchroedingerReal1D )
    WfSchroedingerReal1D(m) 
end

function wave_function(m::SchroedingerReal2D )
    WfSchroedingerReal2D(m) 
end

function wave_function(m::SchroedingerReal3D )
    WfSchroedingerReal3D(m) 
end

## method: is_real_space ##########################################################################

function is_real_space(psi::WfSchroedinger1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_1d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedinger2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedinger3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_3d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerReal1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_real_1d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_real_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_real_3d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

## method: is_frequency_space###########################################################################

function is_frequency_space(psi::WfSchroedinger1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_1d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedinger2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedinger3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_3d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerReal1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_real_1d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_real_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_real_3d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end


## method: to_real_space! #########################################################################

function to_real_space!(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_real_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_real_3d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: to_frequency_space! ####################################################################

function to_frequency_space!(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_real_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_real_3d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: propagate_A! ###########################################################################

function propagate_A!(psi::WfSchroedinger1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfSchroedinger2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfSchroedinger3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_B! ###########################################################################

function propagate_B!(psi::WfSchroedinger1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfSchroedinger2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfSchroedinger3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: add_apply_A! ###########################################################################

function add_apply_A!(this::WfSchroedinger1D, other::WfSchroedinger1D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedinger2D, other::WfSchroedinger2D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedinger3D, other::WfSchroedinger3D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerReal1D, other::WfSchroedingerReal1D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerReal2D, other::WfSchroedingerReal2D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerReal3D, other::WfSchroedingerReal3D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

## method: save ###################################################################################

function save(psi::WfSchroedinger1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedinger2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedinger3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: load! ##################################################################################

function load!(psi::WfSchroedinger1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedinger2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedinger3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: norm ###################################################################################

function norm(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: norm_in_frequency_space ###############################################################

function norm_in_frequency_space(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end


## method: distance ##############################################################################

function distance(psi1::WfSchroedinger1D, psi2::WfSchroedinger1D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_1d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedinger2D, psi2::WfSchroedinger2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedinger3D, psi2::WfSchroedinger3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_3d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerReal1D, psi2::WfSchroedingerReal1D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_real_1d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerReal2D, psi2::WfSchroedingerReal2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerReal3D, psi2::WfSchroedingerReal3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_real_3d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

## method: normalize! #############################################################################

function normalize!(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: scale! #################################################################################

function scale!(psi::WfSchroedinger1D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wfs_1d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfSchroedinger2D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wfs_2d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfSchroedinger3D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wfs_3d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfSchroedingerReal1D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wfs_real_1d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

function scale!(psi::WfSchroedingerReal2D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wfs_real_2d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

function scale!(psi::WfSchroedingerReal3D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wfs_real_3d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

## method: axpy! ##################################################################################

function axpy!(this::WfSchroedinger1D, other::WfSchroedinger1D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedinger2D, other::WfSchroedinger2D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedinger3D, other::WfSchroedinger3D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerReal1D, other::WfSchroedingerReal1D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerReal2D, other::WfSchroedingerReal2D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerReal3D, other::WfSchroedingerReal3D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

## method: get_data ###########################################################################

function get_data(psi::WfSchroedinger1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_1d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, dims[1], false)   
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedinger2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_2d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedinger3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_3d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerReal1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, dims[1], false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_real_3d"), Ptr{Float64},
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

function get_eigenvalues(m::Schroedinger1D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev = pointer_to_array(evp, dim[1], false)     
   if unsafe_access
       return ev
   else
       return copy(ev)
   end
end

function get_eigenvalues(m::Schroedinger2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::Schroedinger3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp3, dim[1], false)     
   if unsafe_access
       return ev1, ev2, ev3
   else
       return copy(ev1), copy(ev2), copy(ev3)
   end
end

function get_eigenvalues(m::SchroedingerReal1D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev = pointer_to_array(evp, dim[1], false)     
   if unsafe_access
       return ev
   else
       return copy(ev)
   end
end

function get_eigenvalues(m::SchroedingerReal2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       ev1, ev2
   else
       copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::SchroedingerReal3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp3, dim[1], false)     
   if unsafe_access
       ev1, ev2, ev3
   else
       copy(ev1), copy(ev2), copy(ev3)
   end
end

## method: get_nodes #############################################################################

function get_nodes(m::Schroedinger1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_nodes(m::Schroedinger2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_nodes(m::Schroedinger3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

function get_nodes(m::SchroedingerReal1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_nodes(m::SchroedingerReal2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_nodes(m::SchroedingerReal3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

## method: set! ##################################################################################

function set!(psi::WfSchroedinger1D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble,))
       ccall( dlsym(tssm_handle, "c_rset_wfs_1d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble,))
       ccall( dlsym(tssm_handle, "c_set_wfs_1d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end
end

function set!(psi::WfSchroedinger2D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wfs_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wfs_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end
end

function set!(psi::WfSchroedinger3D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wfs_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wfs_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end
end

function set!(psi::WfSchroedingerReal1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_wfs_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

function set!(psi::WfSchroedingerReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wfs_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

function set!(psi::WfSchroedingerReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wfs_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

## method: set! with t argument ###########################################################################

function set!(psi::WfSchroedinger1D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble,))
       ccall( dlsym(tssm_handle, "c_rset_t_wfs_1d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble,))
       ccall( dlsym(tssm_handle, "c_set_t_wfs_1d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfSchroedinger2D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wfs_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wfs_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfSchroedinger3D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wfs_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wfs_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfSchroedingerReal1D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_t_wfs_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

function set!(psi::WfSchroedingerReal2D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wfs_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

function set!(psi::WfSchroedingerReal3D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wfs_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end


## method: copy! #################################################################################

function copy!(target::WfSchroedinger1D, source::WfSchroedinger1D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_1d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedinger2D, source::WfSchroedinger2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedinger3D, source::WfSchroedinger3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerReal1D, source::WfSchroedingerReal1D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerReal2D, source::WfSchroedingerReal2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerReal3D, source::WfSchroedingerReal3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

## method: get_nx #############################################################################

function get_nx(m::Schroedinger1D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_1d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_real_1d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_real_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_real_3d"), Int32,
         (Ptr{Void}, ), m.m )
end


## method: get_xmin #############################################################################

function get_xmin(m::Schroedinger1D)
   ccall( dlsym(tssm_handle, "c_get_xmin_schroedinger_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_xmin_schroedinger_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_xmin_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::SchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_get_xmin_schroedinger_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_xmin_schroedinger_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmin(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_xmin_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_xmax #############################################################################

function get_xmax(m::Schroedinger1D)
   ccall( dlsym(tssm_handle, "c_get_xmax_schroedinger_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_xmax_schroedinger_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_xmax_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::SchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_get_xmax_schroedinger_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_xmax_schroedinger_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_xmax(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_xmax_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_ny #############################################################################

function get_ny(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_ny_schroedinger_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_ny(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_ny_schroedinger_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_ny(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_ny_schroedinger_real_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_ny(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_ny_schroedinger_real_3d"), Int32,
         (Ptr{Void}, ), m.m )
end


## method: get_ymin #############################################################################

function get_ymin(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_ymin_schroedinger_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymin(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_ymin_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymin(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_ymin_schroedinger_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymin(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_ymin_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_ymax #############################################################################

function get_ymax(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_ymax_schroedinger_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymax(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_ymax_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymax(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_ymax_schroedinger_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_ymax(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_ymax_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_nz #############################################################################

function get_nz(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_nz_schroedinger_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nz(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_nz_schroedinger_real_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

## method: get_zmin #############################################################################
function get_zmin(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_zmin_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_zmin(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_zmin_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_zmax #############################################################################

function get_zmax(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_zmax_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_zmax(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_zmax_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! schroedinger specific methods !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## method: get_hbar #############################################################################

function get_hbar(m::Schroedinger1D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_mass #############################################################################

function get_mass(m::Schroedinger1D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_cubic_coupling #############################################################################

function get_cubic_coupling(m::Schroedinger1D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::Schroedinger2D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::Schroedinger3D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end



## method: set_potential! ##########################################################################

function set_potential!(m::Schroedinger1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_1d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::Schroedinger2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_2d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::Schroedinger3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_3d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerReal1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

## method: get_potential ###########################################################################

function get_potential(m::Schroedinger1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::Schroedinger2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::Schroedinger3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], dims[3], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerReal1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], dims[3], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

## method: load_potential! #########################################################################

function load_potential!(m::Schroedinger1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::Schroedinger2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::Schroedinger3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

## method: save_potential #########################################################################

function save_potential(m::Schroedinger1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::Schroedinger2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::Schroedinger3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

## method: imaginary_time_propagate_A! ############################################################

function imaginary_time_propagate_A!(psi::WfSchroedinger1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedinger2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedinger3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: imaginary_time_propagate_B! ############################################################

function imaginary_time_propagate_B!(psi::WfSchroedinger1D, dt::Number, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_1d"), Void,
                    (Ptr{Void}, Complex{Float64}, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedinger2D, dt::Number, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_2d"), Void,
                    (Ptr{Void}, Complex{Float64}, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedinger3D, dt::Number, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_3d"), Void,
                    (Ptr{Void}, Complex{Float64}, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerReal1D, dt::Real, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_real_1d"), Void,
                    (Ptr{Void}, Float64, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerReal2D, dt::Real, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_real_2d"), Void,
                    (Ptr{Void}, Float64, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerReal3D, dt::Real, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_real_3d"), Void,
                    (Ptr{Void}, Float64, Int32), psi.p, dt, method_for_B)
end

## method: add_apply_B! ###########################################################################

function add_apply_B!(this::WfSchroedinger1D, other::WfSchroedinger1D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedinger2D, other::WfSchroedinger2D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedinger3D, other::WfSchroedinger3D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerReal1D, other::WfSchroedingerReal1D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerReal2D, other::WfSchroedingerReal2D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerReal3D, other::WfSchroedingerReal3D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

## method: kinetic_energy #########################################################################

function kinetic_energy(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: potential_energy #######################################################################

function potential_energy(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: interaction_energy #####################################################################

function interaction_energy(psi::WfSchroedinger1D)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedinger2D)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedinger3D)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerReal1D)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerReal2D)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerReal3D)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: observable #############################################################################

function observable(psi::WfSchroedinger1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, ))
   ccall( dlsym(tssm_handle, "c_observable_wfs_1d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedinger2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_2d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c)
end

function observable(psi::WfSchroedinger3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_3d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedingerReal1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, ))
   ccall( dlsym(tssm_handle, "c_observable_wfs_real_1d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedingerReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedingerReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_real_3d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

## method: get_energy_expectation_deviation #######################################################

function get_energy_expectation_deviation(psi::WfSchroedinger1D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedinger2D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedinger3D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerReal1D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_real_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerReal2D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_real_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerReal3D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_real_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

## method: get_realspace_observables ##############################################################

function get_realspace_observables(psi::WfSchroedinger1D)
   ans =Array(Float64, 4)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedinger2D)
   ans =Array(Float64, 6)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedinger3D)
   ans =Array(Float64, 8)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerReal1D)
   ans =Array(Float64, 4)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_real_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerReal2D)
   ans =Array(Float64, 6)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_real_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerReal3D)
   ans =Array(Float64, 8)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_real_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end


## method: selfconsistent_nonlinear_step! ##########################################################

function selfconsistent_nonlinear_step!(psi::WfSchroedinger1D, dt::Number, 
                  dt1::Number, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_1d"), Void,
                    (Ptr{Void}, Complex{Float64}, Complex{Float64}, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedinger2D, dt::Number,
                  dt1::Number, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_2d"), Void,
                    (Ptr{Void}, Complex{Float64}, Complex{Float64}, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedinger3D, dt::Number,
                  dt1::Number, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_3d"), Void,
                    (Ptr{Void}, Complex{Float64}, Complex{Float64}, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerReal1D, dt::Real,
                  dt1::Real, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_real_1d"), Void,
                    (Ptr{Void}, Float64, Float64, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerReal2D, dt::Real,
                  dt1::Real, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_real_2d"), Void,
                    (Ptr{Void}, Float64, Float64, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerReal3D, dt::Real,
                  dt1::Real, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_real_3d"), Void,
                    (Ptr{Void}, Float64, Float64, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end



# end
