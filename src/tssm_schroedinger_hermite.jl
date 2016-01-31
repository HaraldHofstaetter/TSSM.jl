type SchroedingerHermite1D <: TSSM
    m::Ptr{Void}
    function SchroedingerHermite1D(nx::Integer, omega_x::Real; 
                        hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                        cubic_coupling::Real=0.0)
        with_potential = potential!=none_1D
        V_c = cfunction(potential, Cdouble, (Cdouble,))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_hermite_1d"), Ptr{Void}, 
               (Int32, Float64, Float64, Float64, Ptr{Void}, Bool, Float64), 
               nx, omega_x, hbar, mass, V_c, with_potential, cubic_coupling) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_hermite_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m     
    end
end

type SchroedingerHermite2D <: TSSM
    m::Ptr{Void}
    function SchroedingerHermite2D(nx::Integer, omega_x::Real, ny::Integer, omega_y::Real;
                        hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                        cubic_coupling::Real=0.0)
        with_potential = potential!=none_2D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_hermite_2d"), Ptr{Void}, 
               (Int32, Float64, Int32, Float64, Float64, Float64, Ptr{Void}, Bool, Float64), 
               nx, omega_x, ny, omega_y, hbar, mass, V_c, with_potential, cubic_coupling) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_hermite_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m     
    end
end

type SchroedingerHermite3D <: TSSM
    m::Ptr{Void}
    function SchroedingerHermite3D(nx::Integer, omega_x::Real,  ny::Integer, omega_y::Real,
                       nz::Integer, omega_z::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_3D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_hermite_3d"), Ptr{Void}, 
               (Int32, Float64, Int32, Float64, Int32, Float64, Float64, Float64, Ptr{Void}, Bool, Float64), 
               nx, omega_x, ny, omega_y, nz, omega_z, hbar, mass, V_c, with_potential, cubic_coupling) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_hermite_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m     
    end
end

type SchroedingerHermiteReal1D <: TSSM
    m::Ptr{Void}
    function SchroedingerHermiteReal1D(nx::Integer, omega_x::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_1D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_1D
        V_c = cfunction(potential, Cdouble, (Cdouble,))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_hermite_real_1d"), Ptr{Void}, 
               (Int32, Float64, Float64, Float64, Ptr{Void}, Bool, Float64), 
               nx, omega_x, hbar, mass, V_c, with_potential, cubic_coupling) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_hermite_real_1d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m     
    end
end

type SchroedingerHermiteReal2D <: TSSM
    m::Ptr{Void}
    function SchroedingerHermiteReal2D(nx::Integer, omega_x::Real, ny::Integer, omega_y::Real;
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_2D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_2D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_hermite_real_2d"), Ptr{Void}, 
               (Int32, Float64, Int32, Float64, Float64, Float64, Ptr{Void}, Bool, Float64), 
               nx, omega_x, ny, omega_y, hbar, mass, V_c, with_potential, cubic_coupling) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_hermite_real_2d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m     
    end
end

type SchroedingerHermiteReal3D <: TSSM
    m::Ptr{Void}
    function SchroedingerHermiteReal3D(nx::Integer, omega_x::Real,  ny::Integer, omega_y::Real,
                       nz::Integer, omega_z::Real; 
                       hbar::Real=1.0, mass::Real=1.0, potential::Function=none_3D,
                       cubic_coupling::Real=0.0)
        with_potential = potential!=none_3D
        V_c = cfunction(potential, Cdouble, (Cdouble,Cdouble,Cdouble))
        m = new( ccall( dlsym(tssm_handle, "c_new_schroedinger_hermite_real_3d"), Ptr{Void}, 
               (Int32, Float64, Int32, Float64, Int32, Float64, Float64, Float64, Ptr{Void}, Bool, Float64), 
               nx, omega_x, ny, omega_y, nz, omega_z, hbar, mass, V_c, with_potential, cubic_coupling) )
        finalizer(m, x -> ccall( dlsym(tssm_handle, "c_finalize_schroedinger_hermite_real_3d"), Void, (Ptr{Ptr{Void}},), &x.m) )
        m     
    end
end

## Wave function types ############################################################################

type WfSchroedingerHermite1D <: WaveFunction1D
    p::Ptr{Void}
    m::SchroedingerHermite1D
    function WfSchroedingerHermite1D( m::SchroedingerHermite1D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_hermite_1d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_hermite_1d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerHermite2D <: WaveFunction2D
    p::Ptr{Void}
    m::SchroedingerHermite2D
    function WfSchroedingerHermite2D( m::SchroedingerHermite2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_hermite_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_hermite_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerHermite3D <: WaveFunction3D
    p::Ptr{Void}
    m::SchroedingerHermite3D
    function WfSchroedingerHermite3D( m::SchroedingerHermite3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_hermite_3d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_hermite_3d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerHermiteReal1D <: WaveFunction1D
    p::Ptr{Void}
    m::SchroedingerHermiteReal1D
    function WfSchroedingerHermiteReal1D( m::SchroedingerHermiteReal1D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_hermite_real_1d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_hermite_real_1d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerHermiteReal2D <: WaveFunction2D
    p::Ptr{Void}
    m::SchroedingerHermiteReal2D
    function WfSchroedingerHermiteReal2D( m::SchroedingerHermiteReal2D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_hermite_real_2d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_hermite_real_2d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

type WfSchroedingerHermiteReal3D <: WaveFunction3D
    p::Ptr{Void}
    m::SchroedingerHermiteReal3D
    function WfSchroedingerHermiteReal3D( m::SchroedingerHermiteReal3D )
        wf = new( ccall( dlsym(tssm_handle, "c_new_wf_schroedinger_hermite_real_3d"), Ptr{Void},
                    (Ptr{Void},), m.m) , m)   
        finalizer(wf, x -> ccall( dlsym(tssm_handle, "c_finalize_wfs_hermite_real_3d"), Void, (Ptr{Ptr{Void}},), &x.p) )
        wf
    end
end    

## method: wave_function ##########################################################################

function wave_function(m::SchroedingerHermite1D )
    WfSchroedingerHermite1D(m) 
end

function wave_function(m::SchroedingerHermite2D )
    WfSchroedingerHermite2D(m) 
end

function wave_function(m::SchroedingerHermite3D )
    WfSchroedingerHermite3D(m) 
end

function wave_function(m::SchroedingerHermiteReal1D )
    WfSchroedingerHermiteReal1D(m) 
end

function wave_function(m::SchroedingerHermiteReal2D )
    WfSchroedingerHermiteReal2D(m) 
end

function wave_function(m::SchroedingerHermiteReal3D )
    WfSchroedingerHermiteReal3D(m) 
end

## method: is_real_space ##########################################################################

function is_frequency_space(psi::WfSchroedingerHermite1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_1d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerHermite2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerHermite3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_3d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerHermiteReal1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_real_1d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerHermiteReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_real_2d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end

function is_frequency_space(psi::WfSchroedingerHermiteReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_real_3d"), Int32,
                    (Ptr{Void},), psi.p) != 1
end


## method: is_real_space ##########################################################################

function is_real_space(psi::WfSchroedingerHermite1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_1d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerHermite2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerHermite3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_3d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerHermiteReal1D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_real_1d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerHermiteReal2D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_real_2d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

function is_real_space(psi::WfSchroedingerHermiteReal3D)
    ccall( dlsym(tssm_handle, "c_is_real_space_wfs_hermite_real_3d"), Int32,
                    (Ptr{Void},), psi.p) == 1
end

## method: to_real_space! #########################################################################

function to_real_space!(psi::WfSchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_hermite_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_hermite_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_hermite_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_hermite_real_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_hermite_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_real_space!(psi::WfSchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_to_real_space_wfs_hermite_real_3d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: to_frequency_space! ####################################################################

function to_frequency_space!(psi::WfSchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_hermite_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_hermite_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_hermite_3d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_hermite_real_1d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_hermite_real_2d"), Void,
                    (Ptr{Void},), psi.p)
end

function to_frequency_space!(psi::WfSchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_to_frequency_space_wfs_hermite_real_3d"), Void,
                    (Ptr{Void},), psi.p)
end

## method: propagate_A! ###########################################################################

function propagate_A!(psi::WfSchroedingerHermite1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerHermite2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerHermite3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerHermiteReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_hermite_real_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerHermiteReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_A!(psi::WfSchroedingerHermiteReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_a_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: propagate_B! ###########################################################################

function propagate_B!(psi::WfSchroedingerHermite1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerHermite2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerHermite3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerHermiteReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_hermite_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerHermiteReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function propagate_B!(psi::WfSchroedingerHermiteReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_propagate_b_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: add_apply_A! ###########################################################################

function add_apply_A!(this::WfSchroedingerHermite1D, other::WfSchroedingerHermite1D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerHermite2D, other::WfSchroedingerHermite2D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerHermite3D, other::WfSchroedingerHermite3D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerHermiteReal1D, other::WfSchroedingerHermiteReal1D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_hermite_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerHermiteReal2D, other::WfSchroedingerHermiteReal2D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_A!(this::WfSchroedingerHermiteReal3D, other::WfSchroedingerHermiteReal3D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_a_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

## method: save ###################################################################################

function save(psi::WfSchroedingerHermite1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_hermite_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerHermite2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_hermite_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerHermite3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerHermiteReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerHermiteReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function save(psi::WfSchroedingerHermiteReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_wfs_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: load! ##################################################################################

function load!(psi::WfSchroedingerHermite1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_hermite_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerHermite2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_hermite_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerHermite3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerHermiteReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerHermiteReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

function load!(psi::WfSchroedingerHermiteReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_wfs_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), psi.p, filename, length(filename))
end    

## method: norm ###################################################################################

function norm(psi::WfSchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_hermite_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_hermite_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm(psi::WfSchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: norm_in_frequency_space ###############################################################

function norm_in_frequency_space(psi::WfSchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_hermite_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_hermite_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function norm_in_frequency_space(psi::WfSchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_norm2_in_frequency_space_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end


## method: distance ##############################################################################

function distance(psi1::WfSchroedingerHermite1D, psi2::WfSchroedingerHermite1D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_hermite_1d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerHermite2D, psi2::WfSchroedingerHermite2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_hermite_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerHermite3D, psi2::WfSchroedingerHermite3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_hermite_3d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerHermiteReal1D, psi2::WfSchroedingerHermiteReal1D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerHermiteReal2D, psi2::WfSchroedingerHermiteReal2D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

function distance(psi1::WfSchroedingerHermiteReal3D, psi2::WfSchroedingerHermiteReal3D)
   if psi1.m ≠ psi2.m
       error("psi1 and psi2 must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_distance_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
end

## method: normalize! #############################################################################

function normalize!(psi::WfSchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_hermite_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_hermite_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function normalize!(psi::WfSchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_normalize_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: scale! #################################################################################

function scale!(psi::WfSchroedingerHermite1D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wfs_hermite_1d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfSchroedingerHermite2D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wfs_hermite_2d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfSchroedingerHermite3D, factor::Number)
   ccall( dlsym(tssm_handle, "c_scale_wfs_hermite_3d"), Void,
         (Ptr{Void}, Complex{Float64} ), psi.p, factor )
end

function scale!(psi::WfSchroedingerHermiteReal1D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wfs_hermite_real_1d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

function scale!(psi::WfSchroedingerHermiteReal2D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wfs_hermite_real_2d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

function scale!(psi::WfSchroedingerHermiteReal3D, factor::Real)
   ccall( dlsym(tssm_handle, "c_scale_wfs_hermite_real_3d"), Void,
         (Ptr{Void}, Float64 ), psi.p, factor )
end

## method: axpy! ##################################################################################

function axpy!(this::WfSchroedingerHermite1D, other::WfSchroedingerHermite1D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerHermite2D, other::WfSchroedingerHermite2D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerHermite3D, other::WfSchroedingerHermite3D,
                     factor::Number)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerHermiteReal1D, other::WfSchroedingerHermiteReal1D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_hermite_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerHermiteReal2D, other::WfSchroedingerHermiteReal2D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

function axpy!(this::WfSchroedingerHermiteReal3D, other::WfSchroedingerHermiteReal3D,
                     factor::Real)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_axpy_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, factor)
end

## method: get_data ###########################################################################

function get_data(psi::WfSchroedingerHermite1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_hermite_1d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, dims[1], false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerHermite2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_hermite_2d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerHermite3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_hermite_3d"), Ptr{Complex{Float64}},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerHermiteReal1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_hermite_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, dims[1], false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerHermiteReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
   data = pointer_to_array(up, (dims[1], dims[2]), false)     
   if unsafe_access
      return data
   else
      return copy(data)
   end
end

function get_data(psi::WfSchroedingerHermiteReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   up = ccall( dlsym(tssm_handle, "c_get_data_wfs_hermite_real_3d"), Ptr{Float64},
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

function get_eigenvalues(m::SchroedingerHermite1D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev = pointer_to_array(evp, dim[1], false)  
   if unsafe_access
       return ev
   else
       return copy(ev)
   end
end

function get_eigenvalues(m::SchroedingerHermite2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::SchroedingerHermite3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp3, dim[1], false)     
   if unsafe_access
       return ev1, ev2, ev3
   else
       return copy(ev1), copy(ev2), copy(ev3)
   end
end

function get_eigenvalues(m::SchroedingerHermiteReal1D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev = pointer_to_array(evp, dim[1], false)     
   if unsafe_access
       return ev
   else
       return copy(ev)
   end
end

function get_eigenvalues(m::SchroedingerHermiteReal2D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   if unsafe_access
       return ev1, ev2
   else
       return copy(ev1), copy(ev2)
   end
end

function get_eigenvalues(m::SchroedingerHermiteReal3D, unsafe_access::Bool=false)
   dim =Array(Int32, 1)
   evp1 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   ev1 = pointer_to_array(evp1, dim[1], false)     
   evp2 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   ev2 = pointer_to_array(evp2, dim[1], false)     
   evp3 = ccall( dlsym(tssm_handle, "c_get_eigenvalues_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   ev3 = pointer_to_array(evp3, dim[1], false)     
   if unsafe_access
       return ev1, ev2, ev3
   else
       return copy(ev1), copy(ev2), copy(ev3)
   end
end

## method: get_nodes #############################################################################

function get_nodes(m::SchroedingerHermite1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_nodes(m::SchroedingerHermite2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_nodes(m::SchroedingerHermite3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

function get_nodes(m::SchroedingerHermiteReal1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_nodes(m::SchroedingerHermiteReal2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_nodes(m::SchroedingerHermiteReal3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_nodes_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

## method: get_weights #############################################################################

function get_weights(m::SchroedingerHermite1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_weights(m::SchroedingerHermite2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_weights(m::SchroedingerHermite3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

function get_weights(m::SchroedingerHermiteReal1D)
   dim =Array(Int32, 1)
   np = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n = pointer_to_array(np, dim[1], false)     
   copy(n)
end

function get_weights(m::SchroedingerHermiteReal2D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   copy(n1), copy(n2)
end

function get_weights(m::SchroedingerHermiteReal3D)
   dim =Array(Int32, 1)
   np1 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
   n1 = pointer_to_array(np1, dim[1], false)     
   np2 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
   n2 = pointer_to_array(np2, dim[1], false)     
   np3 = ccall( dlsym(tssm_handle, "c_get_weights_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
   n3 = pointer_to_array(np3, dim[1], false)     
   copy(n1), copy(n2), copy(n3)
end

## method: get_H #######################################################################

function get_H(m::SchroedingerHermite1D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
   H = pointer_to_array(Hp, (dims[1], dims[2]), false)  
   if unsafe_access
       return H
   else
       return copy(H)
   end
end

function get_H(m::SchroedingerHermite2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp1 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
   H1 = pointer_to_array(Hp1, (dims[1], dims[2]), false)  
   Hp2 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 2 )
   H2 = pointer_to_array(Hp2, (dims[1], dims[2]), false)  
   if unsafe_access
       return H1, H2
   else
       return copy(H1), copy(H2)
   end
end

function get_H(m::SchroedingerHermite3D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp1 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
   H1 = pointer_to_array(Hp1, (dims[1], dims[2]), false)  
   Hp2 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 2 )
   H2 = pointer_to_array(Hp2, (dims[1], dims[2]), false)  
   Hp3 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 3 )
   H3 = pointer_to_array(Hp3, (dims[1], dims[2]), false)  
   if unsafe_access
       return H1, H2, H3
   else
       return copy(H1), copy(H2), copy(H3)
   end
end

function get_H(m::SchroedingerHermiteReal1D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
   H = pointer_to_array(Hp, (dims[1], dims[2]), false)  
   if unsafe_access
       return H
   else
       return copy(H)
   end
end

function get_H(m::SchroedingerHermiteReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp1 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
   H1 = pointer_to_array(Hp1, (dims[1], dims[2]), false)  
   Hp2 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 2 )
   H2 = pointer_to_array(Hp2, (dims[1], dims[2]), false)  
   if unsafe_access
       return H1, H2
   else
       return copy(H1), copy(H2)
   end
end

function get_H(m::SchroedingerHermiteReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Hp1 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 1 )
   H1 = pointer_to_array(Hp1, (dims[1], dims[2]), false)  
   Hp2 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 2 )
   H2 = pointer_to_array(Hp2, (dims[1], dims[2]), false)  
   Hp3 = ccall( dlsym(tssm_handle, "c_get_h_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dims, 3 )
   H3 = pointer_to_array(Hp3, (dims[1], dims[2]), false)  
   if unsafe_access
       return H1, H2, H3
   else
       return copy(H1), copy(H2), copy(H3)
   end
end



## method: set! ###################################################################################

function set!(psi::WfSchroedingerHermite1D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble,))
       ccall( dlsym(tssm_handle, "c_rset_wfs_hermite_1d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble,))
       ccall( dlsym(tssm_handle, "c_set_wfs_hermite_1d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end
end

function set!(psi::WfSchroedingerHermite2D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wfs_hermite_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wfs_hermite_2d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end
end

function set!(psi::WfSchroedingerHermite3D, f::Function)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_wfs_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_wfs_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
   end
end

function set!(psi::WfSchroedingerHermiteReal1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_wfs_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

function set!(psi::WfSchroedingerHermiteReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wfs_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

function set!(psi::WfSchroedingerHermiteReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_wfs_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
end

## method: set! with t argument ###########################################################################

function set!(psi::WfSchroedingerHermite1D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble,))
       ccall( dlsym(tssm_handle, "c_rset_t_wfs_hermite_1d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble,))
       ccall( dlsym(tssm_handle, "c_set_t_wfs_hermite_1d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfSchroedingerHermite2D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wfs_hermite_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wfs_hermite_2d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfSchroedingerHermite3D, f::Function, t::Real)
   try
       f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_rset_t_wfs_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   catch
       f_c = cfunction(f, Complex{Cdouble}, (Cdouble, Cdouble, Cdouble, Cdouble))
       ccall( dlsym(tssm_handle, "c_set_t_wfs_hermite_3d"), Void,
             (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
   end      
end

function set!(psi::WfSchroedingerHermiteReal1D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_t_wfs_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

function set!(psi::WfSchroedingerHermiteReal2D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wfs_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

function set!(psi::WfSchroedingerHermiteReal3D, f::Function, t::Real)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_t_wfs_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}, Float64), psi.p, f_c, t )
end

## method: copy! ##################################################################################

function copy!(target::WfSchroedingerHermite1D, source::WfSchroedingerHermite1D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_hermite_1d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerHermite2D, source::WfSchroedingerHermite2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_hermite_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerHermite3D, source::WfSchroedingerHermite3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerHermiteReal1D, source::WfSchroedingerHermiteReal1D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerHermiteReal2D, source::WfSchroedingerHermiteReal2D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

function copy!(target::WfSchroedingerHermiteReal3D, source::WfSchroedingerHermiteReal3D)
   if target.m ≠ source.m
       error("source and target must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_copy_wfs_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), target.p, source.p )
end

## method: get_nx #############################################################################

function get_nx(m::SchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_hermite_1d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_hermite_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_hermite_3d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_hermite_real_1d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_hermite_real_2d"), Int32,
         (Ptr{Void}, ), m.m )
end

function get_nx(m::SchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_get_nx_schroedinger_hermite_real_3d"), Int32,
         (Ptr{Void}, ), m.m )
end


## method: get_omega_x #############################################################################

function get_omega_x(m::SchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_get_omega_x_schroedinger_hermite_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_x(m::SchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_get_omega_x_schroedinger_hermite_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_x(m::SchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_get_omega_x_schroedinger_hermite_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_x(m::SchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_get_omega_x_schroedinger_hermite_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_x(m::SchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_get_omega_x_schroedinger_hermite_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_x(m::SchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_get_omega_x_schroedinger_hermite_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_omega_y #############################################################################

function get_omega_y(m::SchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_get_omega_y_schroedinger_hermite_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_y(m::SchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_get_omega_y_schroedinger_hermite_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_y(m::SchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_get_omega_y_schroedinger_hermite_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_y(m::SchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_get_omega_y_schroedinger_hermite_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_omega_z #############################################################################

function get_omega_z(m::SchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_get_omega_z_schroedinger_hermite_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_omega_z(m::SchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_get_omega_z_schroedinger_hermite_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! schroedinger specific methods !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## method: get_hbar #############################################################################

function get_hbar(m::SchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_hermite_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_hermite_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_hermite_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_hermite_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_hermite_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_hbar(m::SchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_get_hbar_schroedinger_hermite_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_mass #############################################################################

function get_mass(m::SchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_hermite_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_hermite_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_hermite_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_hermite_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_hermite_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_mass(m::SchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_get_mass_schroedinger_hermite_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

## method: get_cubic_coupling #############################################################################

function get_cubic_coupling(m::SchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_hermite_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_hermite_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_hermite_3d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_hermite_real_1d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_hermite_real_2d"), Float64,
         (Ptr{Void}, ), m.m )
end

function get_cubic_coupling(m::SchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_get_cubic_coupling_schroedinger_hermite_real_3d"), Float64,
         (Ptr{Void}, ), m.m )
end



## method: set_potential! ##########################################################################

function set_potential!(m::SchroedingerHermite1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_hermite_1d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerHermite2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_hermite_2d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerHermite3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerHermiteReal1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble,))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerHermiteReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

function set_potential!(m::SchroedingerHermiteReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_set_potential_schroedinger_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Void}), m.m, f_c )
end

## method: get_potential ###########################################################################

function get_potential(m::SchroedingerHermite1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_hermite_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerHermite2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_hermite_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerHermite3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_hermite_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], dims[3], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerHermiteReal1D, unsafe_access::Bool=false)
   dims =Array(Int32, 1)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_hermite_real_1d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerHermiteReal2D, unsafe_access::Bool=false)
   dims =Array(Int32, 2)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_hermite_real_2d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end

function get_potential(m::SchroedingerHermiteReal3D, unsafe_access::Bool=false)
   dims =Array(Int32, 3)
   Vp = ccall( dlsym(tssm_handle, "c_get_potential_schroedinger_hermite_real_3d"), Ptr{Float64},
         (Ptr{Void}, Ptr{Int32}), m.m, dims )
   V = pointer_to_array(Vp, dims[1], dims[2], dims[3], false)   
   if unsafe_access
      return V 
   else
      return copy(V)
   end
end


## method: load_potential! #########################################################################

function load_potential!(m::SchroedingerHermite1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_hermite_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerHermite2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_hermite_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerHermite3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerHermiteReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerHermiteReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function load_potential!(m::SchroedingerHermiteReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_load_potential_schroedinger_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

## method: save_potential #########################################################################

function save_potential(m::SchroedingerHermite1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_hermite_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerHermite2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_hermite_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerHermite3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_hermite_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerHermiteReal1D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_hermite_real_1d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerHermiteReal2D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_hermite_real_2d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

function save_potential(m::SchroedingerHermiteReal3D, filename::ASCIIString)
   ccall( dlsym(tssm_handle, "c_save_potential_schroedinger_hermite_real_3d"), Void,
         (Ptr{Void}, Ptr{Uint8}, Int32,), m.m, filename, length(filename))
end    

## method: imaginary_time_propagate_A! ############################################################

function imaginary_time_propagate_A!(psi::WfSchroedingerHermite1D, dt::Number)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerHermite2D, dt::Number)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerHermite3D, dt::Number)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64},), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerHermiteReal1D, dt::Real)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_hermite_real_real_1d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerHermiteReal2D, dt::Real)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

function imaginary_time_propagate_A!(psi::WfSchroedingerHermiteReal3D, dt::Real)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_a_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Float64,), psi.p, dt)
end

## method: imaginary_time_propagate_B! ############################################################

function imaginary_time_propagate_B!(psi::WfSchroedingerHermite1D, dt::Number, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Complex{Float64}, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerHermite2D, dt::Number, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Complex{Float64}, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerHermite3D, dt::Number, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64}, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerHermiteReal1D, dt::Real, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_hermite_real_real_1d"), Void,
                    (Ptr{Void}, Float64, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerHermiteReal2D, dt::Real, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64, Int32), psi.p, dt, method_for_B)
end

function imaginary_time_propagate_B!(psi::WfSchroedingerHermiteReal3D, dt::Real, method_for_B::Integer=0)
   ccall( dlsym(tssm_handle, "c_imaginary_time_propagate_b_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Float64, Int32), psi.p, dt, method_for_B)
end

## method: add_apply_B! ###########################################################################

function add_apply_B!(this::WfSchroedingerHermite1D, other::WfSchroedingerHermite1D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerHermite2D, other::WfSchroedingerHermite2D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerHermite3D, other::WfSchroedingerHermite3D,
                     coefficient::Number=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Complex{Float64}), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerHermiteReal1D, other::WfSchroedingerHermiteReal1D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_hermite_real_1d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerHermiteReal2D, other::WfSchroedingerHermiteReal2D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

function add_apply_B!(this::WfSchroedingerHermiteReal3D, other::WfSchroedingerHermiteReal3D,
                     coefficient::Real=1.0)
   if this.m ≠ other.m
       error("this and other must belong to the same method")
   end
   ccall( dlsym(tssm_handle, "c_add_apply_b_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Ptr{Void}, Float64), 
                     this.p, other.p, coefficient)
end

## method: kinetic_energy #########################################################################

function kinetic_energy(psi::WfSchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_hermite_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_hermite_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function kinetic_energy(psi::WfSchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_kinetic_energy_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: potential_energy #######################################################################

function potential_energy(psi::WfSchroedingerHermite1D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_hermite_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerHermite2D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_hermite_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerHermite3D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerHermiteReal1D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerHermiteReal2D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function potential_energy(psi::WfSchroedingerHermiteReal3D)
   ccall( dlsym(tssm_handle, "c_potential_energy_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: interaction_energy #####################################################################

function interaction_energy(psi::WfSchroedingerHermite1D, f::Function)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_hermite_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerHermite2D, f::Function)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_hermite_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerHermite3D, f::Function)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_hermite_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerHermiteReal1D, f::Function)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerHermiteReal2D, f::Function)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, ), psi.p )
end

function interaction_energy(psi::WfSchroedingerHermiteReal3D, f::Function)
   ccall( dlsym(tssm_handle, "c_interaction_energy_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, ), psi.p )
end

## method: observable #############################################################################

function observable(psi::WfSchroedingerHermite1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, ))
   ccall( dlsym(tssm_handle, "c_observable_wfs_hermite_1d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedingerHermite2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_hermite_2d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c)
end

function observable(psi::WfSchroedingerHermite3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_hermite_3d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedingerHermiteReal1D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, ))
   ccall( dlsym(tssm_handle, "c_observable_wfs_hermite_real_1d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedingerHermiteReal2D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_hermite_real_2d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

function observable(psi::WfSchroedingerHermiteReal3D, f::Function)
   f_c = cfunction(f, Cdouble, (Cdouble, Cdouble, Cdouble))
   ccall( dlsym(tssm_handle, "c_observable_wfs_hermite_real_3d"), Float64,
         (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
end

## method: get_energy_expectation_deviation #######################################################

function get_energy_expectation_deviation(psi::WfSchroedingerHermite1D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_hermite_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerHermite2D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_hermite_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerHermite3D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_hermite_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerHermiteReal1D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_hermite_real_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerHermiteReal2D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_hermite_real_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_energy_expectation_deviation(psi::WfSchroedingerHermiteReal3D)
   ans =Array(Float64, 2)
   ccall( dlsym(tssm_handle, "c_get_energy_expectation_deviation_wfs_hermite_real_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

## method: get_realspace_observables ##############################################################

function get_realspace_observables(psi::WfSchroedingerHermite1D)
   ans =Array(Float64, 4)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_hermite_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerHermite2D)
   ans =Array(Float64, 6)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_hermite_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerHermite3D)
   ans =Array(Float64, 8)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_hermite_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerHermiteReal1D)
   ans =Array(Float64, 4)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_hermite_real_1d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerHermiteReal2D)
   ans =Array(Float64, 6)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_hermite_real_2d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end

function get_realspace_observables(psi::WfSchroedingerHermiteReal3D)
   ans =Array(Float64, 8)
   ccall( dlsym(tssm_handle, "c_get_realspace_observables_wfs_hermite_real_3d"), Void,
        (Ptr{Void}, Ptr{Float64}), psi.p, ans )
   tuple(ans...)
end


## method: selfconsistent_nonlinear_step! ##########################################################

function selfconsistent_nonlinear_step!(psi::WfSchroedingerHermite1D, dt::Number, 
                  dt1::Number, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_hermite_1d"), Void,
                    (Ptr{Void}, Complex{Float64}, Complex{Float64}, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerHermite2D, dt::Number,
                  dt1::Number, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_hermite_2d"), Void,
                    (Ptr{Void}, Complex{Float64}, Complex{Float64}, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerHermite3D, dt::Number,
                  dt1::Number, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_hermite_3d"), Void,
                    (Ptr{Void}, Complex{Float64}, Complex{Float64}, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerHermiteReal1D, dt::Real,
                  dt1::Real, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_hermite_real_real_1d"), Void,
                    (Ptr{Void}, Float64, Float64, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerHermiteReal2D, dt::Real,
                  dt1::Real, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_hermite_real_2d"), Void,
                    (Ptr{Void}, Float64, Float64, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end

function selfconsistent_nonlinear_step!(psi::WfSchroedingerHermiteReal3D, dt::Real,
                  dt1::Real, eps::Number=100.0_prec*eps(Float64), max_iters::Integer=30)
   ccall( dlsym(tssm_handle, "c_selfconsistent_nonlinear_step_wfs_hermite_real_3d"), Void,
                    (Ptr{Void}, Float64, Float64, Float64, Int32), psi.p, dt, dt1, eps, max_itesr)
end



# end
