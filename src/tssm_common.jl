println("including tssm_common.jl for type ", T)    

function clone(psi::WaveFunction)
    wave_function(psi.m)
end

function eigen_function!(psi::WaveFunction1D, k)
    to_frequency_space!(psi)
    u=get_data(psi,true)
    u[:]=0.0
    u[k]=1
    psi
end

function eigen_function!(psi::WaveFunction2D, k, l)
    to_frequency_space!(psi)
    u=get_data(psi,true)
    u[:,:]=0.0
    u[k,l]=1
    psi
end

function eigen_function!(psi::WaveFunction3D, k, l, m)
    to_frequency_space!(psi)
    u=get_data(psi,true)
    u[:,:,:]=0.0
    u[k,l,m]=1
    psi
end

#NONSEPARATED_EIGENVALUES 0=cartesian, 1=polar/cylindrical, 2=spherical

for (METHOD, SUF, COMPLEX_METHOD, DIM, NONSEPARATED_EIGENVALUES) in (
                 (:Fourier1D, :_fourier_1d, true, 1, 0 ),             
                 (:Fourier2D, :_fourier_2d, true, 2, 0 ),             
                 (:Fourier3D, :_fourier_3d, true, 3, 0 ),             
                 (:FourierReal1D, :_fourier_real_1d, false, 1, 0 ), 
                 (:FourierReal2D, :_fourier_real_2d, false, 2, 0 ),
                 (:FourierReal3D, :_fourier_real_3d, false, 3, 0 ),
                 (:FourierBessel2D, :_fourier_bessel_2d, true, 2, 1 ),             
                 (:FourierBesselReal2D, :_fourier_bessel_real_2d, false, 2, 1 ),
                 (:BesselRotSym1D, :_bessel_rotsym_1d, true, 1,0),             
                 (:BesselRotSymReal1D, :_bessel_rotsym_real_1d, false, 1, 0),
                 (:Schroedinger1D, :_schroedinger_1d, true, 1, 0 ),             
                 (:Schroedinger2D, :_schroedinger_2d, true, 2, 0 ),             
                 (:Schroedinger3D, :_schroedinger_3d, true, 3, 0 ),             
                 (:SchroedingerReal1D, :_schroedinger_real_1d, false, 1, 0 ), 
                 (:SchroedingerReal2D, :_schroedinger_real_2d, false, 2, 0 ),
                 (:SchroedingerReal3D, :_schroedinger_real_3d, false, 3, 0 ),
                 (:SchroedingerRotating2D, :_schroedinger_rotating_2d, true, 2, 0 ),             
                 (:SchroedingerRotating3D, :_schroedinger_rotating_3d, true, 3, 0 ),                           
                 (:SchroedingerHermite1D, :_schroedinger_hermite_1d, true, 1, 0 ),             
                 (:SchroedingerHermite2D, :_schroedinger_hermite_2d, true, 2, 0 ),             
                 (:SchroedingerHermite3D, :_schroedinger_hermite_3d, true, 3, 0 ),             
                 (:SchroedingerHermiteReal1D, :_schroedinger_hermite_real_1d, false, 1, 0 ), 
                 (:SchroedingerHermiteReal2D, :_schroedinger_hermite_real_2d, false, 2, 0 ),
                 (:SchroedingerHermiteReal3D, :_schroedinger_hermite_real_3d, false, 3, 0 ),
                 (:SchroedingerGeneralizedLaguerre2D, :_schroedinger_gen_laguerre_2d, true, 2, 0 ),             
                 (:SchroedingerGeneralizedLaguerreHermite3D, :_schroedinger_gen_laguerre_hermite_3d, true, 3, 0 ),             
                 (:SchroedingerGeneralizedLaguerreReal2D, :_schroedinger_gen_laguerre_real_2d, false, 2, 0 ),
                 (:SchroedingerGeneralizedLaguerreHermiteReal3D, :_schroedinger_gen_laguerre_hermite_real_3d, false, 3, 0 ),
                )
    println("    ", METHOD)
    if T == :Float128
        SUF = symbol(SUF, "_wf128")
    end
    WF = symbol(:Wf,METHOD)
    @eval begin
        # wave function constructor

        #function ($WF){($T)}( m::($METHOD){($T)} )
        function ($WF)( m::($METHOD){($T)} )
            wf = ($WF){($T)}( ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"new_wf",SUF))), Ptr{Void}, (Ptr{Void},), m.m) , m)   
            finalizer(wf, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"finalize_wf",SUF))), Void, (Ptr{Void},), x.p) )
            wf
        end
     end
     @eval begin

        wave_function(m::($METHOD){($T)}) = ($WF)(m) 

        function set_propagate_time_together_with_A!(m::($METHOD){$T}, flag::Bool)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_propagate_time_together_with_A",SUF))), Void,
               (Ptr{Void}, Int32), m.m, flag) 
        end

        function get_propagate_time_together_with_A(m::($METHOD){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_propagate_time_together_with_A",SUF))), Int32,
               (Ptr{Void},), m.m) == 1
        end

        function is_real_space(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_wf",SUF))), Int32,
               (Ptr{Void},), psi.p) == 1
        end

        function is_frequency_space(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"is_real_space_wf",SUF))), Int32,
               (Ptr{Void},), psi.p) != 1
        end

        function to_real_space!(psi::($WF){$T})
             ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_real_space_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end

        function to_frequency_space!(psi::($WF){$T})
             ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"to_frequency_space_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end

        function save(psi::($WF){$T}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"save_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{UInt8}, Int32,), psi.p, filename, length(filename))
        end         

        function load!(psi::($WF){$T}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"load_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{UInt8}, Int32,), psi.p, filename, length(filename))
        end    

        function copy!(target::($WF){$T}, source::($WF){$T})
           if target.m ≠ source.m
               error("source and target must belong to the same method")
           end
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"copy_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{Void}), target.p, source.p )
        end
        
        function norm(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"norm_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function norm_in_frequency_space(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"norm_in_frequency_space_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function distance(psi1::($WF){$T}, psi2::($WF){$T})
           if psi1.m ≠ psi2.m
               error("psi1 and psi2 must belong to the same method")
           end
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"distance_wf",SUF))), ($T),
                 (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
        end

        function normalize!(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"normalize_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function get_eigenvalues(m::($METHOD){$T}, unsafe_access::Bool=false)
           dim =Array(Int32, 1)
           evp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                 (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
           ev = pointer_to_array(evp, dim[1], false)   
           if unsafe_access
               return ev
           else    
               return copy(ev)
           end
        end

    end #eval

    if NONSEPARATED_EIGENVALUES==0 # cartesian
        if DIM==1
            @eval begin
                function get_eigenvalues(m::($METHOD){$T}, unsafe_access::Bool=false)
                   dim =Array(Int32, 1)
                   evp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
                   ev = pointer_to_array(evp, dim[1], false)   
                   if unsafe_access
                       return ev
                   else    
                       return copy(ev)
                   end
                end

            end
        elseif DIM==2
            @eval begin
               function get_eigenvalues(m::($METHOD){$T}, unsafe_access::Bool=false)
                   dim =Array(Int32, 1)
                   evp1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
                   ev1 = pointer_to_array(evp1, dim[1], false)   
                   evp2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
                   ev2 = pointer_to_array(evp2, dim[1], false)   
                   if unsafe_access
                       return ev1, ev2
                   else    
                       return copy(ev1), copy(ev2)
                   end
                end

            end
        elseif DIM==3
            @eval begin
               function get_eigenvalues(m::($METHOD){$T}, unsafe_access::Bool=false)
                   dim =Array(Int32, 1)
                   evp1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
                   ev1 = pointer_to_array(evp1, dim[1], false)   
                   evp2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
                   ev2 = pointer_to_array(evp2, dim[1], false)   
                   evp3 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
                   ev3 = pointer_to_array(evp3, dim[1], false)   
                   if unsafe_access
                       return ev1, ev2, ev3
                   else    
                       return copy(ev1), copy(ev2), copy(ev3)
                   end
                end

            end
        end
    elseif NONSEPARATED_EIGENVALUES==1 # polar/cylindrical
        if DIM==2
            @eval begin
                function get_eigenvalues(m::($METHOD){$T}, unsafe_access::Bool=false)
                   dim =Array(Int32, 2)
                   evp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_eigenvalues",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}), m.m, dim )
                   ev = pointer_to_array(evp, (dim[1], dim[2]), false)   
                   if unsafe_access
                       return ev
                   else    
                       return copy(ev)
                   end
                end
            
            end
        elseif DIM==3
            # tbd
        end
    elseif NONSEPARATED_EIGENVALUES==2 # spherical 
        if DIM==3
            # tbd
        end
    end


    if DIM==1
        @eval begin

            function get_nodes(m::($METHOD){$T})
               dim =Array(Int32, 1)
               np = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
               n = pointer_to_array(np, dim[1], false)     
               copy(n)
            end

        end  
    elseif DIM==2
        @eval begin

            function get_nodes(m::($METHOD){$T})
               dim =Array(Int32, 1)
               np1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
               n1 = pointer_to_array(np1, dim[1], false)     
               np2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
               n2 = pointer_to_array(np2, dim[1], false)     
               copy(n1), copy(n2)
            end

        end  
    elseif DIM==3
        @eval begin

            function get_nodes(m::($METHOD){$T})
               dim =Array(Int32, 1)
               np1 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 1 )
               n1 = pointer_to_array(np1, dim[1], false)     
               np2 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 2 )
               n2 = pointer_to_array(np2, dim[1], false)     
               np3 = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_nodes",SUF))), Ptr{$T},
                     (Ptr{Void}, Ptr{Int32}, Int32), m.m, dim, 3 )
               n3 = pointer_to_array(np3, dim[1], false)     
               copy(n1), copy(n2), copy(n3)
            end

          end  
     end


    if COMPLEX_METHOD
        @eval begin

            function save(psi::($WF){$T}, filename::ASCIIString, 
                         dset_name_real::ASCIIString, dset_name_imag::ASCIIString;
                         append::Bool=false)
   	         ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"save1_wf",SUF))), Void,
                 	(Ptr{Void}, Ptr{UInt8}, Int32, Ptr{UInt8}, Int32, Ptr{UInt8}, Int32, Int32), 
	                  psi.p, filename, length(filename),
        	          dset_name_real, length(dset_name_real),
                	  dset_name_imag, length(dset_name_imag), append)
            end         

            function load!(psi::($WF){$T}, filename::ASCIIString, 
                         dset_name_real::ASCIIString, dset_name_imag::ASCIIString)
   	         ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"load1_wf",SUF))), Void,
                 	(Ptr{Void}, Ptr{UInt8}, Int32, Ptr{UInt8}, Int32, Ptr{UInt8}, Int32), 
	                  psi.p, filename, length(filename),
        	          dset_name_real, length(dset_name_real),
                	  dset_name_imag, length(dset_name_imag))
            end         

            function set_time!(psi::($WF){$T}, t::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_time_wf",SUF))), Void,
                        (Ptr{Void}, Complex{$T},), psi.p, t)
            end

            function get_time(psi::($WF){$T})
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_time_wf",SUF))), ($T),
                        (Ptr{Void},), psi.p)
            end

            function propagate_time!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_time_wf",SUF))), Void,
                        (Ptr{Void}, Complex{$T},), psi.p, dt)
            end

            function propagate_A!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_A_wf",SUF))), Void,
                        (Ptr{Void}, Complex{$T},), psi.p, dt)
            end

            function propagate_A_derivative!(this::($WF){$T}, other::($WF){$T},
                                 dt::Number)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_A_derivative_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, Complex{$T}), 
                         this.p, other.p, dt)
            end
    
            function add_apply_A!(this::($WF){$T}, other::($WF){$T},
                                 coefficient::Number=1.0)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"add_apply_A_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, Complex{$T}), 
                         this.p, other.p, coefficient)
            end
    
            function scale!(psi::($WF){$T}, factor::Number)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"scale_wf",SUF))), Void,
                     (Ptr{Void}, Complex{$T} ), psi.p, factor )
            end

            function axpy!(this::($WF){$T}, other::($WF){$T},
                                 factor::Number)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"axpy_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, Complex{$T}), 
                         this.p, other.p, factor)
            end

            function inner_product(psi1::($WF){$T}, psi2::($WF){$T})
                if psi1.m ≠ psi2.m
                    error("psi1 and psi2 must belong to the same method")
                end
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"inner_product_wf",SUF))), (Complex{$T}),
                    (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
            end

        end # eval    
            

        if DIM==1
            @eval begin
                function set!(psi::($WF){$T}, x::Number)
                    u = get_data(psi, true)
                    u[:] = x
                end

                function set!(psi::($WF){$T}, f::Function)
                   try
                       f_c = cfunction(f, ($T), (($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"rset_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   end      
                end
        
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   try
                       f_c = cfunction(f, ($T), (($T), ($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"rset_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T), ($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   end      
                end
            
                function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 1)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_data_wf",SUF))), Ptr{Complex{$T}},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, dims[1], false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

                function evaluate(psi::($WF){$T}, x::Real)
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"evaluate_wf",SUF))), Complex{$T},
                             (Ptr{Void}, ($T)), psi.p, x)
                end
                

            end # eval    
        elseif DIM==2
            @eval begin

                function set!(psi::($WF){$T}, x::Number)
                    u = get_data(psi, true)
                    u[:,:] = x
                end

                function set!(psi::($WF){$T}, f::Function)
                   try
                       f_c = cfunction(f, ($T), (($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"rset_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   end      
                end
        
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   try
                       f_c = cfunction(f, ($T), (($T), ($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"rset_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T), ($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   end      
                end

                function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 2)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_data_wf",SUF))), Ptr{Complex{$T}},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

                function evaluate(psi::($WF){$T}, x::Real, y::Real)
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"evaluate_wf",SUF))), Complex{$T},
                             (Ptr{Void}, ($T),($T)), psi.p, x, y)
                end

            end # eval    
        elseif DIM==3
            @eval begin
                function set!(psi::($WF){$T}, x::Number)
                    u = get_data(psi, true)
                    u[:,:,:] = x
                end

                function set!(psi::($WF){$T}, f::Function)
                   try
                       f_c = cfunction(f, ($T), (($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"rset_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   end      
                end
        
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   try
                       f_c = cfunction(f, ($T), (($T), ($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"rset_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T), ($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   end      
                end

               function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 3)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_data_wf",SUF))), Ptr{Complex{$T}},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

                function evaluate(psi::($WF){$T}, x::Real, y::Real, z::Real)
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"evaluate_wf",SUF))), Complex{$T},
                             (Ptr{Void}, ($T),($T),($T)), psi.p, x, y,z)
                end
            
            end # eval    
        end

    else    
        @eval begin

            function save(psi::($WF){$T}, filename::ASCIIString, 
                         dset_name::ASCIIString;
			 append::Bool=false)
   	         ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"save1_wf",SUF))), Void,
                 	(Ptr{Void}, Ptr{UInt8}, Int32, Ptr{UInt8}, Int32, Int32), 
	                  psi.p, filename, length(filename), dset_name, length(dset_name))
            end         

            function load!(psi::($WF){$T}, filename::ASCIIString, 
                         dset_name::ASCIIString)
   	         ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"load1_wf",SUF))), Void,
                 	(Ptr{Void}, Ptr{UInt8}, Int32, Ptr{UInt8}, Int32), 
	                  psi.p, filename, length(filename), dset_name, length(dset_name))
            end         

            function set_time!(psi::($WF){$T}, t::Real)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_time_wf",SUF))), Void,
                        (Ptr{Void}, $T,), psi.p, t)
            end

            function get_time(psi::($WF){$T}, t::Real)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_time_wf",SUF))), ($T),
                        (Ptr{Void}, $T,), psi.p, t)
            end

            function propagate_time!(psi::($WF){$T}, dt::Real)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_time_wf",SUF))), Void,
                        (Ptr{Void}, $T,), psi.p, dt)
            end
        
            function propagate_A!(psi::($WF){$T}, dt::Real)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_A_wf",SUF))), Void,
                        (Ptr{Void}, ($T),), psi.p, dt)
            end

            function propagate_A_derivative!(this::($WF){$T}, other::($WF){$T},
                                 dt::Real)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_A_derivative_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, ($T)), 
                         this.p, other.p, dt)
            end
    
            function add_apply_A!(this::($WF){$T}, other::($WF){$T},
                                 coefficient::Real=1.0)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"add_apply_A_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, ($T)), 
                         this.p, other.p, coefficient)
            end

            function scale!(psi::($WF){$T}, factor::Real)
               ccall(  Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"scale_wf",SUF))) , Void,
                     (Ptr{Void}, ($T) ), psi.p, factor )
            end
    
            function axpy!(this::($WF){$T}, other::($WF){$T},
                                 factor::Real)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"axpy_wf",SUF))), Void,
                        (Ptr{Void}, Ptr{Void}, ($T)), 
                         this.p, other.p, factor)
            end

            function inner_product(psi1::($WF){$T}, psi2::($WF){$T})
                if psi1.m ≠ psi2.m
                    error("psi1 and psi2 must belong to the same method")
                end
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"inner_product_wf",SUF))), ($T),
                    (Ptr{Void}, Ptr{Void}), psi1.p, psi2.p )
            end

        end # eval    

        if DIM==1
            @eval begin

                function set!(psi::($WF){$T}, x::Real)
                    u = get_data(psi, true)
                    u[:] = x
                end

                function set!(psi::($WF){$T}, f::Function)
                   f_c = cfunction(f, ($T), (($T),))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                end
    
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   f_c = cfunction(f, ($T), (($T), ($T),))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                end

                function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 1)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_data_wf",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, dims[1], false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

                function evaluate(psi::($WF){$T}, x::Real)
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"evaluate_wf",SUF))), ($T),
                             (Ptr{Void}, ($T)), psi.p, x)
                end

            end # eval    
        elseif DIM==2
            @eval begin
                function set!(psi::($WF){$T}, x::Real)
                    u = get_data(psi, true)
                    u[:,:] = x
                end

                function set!(psi::($WF){$T}, f::Function)
                   f_c = cfunction(f, ($T), (($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                end
    
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   f_c = cfunction(f, ($T), (($T), ($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                end

                function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 2)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_data_wf",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

                function evaluate(psi::($WF){$T}, x::Real, y::Real)
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"evaluate_wf",SUF))), ($T),
                             (Ptr{Void}, ($T),($T)), psi.p, x, y)
                end

            end # eval    
        elseif DIM==3
            @eval begin
                function set!(psi::($WF){$T}, x::Real)
                    u = get_data(psi, true)
                    u[:,:,:] = x
                end

                function set!(psi::($WF){$T}, f::Function)
                   f_c = cfunction(f, ($T), (($T),($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                end
    
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   f_c = cfunction(f, ($T), (($T), ($T),($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                end
            
                function get_data(psi::($WF){$T}, unsafe_access::Bool=false)
                   dims =Array(Int32, 3)
                   up = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_data_wf",SUF))), Ptr{$T},
                         (Ptr{Void}, Ptr{Int32}), psi.p, dims )
                   data = pointer_to_array(up, (dims[1], dims[2], dims[3]), false)     
                   if unsafe_access
                      return data
                   else
                      return copy(data)
                   end
                end

                function evaluate(psi::($WF){$T}, x::Real, y::Real, z::Real)
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"evaluate_wf",SUF))), ($T),
                             (Ptr{Void}, ($T),($T),($T)), psi.p, x, y, z)
                end
            
            end # eval    
        end
    end        

end #for 

