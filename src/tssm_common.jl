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

for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:Fourier1D, :_fourier_1d, true, 1 ),             
                 (:Fourier2D, :_fourier_2d, true, 2 ),             
                 (:Fourier3D, :_fourier_3d, true, 3 ),             
                 (:FourierReal1D, :_fourier_real_1d, false, 1 ), 
                 (:FourierReal2D, :_fourier_real_2d, false, 2 ),
                 (:FourierReal3D, :_fourier_real_3d, false, 3 ),
                 (:Schroedinger1D, :_schroedinger_1d, true, 1 ),             
                 (:Schroedinger2D, :_schroedinger_2d, true, 2 ),             
                 (:Schroedinger3D, :_schroedinger_3d, true, 3 ),             
                 (:SchroedingerReal1D, :_schroedinger_real_1d, false, 1 ), 
                 (:SchroedingerReal2D, :_schroedinger_real_2d, false, 2 ),
                 (:SchroedingerReal3D, :_schroedinger_real_3d, false, 3 ),
                 (:SchroedingerHermite1D, :_schroedinger_hermite_1d, true, 1 ),             
                 (:SchroedingerHermite2D, :_schroedinger_hermite_2d, true, 2 ),             
                 (:SchroedingerHermite3D, :_schroedinger_hermite_3d, true, 3 ),             
                 (:SchroedingerHermiteReal1D, :_schroedinger_hermite_real_1d, false, 1 ), 
                 (:SchroedingerHermiteReal2D, :_schroedinger_hermite_real_2d, false, 2 ),
                 (:SchroedingerHermiteReal3D, :_schroedinger_hermite_real_3d, false, 3 ),
                 
                )
    println("    ", METHOD)    
    WF = symbol(:Wf,METHOD)
    @eval begin
        # wave function constructor

        function ($WF){($T)}( m::($METHOD){($T)} )
            wf = ($WF){($T)}( ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_new_wf",SUF))), Ptr{Void}, (Ptr{Void},), m.m) , m)   
            finalizer(wf, x -> ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_finalize_wf",SUF))), Void, (Ptr{Ptr{Void}},), &x.p) )
            wf
        end
     end
     @eval begin

        wave_function(m::($METHOD){($T)}) = ($WF)(m) 

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

        function to_frequency_space!(psi::($WF){$T})
             ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_to_frequency_space_wf",SUF))), Void,
                    (Ptr{Void},), psi.p)
        end

        
        function save(psi::($WF){$T}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_save_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{UInt8}, Int32,), psi.p, filename, length(filename))
        end         

        function load!(psi::($WF){$T}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_load_wf",SUF))), Void,
                 (Ptr{Void}, Ptr{UInt8}, Int32,), psi.p, filename, length(filename))
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
               dim =Array(Int32, 1)
               np = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_get_nodes",SUF))), Ptr{$T},
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

        end # eval    
            

        if DIM==1
            @eval begin

                function set!(psi::($WF){$T}, f::Function)
                   try
                       f_c = cfunction(f, ($T), (($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   end      
                end
        
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   try
                       f_c = cfunction(f, ($T), (($T), ($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T), ($T),))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   end      
                end
            
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

                function set!(psi::($WF){$T}, f::Function)
                   try
                       f_c = cfunction(f, ($T), (($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   end      
                end
        
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   try
                       f_c = cfunction(f, ($T), (($T), ($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T), ($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   end      
                end

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

                function set!(psi::($WF){$T}, f::Function)
                   try
                       f_c = cfunction(f, ($T), (($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                   end      
                end
        
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   try
                       f_c = cfunction(f, ($T), (($T), ($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_rset_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   catch
                       f_c = cfunction(f, Complex{$T}, (($T), ($T),($T),($T)))
                       ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                             (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                   end      
                end

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
                        (Ptr{Void}, ($T),), psi.p, dt)
            end
    
            function propagate_B!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_b_wf",SUF))), Void,
                        (Ptr{Void}, ($T),), psi.p, dt)
            end
    
            function propagate_C!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_propagate_c_wf",SUF))), Void,
                        (Ptr{Void}, ($T),), psi.p, dt)
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
               ccall(  Libdl.dlsym(($TSSM_HANDLE), $(string(:c_scale_wf,SUF))) , Void,
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

        end # eval    

        if DIM==1
            @eval begin

                function set!(psi::($WF){$T}, f::Function)
                   f_c = cfunction(f, ($T), (($T),))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                end
    
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   f_c = cfunction(f, ($T), (($T), ($T),))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                end

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

                function set!(psi::($WF){$T}, f::Function)
                   f_c = cfunction(f, ($T), (($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                end
    
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   f_c = cfunction(f, ($T), (($T), ($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                end

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

                function set!(psi::($WF){$T}, f::Function)
                   f_c = cfunction(f, ($T), (($T),($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}), psi.p, f_c )
                end
    
                function set!(psi::($WF){$T}, f::Function, t::Real)
                   f_c = cfunction(f, ($T), (($T), ($T),($T),($T)))
                   ccall( Libdl.dlsym(($TSSM_HANDLE), $(string("c_set_t_wf",SUF))), Void,
                         (Ptr{Void}, Ptr{Void}, ($T)), psi.p, f_c, t )
                end
            
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

end #for 

