println("including tssm_schroedinger.jl for type ", T)
for (METHOD, SUF, COMPLEX_METHOD, DIM) in (
                 (:Schroedinger1D, :_schroedinger_1d, true, 1 ),             
                 (:Schroedinger2D, :_schroedinger_2d, true, 2 ),             
                 (:Schroedinger3D, :_schroedinger_3d, true, 3 ),             
                 (:SchroedingerReal1D, :_schroedinger_real_1d, false, 1 ), 
                 (:SchroedingerReal2D, :_schroedinger_real_2d, false, 2 ),
                 (:SchroedingerReal3D, :_schroedinger_real_3d, false, 3 ),
                 (:SchroedingerRotating2D, :_schroedinger_rotating_2d, true, 2 ),             
                 (:SchroedingerRotating3D, :_schroedinger_rotating_3d, true, 3 ),                       
                 (:SchroedingerHermite1D, :_schroedinger_hermite_1d, true, 1 ),             
                 (:SchroedingerHermite2D, :_schroedinger_hermite_2d, true, 2 ),             
                 (:SchroedingerHermite3D, :_schroedinger_hermite_3d, true, 3 ),             
                 (:SchroedingerHermiteReal1D, :_schroedinger_hermite_real_1d, false, 1 ), 
                 (:SchroedingerHermiteReal2D, :_schroedinger_hermite_real_2d, false, 2 ),
                 (:SchroedingerHermiteReal3D, :_schroedinger_hermite_real_3d, false, 3 ),
                 (:SchroedingerGeneralizedLaguerre2D, :_schroedinger_gen_laguerre_2d, true, 2 ),             
                 (:SchroedingerGeneralizedLaguerreHermite3D, :_schroedinger_gen_laguerre_hermite_3d, true, 3 ),             
                 (:SchroedingerGeneralizedLaguerreReal2D, :_schroedinger_gen_laguerre_real_2d, false, 2 ),
                 (:SchroedingerGeneralizedLaguerreHermiteReal3D, :_schroedinger_gen_laguerre_hermite_real_3d, false, 3 ),                 
                )
    println("    ", METHOD)   
    if T == :Float128
        SUF = symbol(SUF, "_wf128")
    end
    WF = symbol(:Wf,METHOD)
    @eval begin

        function get_hbar(m::($METHOD){($T)})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_hbar",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end

        function get_hbar(m::($METHOD){($T)})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_hbar",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end

        function get_cubic_coupling(m::($METHOD){($T)})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_cubic_coupling",SUF))), ($T),
                 (Ptr{Void}, ), m.m )
        end

        function set_cubic_coupling!(m::($METHOD){($T)}, cc::Real)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_cubic_coupling",SUF))), Void,
                 (Ptr{Void}, ($T)), m.m, cc )
        end

        function load_potential!(m::($METHOD){($T)}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"load_potential",SUF))), Void,
                 (Ptr{Void}, Ptr{UInt8}, Int32,), m.m, filename, length(filename))
        end    

        function save_potential(m::($METHOD){($T)}, filename::ASCIIString)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"save_potential",SUF))), Void,
                 (Ptr{Void}, Ptr{UInt8}, Int32,), m.m, filename, length(filename))
        end

        function kinetic_energy(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"kinetic_energy_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function potential_energy(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"potential_energy_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function interaction_energy(psi::($WF){$T})
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"interaction_energy_wf",SUF))), ($T),
                 (Ptr{Void}, ), psi.p )
        end

        function get_energy_expectation_deviation(psi::($WF){$T})
           ans =Array(($T), 2)
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_energy_expectation_deviation_wf",SUF))), Void,
                (Ptr{Void}, Ptr{($T)}), psi.p, ans )
           tuple(ans...)
        end

        function get_realspace_observables(psi::($WF){$T})
           ans =Array(($T), 2+2*($DIM))
           ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_realspace_observables_wf",SUF))), Void,
                (Ptr{Void}, Ptr{($T)}), psi.p, ans )
           tuple(ans...)
        end

    end #eval

    if COMPLEX_METHOD
        @eval begin
            function propagate_B!(psi::($WF){$T}, dt::Number)
                ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_B_wf",SUF))), Void,
                        (Ptr{Void}, Complex{$T},), psi.p, dt)
            end

            function propagate_B_derivative!(this::($WF){$T}, other::($WF){$T}, dt::Number)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"propagate_B_derivative_wf",SUF))), Void,
                                (Ptr{Void}, Ptr{Void}, Complex{($T)}), 
                                 this.p, other.p, dt)
            end

            function imaginary_time_propagate_A!(psi::($WF){$T}, dt::Number)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"imaginary_time_propagate_A_wf",SUF))), Void,
                                (Ptr{Void}, Complex{($T)},), psi.p, dt)
            end

            function imaginary_time_propagate_B!(psi::($WF){$T}, dt::Number, method_for_B::Integer=0)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"imaginary_time_propagate_B_wf",SUF))), Void,
                                (Ptr{Void}, Complex{($T)}, Int32), psi.p, dt, method_for_B)
            end

            function add_apply_B!(this::($WF){$T}, other::($WF){$T},
                                 coefficient::Number=1.0)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"add_apply_B_wf",SUF))), Void,
                                (Ptr{Void}, Ptr{Void}, Complex{($T)}), 
                                 this.p, other.p, coefficient)
            end

            function selfconsistent_nonlinear_step!(psi::($WF){$T}, dt::Number, 
                           dt1::Number, eps::Number=100.0_prec*eps(($T)), max_iters::Integer=30)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"selfconsistent_nonlinear_step_wf",SUF))), Void,
                          (Ptr{Void}, Complex{($T)}, Complex{($T)}, ($T), Int32), psi.p, dt, dt1, eps, max_itesr)
            end
            
         end #eval   
    else
        @eval begin
            function imaginary_time_propagate_A!(psi::($WF){$T}, dt::Real)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"imaginary_time_propagate_A_wf",SUF))), Void,
                                (Ptr{Void}, ($T),), psi.p, dt)
            end

            function imaginary_time_propagate_B!(psi::($WF){$T}, dt::Real, method_for_B::Integer=0)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"imaginary_time_propagate_B_wf",SUF))), Void,
                                (Ptr{Void}, ($T), Int32), psi.p, dt, method_for_B)
            end

            function add_apply_B!(this::($WF){$T}, other::($WF){$T},
                                 coefficient::Number=1.0)
               if this.m ≠ other.m
                   error("this and other must belong to the same method")
               end
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"add_apply_B_wf",SUF))), Void,
                                (Ptr{Void}, Ptr{Void}, ($T)), 
                                 this.p, other.p, coefficient)
            end

            function selfconsistent_nonlinear_step!(psi::($WF){$T}, dt::Number, 
                           dt1::Number, eps::Number=100.0_prec*eps(($T)), max_iters::Integer=30)
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"selfconsistent_nonlinear_step_wf",SUF))), Void,
                          (Ptr{Void}, ($T), ($T), ($T), Int32), psi.p, dt, dt1, eps, max_itesr)
            end
        end #eval   
    end #if

    if DIM==1
        @eval begin
            function set_potential!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T),))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}), m.m, f_c )
            end

 
            function get_potential(m::($METHOD){($T)}, unsafe_access::Bool=false)
               dims =Array(Int32, 1)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential",SUF))), Ptr{($T)},
                     (Ptr{Void}, Ptr{Int32}), m.m, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1])
               end
               V = pointer_to_array(Vp, dims[1], false)   
               if unsafe_access
                  return V 
               else
                  return copy(V)
               end
            end

            function observable(psi::($WF){$T}, f::Function)
               f_c = cfunction(f, ($T), (($T), ))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"observable_wf",SUF))), ($T),
                     (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
            end

        end #eval   
        if COMPLEX_METHOD
        @eval begin
            function set_potential_t!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T),($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential_t",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}, Bool), m.m, f_c, false )
            end

            function set_potential_t_derivative!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T),($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential_t",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}, Bool), m.m, f_c, true )
            end

            function get_potential_t(m::($METHOD){($T)}, t::Real)
               dims =Array(Int32, 1)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential_t",SUF))), Ptr{($T)},
                     (Ptr{Void}, ($T), Bool, Ptr{Int32}), m.m, t, false, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1])
               end
               V = pointer_to_array(Vp, dims[1], false)   
               return copy(V)
            end

            function get_potential_t_derivative(m::($METHOD){($T)}, t::Real)
               dims =Array(Int32, 1)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential_t",SUF))), Ptr{($T)},
                     (Ptr{Void}, ($T), Bool, Ptr{Int32}), m.m, t, true, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1])
               end
               V = pointer_to_array(Vp, dims[1], false)   
               return copy(V)
            end

        end #eval
        end #if
    elseif DIM==2
        @eval begin
            function set_potential!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T), ($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}), m.m, f_c )
            end
            
            function get_potential(m::($METHOD){($T)}, unsafe_access::Bool=false)
               dims =Array(Int32, 2)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential",SUF))), Ptr{($T)},
                     (Ptr{Void}, Ptr{Int32}), m.m, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1], dims[2])
               end
               V = pointer_to_array(Vp, dims[1], dims[2], false)   
               if unsafe_access
                  return V 
               else
                  return copy(V)
               end
            end

            function observable(psi::($WF){$T}, f::Function)
               f_c = cfunction(f, ($T), (($T), ($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"observable_wf",SUF))), ($T),
                     (Ptr{Void}, Ptr{Void} ), psi.p, f_c)
            end

        end #eval   
        if COMPLEX_METHOD
        @eval begin
            function set_potential_t!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T),($T),($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential_t",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}, Bool), m.m, f_c, false )
            end

            function set_potential_t_derivative!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T),($T),($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential_t",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}, Bool), m.m, f_c, true )
            end

            function get_potential_t(m::($METHOD){($T)}, t::Real)
               dims =Array(Int32, 2)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential_t",SUF))), Ptr{($T)},
                     (Ptr{Void}, ($T), Bool, Ptr{Int32}), m.m, t, false, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1], dims[2])
               end
               V = pointer_to_array(Vp, dims[2], false)   
               return copy(V)
            end

            function get_potential_t_derivative(m::($METHOD){($T)}, t::Real)
               dims =Array(Int32, 2)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential_t",SUF))), Ptr{($T)},
                     (Ptr{Void}, ($T), Bool, Ptr{Int32}), m.m, t, true, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1], dims[2])
               end
               V = pointer_to_array(Vp, dims[2], false)   
               return copy(V)
            end
            
            
        end #eval
        end #if
    elseif DIM==3
        @eval begin
            function set_potential!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T), ($T), ($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}), m.m, f_c )
            end

            function get_potential(m::($METHOD){($T)}, unsafe_access::Bool=false)
               dims =Array(Int32, 3)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential",SUF))), Ptr{($T)},
                     (Ptr{Void}, Ptr{Int32}), m.m, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1], dims[2], dims[3])
               end
               V = pointer_to_array(Vp, dims[1], dims[2], dims[3], false)   
               if unsafe_access
                  return V 
               else
                  return copy(V)
               end
            end

            function observable(psi::($WF){$T}, f::Function)
               f_c = cfunction(f, ($T), (($T), ($T), ($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"observable_wf",SUF))), ($T),
                     (Ptr{Void}, Ptr{Void} ), psi.p, f_c )
            end

        end #eval   
        if COMPLEX_METHOD
        @eval begin
            function set_potential_t!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T),($T),($T)($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential_t",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}, Bool), m.m, f_c, false )
            end

            function set_potential_t_derivative!(m::($METHOD){($T)}, f::Function)
               f_c = cfunction(f, ($T), (($T),($T),($T)($T)))
               ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"set_potential_t",SUF))), Void,
                     (Ptr{Void}, Ptr{Void}, Bool), m.m, f_c, true )
            end

            function get_potential_t(m::($METHOD){($T)}, t::Real)
               dims =Array(Int32, 3)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential_t",SUF))), Ptr{($T)},
                     (Ptr{Void}, ($T), Bool, Ptr{Int32}), m.m, t, false, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1], dims[2], dims[3])
               end
               V = pointer_to_array(Vp, dims[3], false)   
               return copy(V)
            end 

            function get_potential_t_derivative(m::($METHOD){($T)}, t::Real)
               dims =Array(Int32, 3)
               Vp = ccall( Libdl.dlsym(($TSSM_HANDLE), $(string(PRE,"get_potential_t",SUF))), Ptr{($T)},
                     (Ptr{Void}, ($T), Bool, Ptr{Int32}), m.m, t, true, dims )
               if Vp==C_NULL
                   return zeros(($T), dims[1], dims[2], dims[3])
               end
               V = pointer_to_array(Vp, dims[3], false)   
               return copy(V)
            end 
            

        end #eval
        end #if
    end #if

end # for 

