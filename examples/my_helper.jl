function equidistant_time_stepper(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, es::EmbeddedScheme, 
                         operator_sequence="AB")
    equidistant_time_stepper(psi, t0, tend, dt, es.scheme2,
                         operator_sequence)
end 
function equidistant_time_stepper(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, es::DefectBasedScheme, 
                         operator_sequence="AB")
    equidistant_time_stepper(psi, t0, tend, dt, es.scheme,
                         operator_sequence)
end 
function equidistant_time_stepper(psi::WaveFunction, t0::Real, tend::Real, 
                         dt::Real, ps::PalindromicScheme, 
                         operator_sequence="AB")
    equidistant_time_stepper(psi, t0, tend, dt, ps.scheme,
                         operator_sequence)
end 
function step!(psi::WaveFunction, dt::Real, scheme::EmbeddedScheme, operator_sequence="AB")
    step!(psi,dt,scheme.scheme2,operator_sequence)
end
function step!(psi::WaveFunction, dt::Real, scheme::PalindromicScheme, operator_sequence="AB")
    step!(psi,dt,scheme.scheme,operator_sequence)
end
function step!(psi::WaveFunction, dt::Real, scheme::DefectBasedScheme, operator_sequence="AB")
    step!(psi,dt,scheme.scheme,operator_sequence)
end

function savefig!(scheme,schemes,mytime,steps_a,nsteps_a)
    figure(1)
    hold(false)
            plot(mytime[1:nsteps_a-1], steps_a[1:nsteps_a-1])
            xlabel("t")
            ylabel("stepsize")
            if scheme==schemes[1]
                savefig("Emb43AKp_$(tol)_step.png", bbox_inches="tight")
            end
            if scheme==schemes[2]
                savefig("PP34A$(tol)_step.png", bbox_inches="tight")
            end
            if scheme==schemes[3]
                savefig("PP56A_$(tol)_step.png", bbox_inches="tight")
            end
            if scheme==schemes[4]
                savefig("Emb43AKp_defect_$(tol)_step.png", bbox_inches="tight")
            end
            if scheme==schemes[5]
                savefig("PP34A_defect$(tol)_step.png", bbox_inches="tight")
            end
            if scheme==schemes[6]
                savefig("PP56A_defect_$(tol)_step.png", bbox_inches="tight")
            end
end
function calc_stepmin(psi_ref,psia,psi,scheme,soliton)
    h=1.0
    chg=h/2
    #@printf("||psi_a-psi_ref||=%5.3E ,h=%5.3E \n",distance(psia,psi_ref),h)
    set!(psi, soliton) 
    while chg>1e-3 || distance(psia,psi_ref) < distance(psi,psi_ref)
        set!(psi, soliton) 
        time = @elapsed for t in equidistant_time_stepper(psi, t0, tend, h, scheme, "AB")
        end
        #@printf("||psi-psi_ref||=%5.3E,  h=%5.3E \n",distance(psi,psi_ref),h)
        if distance(psia,psi_ref) < distance(psi,psi_ref)
            h=h-chg
        else
            h=h+chg
        end
        chg=chg/2
    end
    #@printf("||psi-psi_ref||=%5.3E,  h=%20.18E \n",distance(psi,psi_ref),h)
    return h
end
function latextable(scheme,schemes,tol,nsteps_a,nsteps,time_a,time)
            if scheme==schemes[1]
            @printf("\\texttt{Emb 4/3 AK p}, TOL = \$ 10^{%d} \$ & \$ %d \$ & \$ %d \$ & \$ %5.3f \$ & \$ %5.3f \$ \\\\ \n",tol,nsteps_a,nsteps,time_a,time)
        end
        if scheme==schemes[2]
            @printf(    "\\texttt{PP 3/4 A}, TOL = \$ 10^{%d} \$ & \$ %d \$ & \$ %d \$ & \$ %5.3f \$ & \$ %5.3f \$ \\\\ \n",tol,nsteps_a,nsteps,time_a,time)
        end
        if scheme==schemes[3]
            @printf(    "\\texttt{PP 5/6 A}, TOL = \$ 10^{%d} \$ & \$ %d \$ & \$ %d \$ & \$ %5.3f \$ & \$ %5.3f \$ \\\\ \\hline \n",tol,nsteps_a,nsteps,time_a,time)
        end
        if scheme==schemes[4]
            @printf("\\texttt{Emb 4/3 AK p (defect)}, TOL = \$ 10^{%d} \$ & \$ %d \$ & \$ %d \$ & \$ %5.3f \$ & \$ %5.3f \$ \\\\ \n",tol,nsteps_a,nsteps,time_a,time)
        end
        if scheme==schemes[5]
            @printf(    "\\texttt{PP 3/4 A (defect)}, TOL = \$ 10^{%d} \$ & \$ %d \$ & \$ %d \$ & \$ %5.3f \$ & \$ %5.3f \$ \\\\ \n",tol,nsteps_a,nsteps,time_a,time)
        end
        if scheme==schemes[6]
            @printf(    "\\texttt{PP 5/6 A (defect)}, TOL = \$ 10^{%d} \$ & \$ %d \$ & \$ %d \$ & \$ %5.3f \$ & \$ %5.3f \$ \\\\ \\hline \n",tol,nsteps_a,nsteps,time_a,time)
        end
end