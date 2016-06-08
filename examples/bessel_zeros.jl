function besselj_zero_approx1(ν::Integer, m::Integer)
    μ = 4.0*ν^2
    a = (m+ν/2-1/4)*π
    a8 = 1/(8a)
    a82 = a8^2
    j = (a -(μ-1)*a8*(1 + a82*(4*(7μ-31)/3 +a82*( 32*(83*μ^2-982*μ+3779)/15 +
         a82*(64*(6949*μ^3-153855*μ^2+1585743*μ-6277237)/105))))
        )
end

function besseljprime_zero_approx1(ν::Integer, m::Integer)
    #if (ν==0 && m==1)
    #    return 0.0
    #end
    if ν==0
        m += 1
    end
    μ = 4*ν^2
    b = (m+ν/2-3/4)*π
    b8 = 1/(8*b)
    b82 = b8^2
    j = (b - b8*(μ+3 + b82*(4*(7*μ^2+82*μ-9)/3 +b82*(32*(83*μ^3+2075*μ^2-3039*μ+3537)/15 + 
             b82*(64*(6949*μ^4+296492*μ^3-1248002*μ^2+7414380*μ-5853627)/105)))))
end

function airy_zero_approx(m::Integer)
    if m<=10
        z =[-2.338107410459767,
            -4.087949444130970,
            -5.520559828095552,
            -6.786708090071759,
            -7.944133587120853,
            -9.022650853340980,
           -10.04017434155809,
           -11.00852430373326,
           -11.93601556323626,
           -12.82877675286576][m]
    else    
        t = 3/8*π*(4*m-1)
        t2 = t^(-2)
        z = -t^(2/3)*(1+t2*(5/48+t2*(-5/36+t2*(77125/82944+t2*(
            -108056875/6967296+ t2*162375596875/334430208)))))
    end   
    z
end

function airyprime_zero_approx(m::Integer)
    if m<=10
        z = [-1.0187929716474710890,
             -3.2481975821798365379,
             -4.8200992111787356394,
             -6.1633073556394865476,
             -7.3721772550477701771,
             -8.4884867340197221329,
             -9.5354490524335474707,
            -10.527660396957407282,
            -11.475056633480245295,
            -12.384788371845747325][m]
    else    
        t = 3/8*π*(4*m-3)
        t2 = t^(-2)
        z = -t^(2/3)*(1+t2*(-7/48+t2*(35/288+t2*(-181223/207360+t2*(
            18683371/1244160 - t2*9114588436/191102976)))))
    end   
    z
end

function besselj_zero_approx2(ν::Integer, m::Integer; prime::Bool=false)
    if prime 
        ζ = ν^(-2/3)*airyprime_zero_approx(m)
    else
         ζ = ν^(-2/3)*airy_zero_approx(m)
    end
    y = 2/3*(-ζ)^(3/2)
    if y>100000
        x = π/2
    elseif y<1
        p = (3*y)^(1/3)
        p2 = p^2
        x = p*(1+p2*(-2/15+p2*(3/175+p2*(-2/1575))))
    else
        p = 1/(y+π/2)
        p2 = p^2
        x = π/2 - p*(1+p2*(2/3+p2*(13/15+p2*(146/105+p2*(781/315+p2*16328/3465)))))        
    end
    x2 = (y+x)^2
    r = (x-atan(x+y))/x2
    x = x - (1+x2)*r*(1+r/(x+y))

    z = 1/cos(x)
    h = sqrt(ζ*(1-z^2))
    if prime
        g1 = z/ζ*h*(7/(48*ζ) + h*(7/(z^2-1)+9
        )/24)
        j = ν*z + g1/ν
    else
        f1 = -z/ζ*h*(5/(48*ζ) + h*(5/(z^2-1)+3)/24)
        j = ν*z + f1/ν
    end
    j
end

function besselj_zero_approx(ν::Integer, m::Integer)
  if m>=ν
     j = besselj_zero_approx1(ν, m)
  else
     j = besselj_zero_approx2(ν, m)
  end      
  j
end

function besseljprime_zero_approx(ν::Integer, m::Integer)
  if m>=ν
        j = besseljprime_zero_approx1(ν, m)
  else
        j = besselj_zero_approx2(ν, m, prime=true)
  end      
  j
end

function besselj_zero_iter(ν::Integer, z::AbstractFloat)
    T = typeof(z)
    ɛ = eps(T)
    for i = 1:200
        J = besselj(ν, z)
        if abs(J) < 1000ɛ 
            break
        end        
        Jprime = besselj(ν-1, z) - ν*J/z
        z -= J/Jprime   
        if i==200
            println("200 iterations, res=",abs(J))
        end
    end
    return z
end

function besseljprime_zero_iter(ν::Integer, z::AbstractFloat)    
    T = typeof(z)
    ɛ = eps(T)
    for i = 1:200
        J = besselj(ν, z)
        J1 = besselj(ν-1, z) - ν*J/z
        if abs(J1) < 1000ɛ 
            break
        end       
        J2 = -J1/z + ((ν/z)^2-1)*J
        z -= J1/J2   
        if i==200
            println("200 iterations, res=",abs(J1))
        end
    end
    return z
end

besselj_zero(ν, m; T=Float64) = besselj_zero_iter(ν, convert(T, besselj_zero_approx(ν, m)))

besseljprime_zero(ν::Integer, m::Integer; T=Float64) =
    besseljprime_zero_iter(ν, convert(T, besseljprime_zero_approx(ν, m)))


