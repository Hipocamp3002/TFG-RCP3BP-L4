using StaticArrays

function qp2uU(u)
    q1,q2,p1,p2 = u
    
    if q1 == 0 && q2 == 0
	throw(DomainError(u,"First terms (q1 q2) can not be both 0"))
    end

    if q2 != 0
	u = sqrt(q1^2 + q2^2)
	u1sq = 0.5*(q1 + u)

	u1 = sqrt(u1sq)
	u2 = q2/(u1*2)
    else
	if q1 > 0
	    u2 = 0
	    u1 = sqrt(q1)
	elseif q1 < 0
	    u1 = 0
	    u2 = sqrt(-q1)
	elseif q1 == 0
	    u1 = 0
	    u2 = 0
	end
    end

    U1 = 2*(u1*p1 + u2*p2)
    U2 = 2*(u1*p2 - u2*p1)
    
    SA[u1,u2,U1,U2]
end

function uU2qp(u)
    u1,u2,U1,U2 = u
    q1 = u1^2 - u2^2
    q2 = 2*u1*u2
    usq = u1^2 + u2^2
    p1 = 0.5*(U1*u1 - U2*u2)/usq
    p2 = 0.5*(U1*u2 + U2*u1)/usq

    SA[q1,q2,p1,p2]
end

function LCorbitP1(u,p,t)
    mu,E = p
    nu = 1-mu
    u1,u2,U1,U2 = u
    
    usq = u1^2 + u2^2

    du1 = 0.25*(U1 + 2*usq*u2)
    du2 = 0.25*(U2 - 2*usq*u1)
    
    a = (1 - 2*(u1^2-u2^2)+usq^2)
    A = 1/sqrt(a) - u1^2 + u2^2
    a32 = sqrt(a)*a

    Adu1 = -0.5*(4*usq*u1 - 4*u1)/a32 - 2*u1
    Adu2 = -0.5*(4*usq*u2 + 4*u2)/a32 + 2*u2

    E1 = E + (mu^2)/2
    #Negats al retornar
    dU1 = 4*du1*u1*u2 - 2*du2*(usq+2*u1^2) - 3*u1*usq^2 - 2*u1*E1 - mu*(2*u1*A + usq*Adu1)
    dU2 = 2*du1*(usq+2*u2^2) - 4*du2*u1*u2 - 3*u2*usq^2 - 2*u2*E1 - mu*(2*u2*A + usq*Adu2)

    SA[du1,du2,-dU1,-dU2]
end

function LCorbitBackP1(u,p,t)
    mu,E = p
    nu = 1-mu
    u1,u2,U1,U2 = u
    
    usq = u1^2 + u2^2

    du1 = 0.25*(U1 + 2*usq*u2)
    du2 = 0.25*(U2 - 2*usq*u1)
    
    a = (1 - 2*(u1^2-u2^2)+usq^2)
    A = 1/sqrt(a) - u1^2 + u2^2
    a32 = sqrt(a)*a

    Adu1 = -0.5*(4*usq*u1 - 4*u1)/a32 - 2*u1
    Adu2 = -0.5*(4*usq*u2 + 4*u2)/a32 + 2*u2

    E1 = E + (mu^2)/2
    #Negats al retornar
    dU1 = 4*du1*u1*u2 - 2*du2*(usq+2*u1^2) - 3*u1*usq^2 - 2*u1*E1 - mu*(2*u1*A + usq*Adu1)
    dU2 = 2*du1*(usq+2*u2^2) - 4*du2*u1*u2 - 3*u2*usq^2 - 2*u2*E1 - mu*(2*u2*A + usq*Adu2)

    SA[-du1,-du2,dU1,dU2]
end

function LCorbitP2(u,p,t)
    mu,E = p
    nu = 1-mu
    u1,u2,U1,U2 = u
    
    usq = u1^2 + u2^2

    du1 = 0.25*(U1 + 2*usq*u2)
    du2 = 0.25*(U2 - 2*usq*u1)
    
    a = (1 + 2*(u1^2-u2^2)+usq^2)
    A = 1/sqrt(a) + u1^2 - u2^2
    a32 = sqrt(a)*a

    Adu1 = -0.5*(4*u1 + 4*usq*u1)/a32 + 2*u1
    Adu2 = -0.5*(4*usq*u2 - 4*u2)/a32 - 2*u2

    E1 = E + (nu^2)/2
    #Negats al retornar
    dU1 = 4*du1*u1*u2 - 2*du2*(usq+2*u1^2) - 3*u1*usq^2 - 2*u1*E1 - nu*(2*u1*A + usq*Adu1)
    dU2 = 2*du1*(usq+2*u2^2) - 4*du2*u1*u2 - 3*u2*usq^2 - 2*u2*E1 - nu*(2*u2*A + usq*Adu2)

    SA[du1,du2,-dU1,-dU2]
end

function LCorbitBackP2(u,p,t)
    mu,E = p
    nu = 1-mu
    u1,u2,U1,U2 = u
    
    usq = u1^2 + u2^2

    du1 = 0.25*(U1 + 2*usq*u2)
    du2 = 0.25*(U2 - 2*usq*u1)
    
    a = (1 + 2*(u1^2-u2^2)+usq^2)
    A = a^-1 + u1^2 - u2^2
    a32 = sqrt(a)*a

    Adu1 = (4*u1 + 4*usq*u1)/a32 + 2*u1
    Adu2 = (4*usq*u2 - 4*u2)/a32 - 2*u2

    E1 = E+(nu^2)/2
    #Negats al retornar
    dU1 = 4*du1*u1*u2 - 2*du2*(usq+2*u1^2) - 3*u1*usq^2 - 2*u1*E1 - nu*(2*u1*A + usq*dAu1)
    dU2 = 2*du1*(usq+2*u2^2) - 4*du2*u1*u2 - 3*u2*usq^2 - 2*u1*E1 - nu*(2*u2*A + usq*dAu2)

    SA[-du1,-du2,dU1,dU2]
end

function hamiltonianP2(u,mu,E)
    u1,u2,U1,U2 = u
    nu = 1 - mu
    usq = u1^2 + u2^2
    return (1/8)*((U1+2*usq*u2)^2 + (U2 - 2*usq*u1)^2) - 0.5*usq^3 - mu - usq*(E+(nu^2)*0.5) - nu*usq*(1/sqrt(1+2*(u1^2 - u2^2) + usq^2) + u1^2 - u2^2)
end
