include("systems.jl")
include("LeviCivita.jl")

@enum orbit_type E I

H_0 = 0.1

function set_orbit_system(mu)
    L4 = [0.5-mu,sqrt(3)/2,-sqrt(3)/2,0.5-mu]
    global H_0 = hamiltonian(L4,mu)
end

function orbitWithColl(U,p,t)::SVector{5,Float64}
    mu,E = p
    u = SVector{4}([U[1],U[2],U[3],U[4]])
    coll = U[5]
    if coll == 2.0
	res = LCorbitP2(u,(mu,E),t) 
	return [res ;0.0]
    elseif coll == 1.0
	res = LCorbitP1(u,(mu,E),t)
	return [res ;0.0]
    else
	res = orbit(u,mu,t)
	return [res ; 0.0]
    end
    SA[0.0,0.0,0.0,0.0,0.0]
end

function orbitWithCollBack(U,p,t)::SVector{5,Float64}
    mu,E = p
    u = SVector{4}([U[1],U[2],U[3],U[4]])
    coll = U[5]
    if coll == 2.0
	res = LCorbitBackP2(u,(mu,E),t) 
	return [res ;0.0]
    elseif coll == 1.0
	res = LCorbitBackP1(u,(mu,E),t)
	return [res ;0.0]
    else
	res = orbitBack(u,mu,t)
	return [res ; 0.0]
    end
    SA[0.0,0.0,0.0,0.0,0.0]
end

function getCond()
    function condition(out,U,t,integrator)
	c = U[5]
	mu = integrator.p[1]
	u = SVector{4}([U[1],U[2],U[3],U[4]])
	if c == 0.0
	    q1 = u[1]
	    q2 = u[2]
	    out[1] = (q1+mu)^2 + q2^2 - 0.01
	    out[2] = (q1+mu-1)^2 + q2^2 - 0.01
	elseif c == 1.0
	    u1 = u[1]
	    u2 = u[2]
	    out[1] = (u1^2 + u2^2)^2 - 0.01
	elseif c == 2.0
	    u1 = u[1]
	    u2 = u[2]
	    out[2] = (u1^2 + u2^2)^2 - 0.01
	end
    end
end

function getAffect(dir::orbit_type)
    lt = dir == I ? (<) : (>)
    gt = dir == I ? (>) : (<)
    function affect!(integrator,idx)
	U = integrator.u
	c = U[5]
	mu = integrator.p[1]
	t = integrator.t
	u = SVector{4}([U[1],U[2],U[3],U[4]])
	if idx == 1
	    if c == 0
		q1 = u[1] + mu
		q2 = u[2]
		p1 = u[3]
		p2 = u[4] + mu
		if lt(q1*(p1+q2) + q2*(p2-q1), 0)
		   integrator.u = SVector{5}([qp2uU(SA[q1,q2,p1,p2]);1])
		end
	    elseif c == 1
		u1,u2,U1,U2 = u
		if gt((U1*u1 + U2*u2)/2 , 0)
		    u_new = uU2qp(u) + [-mu,0.0,0.0,-mu]
		    integrator.u = SVector{5}([u_new ; 0])
		end
	    end
	elseif idx == 2
	    if c == 0
		q1 = u[1] - 1 + mu
		q2 = u[2]
		p1 = u[3]
		p2 = u[4] - 1 + mu
		if lt(q1*(p1+q2) + q2*(p2-q1), 0)
		    integrator.u = SVector{5}([qp2uU(SA[q1,q2,p1,p2]);2])
		end
	    elseif c == 2
		u1,u2,U1,U2 = u
		if gt((U1*u1 + U2*u2)/2 , 0)
		    u_new = uU2qp(u) + [1-mu,0,0,1-mu]
		    integrator.u = SVector{5}([u_new ; 0])
		end
	    end
	end
    end
end

cbVect(dir) = VectorContinuousCallback(getCond(),getAffect(dir),2,save_positions=(false,false))

function sol2pos(U,mu)
    u = SVector{4}([U[1],U[2],U[3],U[4]])
    c = U[5]
    if c == 1
	return uU2qp(u) - [mu,0,0,mu]
    elseif c== 2
	return uU2qp(u) + [1-mu,0,0,1-mu]
    end
    return u
end

function calc_orbit(u0,mu,time)
    #TODO: check boundary
    q1,q2,p1,p2 = u0
    if (q1+mu)^2 + q2^2 < 0.01
	q1 = u0[1] + mu
	p2 = u0[4] + mu
	u0sa = SVector{5}([qp2uU(SA[q1,q2,p1,p2]);1])
    elseif (q1+mu-1)^2 + q2^2 < 0.01
	q1 = u0[1] - 1 + mu
	p2 = u0[4] - 1 + mu
	u0sa = SVector{5}([qp2uU(SA[q1,q2,p1,p2]);2])
    else
	u0sa = SVector{5,Float64}([u0;0])
    end

    
    cbE = cbVect(E)
    cbI = cbVect(I)

    probE = ODEProblem(orbitWithCollBack,
		      u0sa,(0.0,time),[mu,H_0])
    solE = solve(probE,Vern9(),
	    	abstol = 1e-14,reltol = 1e-14,
		dense = false,
		callback = cbE,
	    	#isoutofdomain = unstable_f(H0,1e-10)
	     )
    probI = ODEProblem(orbitWithColl,
		      u0sa,(0.0,time),[mu,H_0])
    solI = solve(probI,Vern9(),
	    	abstol = 1e-14,reltol = 1e-14,
		dense = false,
		callback = cbI,
	    	#isoutofdomain = unstable_f(H0,1e-10)
	     )
    pathE = sol2pos.(solE.u,mu)
    pathI = sol2pos.(solI.u,mu)
    return ([reverse(pathE[2:end]);pathI], length(pathE))
end

