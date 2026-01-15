using GLMakie
include("systems.jl")
include("LeviCivita.jl")

@enum orbit_type E I

#Layout
fig  = Figure()
ing = fig[1,1] = GridLayout(tellwidth=false)
mu_input = Textbox(ing[1,1], validator=Float64, stored_string="0.5")
time_input = Textbox(ing[1,2],validator=Float64, stored_string="60")
num_input = Textbox(ing[1,3],validator=Int64, stored_string="100")

ax = Axis3(fig[2,1],perspectiveness = 0.8,
	   aspect = :data,
	   #viewmode = :strech
	   )

#data
mu = 0.5
time = 60
num_points = 100

pointsE = Observable(Point3f[])
pointsI = Observable(Point3f[])

#input

on(mu_input.stored_string) do v
    global mu = clamp(parse(Float64,v),0.04,0.5)
    recalc()
end
on(time_input.stored_string) do v
    global time = max(parse(Float64,v),0)
    recalc()
end
on(num_input.stored_string) do v
    global num_points = max(parse(Int64,v),0)
    recalc()
end

#calc

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


function recalc()
    s = sistema_L4(mu)
    for dir in [E,I]
	points = dir == E ? pointsE : pointsI
	empty!(points[])
	vecs = dir == E ? s.Evecs : s.Ivecs
	for angle in LinRange(0,2pi,num_points) |> collect
	    u0 = s.center + s.eps*(cos(angle)*vecs[1] + sin(angle)*vecs[2])
	    H0 = hamiltonian(s.center,mu)
	    u0sa = SVector{5,Float64}([u0;0])
    
	    cb = cbVect(dir)

	    prob = ODEProblem(dir==E ? orbitWithCollBack : orbitWithColl,
			    u0sa,(0.0,time),[mu,H0])
	    sol = solve(prob,Vern9(),
			abstol = 1e-14,reltol = 1e-14,
			dense = true,
			callback = cb)

	    path = sol2pos.(sol.u,mu)
	    #proj = [Point3f(u[1],u[2],atan(u[4]-u[1],u[3]+u[2])|>cos) for u in path]
	    #proj = [Point3f(u[1],u[2],u[4]-u[1]/u[3]+u[2]) for u in path]
	    proj = [begin
		ang = atan(u[4]-u[1],u[3]+u[2])
		Point3f(cos(ang)*u[1],sin(ang)*u[1],u[2])
		end for u in path]
	    append!(points[],proj)
	end
	notify(points)
    end
end

#display

recalc()
scatter!(ax,pointsE,color="blue",depthsorting=true,markersize=5)
scatter!(ax,pointsI,color="red",depthsorting=true,markersize=5)
display(fig) |> wait
