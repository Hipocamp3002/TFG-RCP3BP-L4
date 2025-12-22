using JLD2
using CSV

systems = include("systems.jl")
sections = include("sections.jl")
include("LeviCivita.jl")

sections_keys = []
for k in keys(sections)
    push!(sections_keys,k)
end

function list()
    for (i,k) in pairs(sections_keys)
	println(string(i)*": "*k)
    end
end

@enum orbit_type E I

struct calc_parameters
    tspan
    max_intersections
end

@enum Missing no_reach h_error
struct Intersection
    u::Vector{Float64}
    t::Float64
    theta::Float64
end

function new_line(num_orbits)
    return Array{Union{Nothing,Intersection}}(nothing,num_orbits)
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

function getCond(S)
    eq_cond = S.eq_cond
    cut = S.cut
    function condition(out,U,t,integrator)
	c = U[5]
	mu = integrator.p[1]
	u = SVector{4}([U[1],U[2],U[3],U[4]])
	if c == 0.0
	    q1 = u[1]
	    q2 = u[2]
	    out[1] = (q1+mu)^2 + q2^2 - 0.01
	    out[2] = (q1+mu-1)^2 + q2^2 - 0.01
	    out[3] = eq_cond(u,t,mu) - cut
	elseif c == 1.0
	    u1 = u[1]
	    u2 = u[2]
	    out[1] = (u1^2 + u2^2)^2 - 0.01

	    u = uU2qp(u) - [mu,0,0,mu]
	    out[3] = eq_cond(u,t,mu) - cut
	elseif c == 2.0
	    u1 = u[1]
	    u2 = u[2]
	    out[2] = (u1^2 + u2^2)^2 - 0.01
	    
	    un = uU2qp(u) + [1-mu,0,0,1-mu]
	    out[3] = eq_cond(un,t,mu) - cut
	end
    end
end

function getAffect(S,dir::orbit_type,max_inter::Int64)
    lt = dir == I ? (<) : (>)
    gt = dir == I ? (>) : (<)
    dir_cond = S.dir_cond
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
	elseif idx == 3
	    un = sol2pos(U,mu)
	    if dir_cond(un,t,mu)
		savevalues!(integrator, true)
	    end
	    u_modified!(integrator,false)
	    if max_inter != -1 && length(integrator.sol) >= max_inter
		terminate!(integrator)
	    end
	end
    end
end

cbVect(S,dir,max_inter) = VectorContinuousCallback(getCond(S),getAffect(S,dir,max_inter),3,
				 save_positions = (false,false))

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


function uniformOrbits(system_key::String,section::Union{Int,String}, num_orbits::Int, mu::Float64, 
		       type::orbit_type = I, calc_p = calc_parameters((0.0,100),-1);
		       override = false)
    #get key if integer
    skey::String = ""
    if (section isa Int)
	skey = sections_keys[section]
    else
	skey = section
    end
    S = sections[skey]

    #mu_exists
    mu_dir = "data/"*string(mu)
    if !isdir(mu_dir)
	mkdir(mu_dir)
    end
    
    #system_exists
    sys_dir = mu_dir*"/"*system_key
    if !isdir(sys_dir)
	mkdir(sys_dir)
	sis = systems[system_key](mu)
	save_object(sys_dir*"/params.jld2",sis)
    else
	sis = load_object(sys_dir*"/params.jld2")
    end

    #section_exists
    sec_dir = sys_dir*"/"*skey
    if !isdir(sec_dir)
	mkdir(sec_dir)
    end

    #type_exists
    type_dir = sec_dir*"/"*(type == I ? "I" : "E")
    if !isdir(type_dir)
	mkdir(type_dir)
    elseif override
	rm(type_dir,recursive=true)
	mkdir(type_dir)
    else
	error("Orbit already exists, set override to true to override")
    end

    #calculate
    data = []
    #push!(data,new_line(num_orbits))
    push!(data,[])

    range = collect(LinRange(0.0,2pi,num_orbits))

    cb = cbVect(S,type,calc_p.max_intersections)
    vecs = type == I ? sis.Ivecs : sis.Evecs
    
    H0 = hamiltonian(sis.center,mu)

    for (i,theta) in enumerate(range)
	u0 = sis.center + sis.eps*(vecs[1]*cos(theta) + vecs[2]*sin(theta))
	u0sa = SVector{5,Float64}([u0;0])
	prob = ODEProblem(type == I ? orbitWithColl : orbitWithCollBack,
		   u0sa,calc_p.tspan,[mu,H0])
	sol = solve(prob,Vern9(),
	    	abstol = 1e-14,reltol = 1e-14,
		callback = cb,
		save_everystep = false,
		save_start = false,
		save_end = false
	    	#isoutofdomain = unstable_f(H0,1e-10)
	     )
	path = sol2pos.(sol.u,mu)
	

	if length(sol.u) > length(data)
	    for test in length(data):length(sol)
		push!(data,[])
	    end
	end

	i_inter = 1
	for (u,t) in zip(path,sol.t)
	    push!(data[i_inter], Intersection(u,t,theta))
	    i_inter += 1
	end
    end

    #save data
    for (i,row) in enumerate(data)
	path = type_dir*"/"*string(i)*".csv"
	CSV.write(path,row)
    end
end

function unstable_f(H_ini,mu,err)
    return (dt,u,p,t) -> err < abs(H_ini - hamiltonian(u,mu))
end
