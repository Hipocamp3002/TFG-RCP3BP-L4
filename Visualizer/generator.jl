using JLD2
using CSV

sections = include("sections.jl")

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

function uniformOrbits(section::Union{Int,String}, num_orbits::Int, mu::Float64, 
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

    #check if section exists
    sec_dir = "data/"*skey
    if !isdir(sec_dir)
	mkdir(sec_dir)
    end

    #mu folder exists
    mu_dir = sec_dir*"/"*string(mu)
    if !isdir(mu_dir)
	mkdir(mu_dir)
	#calc sistem parameters
	sis = S.constructor(mu)
	save_object(mu_dir*"/params.jld2", sis)
    else
	sis = load_object(mu_dir*"/params.jld2")
    end

    #check if type orbit exists
    if type == I
	type_dir = mu_dir*"/I"
    else
	type_dir = mu_dir*"/E"
    end
    if !isdir(type_dir) mkdir(type_dir) else
	if override
	    rm(type_dir, recursive=true)
	    mkdir(type_dir)
	else
	    error("This sis is already calculated")
	end
    end

    #calculate
    data = []
    #push!(data,new_line(num_orbits))
    push!(data,[])

    range = collect(LinRange(0.0,2pi,num_orbits))

    cond(u,t,integrator) = S.eq_cond(u,mu)
    function affect!(integrator)
	u = integrator.u
	if S.dir_cond(u,mu)
	    res = savevalues!(integrator, true)
	end
	u_modified!(integrator,false)
    end
    cb = ContinuousCallback(cond,affect!,
			    save_positions = (false,false))
    vecs = type == I ? sis.Ivecs : sis.Evecs

    for (i,theta) in enumerate(range)
	u0 = sis.center + sis.eps*(vecs[1]*cos(theta) + vecs[2]*sin(theta))
	u0sa = SVector{4,Float64}(u0)
	prob = ODEProblem(type == I ? orbit : orbitBack,u0sa,calc_p.tspan,mu)

	integrator = init(prob,Vern9(),
			abstol = 1e-14, reltol = 1e-14, 
			callback = cb,
			save_everystep = false,
			save_start = false,
			save_end = false)
	
	H_ini = hamiltonian(u0,mu)
	for (u,t) in tuples(integrator)
	    H = hamiltonian(u,mu)
	    if 1e-10 < abs(H_ini-H)
		break#stop calculating
	    end
	    if calc_p.max_intersections != -1
		if calc_p.max_intersections <= length(integrator.sol) break end
	    end
	end

	sol = integrator.sol
	
	#sol = solve(prob,Feagin14(),
	#	dtmax=0.01,
		#abstol = 1e-14,reltol = 1e-14,
	#	callback = cb,
	#	save_everystep = false,
	#	save_start = false,
	#	save_end = false,
	#	unstable_check = unstable_f(hamiltonian(u0,mu),mu,1e-10))
#	
	if length(sol.u) > length(data)
	    for test in length(data):length(sol)
		#push!(data,new_line(num_orbits))
		push!(data,[])
	    end
	end

	i_inter = 1
	for (u,t) in zip(sol.u,sol.t)
	    #data[i_inter][i] = Intersection(u,t,theta)
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
