include("../RCP3BP.jl")
include("constrains.jl")
include("definitions.jl")
using JLD2
using DataStructures

#Crear nova seccio
#   Identificador (files system)
#   Definir condicions

#Cambiar seccio

#Generar orbites per una mu determinada

#Guardar les orbites

#Carregar orbites

#Generar orbites intermitges
#   2 punts consecutius molt allunyats
#   aproximar la colisio

#Current state
sectionLoaded = false
dataLoaded = false
dataSaved = true

currSection = ""

#Current data
mu = 0.0
eqCond = eq_dq1
dirCond = gt_dq2
Ln = 4
center = Vector{Float64}(undef,4)
vecs = Array{Vector{Float64}}(undef,2)
backwards = false
eps::Float64 = 0.0

data = SortedDict{Float64,Vector{intersection}}()

#Seccions

function newSection(name::String, def_eqCond, def_dirCond, Lagrange=4, isBackwards=false; override = false)
    directory = "./data/"*name
    if(isdir(directory))
	if override
	    println("Previus data deleted\n")
	    rm(directory,recursive=true,force=true)
	else
	    error("Seccio "*name*" ja existeix")
	end
    end
    mkdir(directory)

    global currSection = name
    global eqCond = def_eqCond
    global dirCond = def_dirCond

    global Ln = Lagrange
    global backwards = isBackwards

    global sectioLoaded = true

    save(directory*"/definition.jld2",
	 "eq", eqCond,
	 "dir", dirCond,
	 "Ln", Ln,
	 "backwards",backwards)
end

function loadSection(name::String)
    directory = "./data/"*name
    if(!isdir(directory))
	error("Seccio "*name*" no existeix")
    end

    definitions = load(directory*"/definition.jld2")
    global currSection = name
    global eqCond = definitions["eq"]
    global dirCond = definitions["dir"]
    global Ln = definitions["Ln"]
    global backwards = definitions["backwards"]

    global sectionLoaded = true
end

muDirectory(mu_val) = "./data/"*currSection*"/"*string(mu_val)*".jld2"

function setMu(new_mu::Float64; discard=false, override=false)
    if dataLoaded && !discard && !dataSaved
	error("Previous data not saved")
    end


    dir = muDirectory(new_mu)
    if isfile(dir) && !override
	println("Loaded existing data\n")
	loadMu(new_mu,override = discard)
    else
	global data = SortedDict{Float64,Vector{intersection}}()
	global dataSaved = true
    end

    global mu = new_mu
    #calcular el resto de valores
    if Ln == 4
	s = Sistem(mu)
	global center = s.L4
	global vecs = backwards ? s.Evecs : s.Ivecs
	global eps = s.eps
    elseif Ln == 5
	s = Sistem(mu)
	global center = s.L4 .* [1,-1,-1,1]
	global vecs = backwards ? s.Ivecs : s.Evecs
	global vecs[1] = vecs[1] .* [-1,1,1,-1]
	global vecs[2] = vecs[2] .* [1,-1,-1,1]
	global eps = s.eps
    end

    global dataLoaded = true
end

function newOrbit(theta::Float64,tspan, min_dist=0.1)
    u0 = center + eps*(cos(theta)*vecs[1] + sin(theta)*vecs[2])
    u0sa = SVector{4,Float64}(u0)
    cb = seccio_callback(eqCond,dirCond)
    
    global dataSaved = false

    prob = ODEProblem(backwards ? orbitBack : orbit,u0sa,tspan,mu)

    integrator = init(prob,Feagin14(),
		  dtmax = 0.01, 
		  callback = cb,
		  save_everystep = false,
		  save_start = false,
		  save_end = false)

    H_ini = hamiltonian((u0[3],u0[4]),(u0[1],u0[2]),mu)
    
    for(u,i) in tuples(integrator)
	H = hamiltonian((u[3],u[4]),(u[1],u[2]),mu)
	if 1e-10 < abs(H-H_ini)
	    break
	end
    end

    sol = integrator.sol
    new_col = Vector{intersection}(undef,0)
    for (t,u) in zip(sol.t,sol.u)
	if norm(center - u) < min_dist
	    continue
	end
	push!(new_col, intersection(u,t))
    end
    #afegir columna nova
    push!(data,theta=>new_col)
    
end

function newOrbitRange(angles::Array{Float64},tspan, min_dist=0.1)
    cb = seccio_callback(eqCond,dirCond)
    global dataSaved = false

    for theta in angles
	u0 = center + eps*(cos(theta)*vecs[1] + sin(theta)*vecs[2])
	u0sa = SVector{4,Float64}(u0)

	prob = ODEProblem(backwards ? orbitBack : orbit,u0sa,tspan,mu)

	integrator = init(prob,Feagin14(),
			dtmax = 0.01, 
			callback = cb,
			save_everystep = false,
			save_start = false,
			save_end = false)

	H_ini = hamiltonian((u0[3],u0[4]),(u0[1],u0[2]),mu)
    
	for(u,i) in tuples(integrator)
	    H = hamiltonian((u[3],u[4]),(u[1],u[2]),mu)
	    if 1e-10 < abs(H-H_ini)
		break
	    end
	end

	sol = integrator.sol
	new_col = Vector{intersection}(undef,0)
	for (t,u) in zip(sol.t,sol.u)
	    if norm(center - u) < min_dist #eliminar orbites mol properes al inici
		continue
	    end
	    push!(new_col, intersection(u,t))
	end
	#afegir columna nova
	push!(data,theta=>new_col)
    end
end

function findHoles(max_size = 1e-1)
    if !dataLoaded
	error("No loaded data")
    end
    global dataSaved = false

    max_length = 0
    for (k,V) in data
	max_length = max(max_length, length(V))
    end

    holes = Array{Float64}(undef,0)

    for step=1:max_length
	(prev_k,prev_V) =  last(data)
	if length(prev_V) >= step
	    prev_u = prev_V[step].u
	end

	for (k,V) in data
	    u = nothing
	    if length(V) >= step
		u = V[step].u
	    end

	    #check distance
	    if u != nothing && prev_u != nothing
		dist = norm(u-prev_u)
		if dist > max_size
		    push!(holes,(k+prev_k)/2)
		end
	    end

	    prev_u = u
	    prev_k = k
	end
   end
    return unique(holes)
end

function aproxCollisions(iters, tspan)
    holes = Array{Float64}(undef,0)
    for i=1:iters
	(prev_k,prev_V) = last(data)
	for (k,V) in data
	    if length(V) != length(prev_V)
		push!(holes,(k+prev_k)/2)
	    end
	    prev_k = k
	    prev_V = V
	end
	#calc new
	newOrbitRange(holes,tspan)
	empty!(holes)
    end
end

function saveMu()
    dir = muDirectory(mu)
    save_object(dir,data)
    global dataSaved = true
end

function loadMu(new_mu;override = false)
    if !dataSaved && !override
	error("Data not saved\n")
    end
    dir = muDirectory(new_mu)
    global data = load_object(dir)
    global dataSaved=true
end
