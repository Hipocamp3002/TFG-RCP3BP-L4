using GLMakie
include("generator.jl")



µ_1 = 0.5*(1-sqrt(69)/9)

#LAYOUT
f = Figure(size = (600,400))
gl = Axis(f[1,1], autolimitaspect = 1)
gr = f[1,2] = GridLayout()


µ_grid = gr[1,1] = GridLayout()
Label(µ_grid[1,1],"µ=")
µ_input = Textbox(µ_grid[1,2], placeholder = "0.5", validator=Float64,tellwidth = false)
µ_slider = Slider(µ_grid[2,1:2], range=LinRange(µ_1,0.5,1000), startvalue = 0.5,update_while_dragging = false)

ang_grid = gr[2,1] = GridLayout()
Label(ang_grid[1,1:2],"Rang de angle inicial d'Estable")
angE_min = Textbox(ang_grid[2,1],placeholder= "0.0", validator=Float64)
angE_max = Textbox(ang_grid[2,2],placeholder= "360.0", validator=Float64)
Label(ang_grid[3,1:2],"Rang de angle inicial d'Inestable")
angI_min = Textbox(ang_grid[4,1],placeholder= "0.0", validator=Float64)
angI_max = Textbox(ang_grid[4,2],placeholder= "360.0", validator=Float64)

E_grid = gr[3,1] = GridLayout()
Label(E_grid[1,1],"Orbita estable")
menu_E = Menu(E_grid[1,2], options = ["L4","L5"],default = "L4")

I_grid = gr[4,1] = GridLayout()
Label(I_grid[1,1],"Orbita inestable")
menu_I = Menu(I_grid[1,2], options = ["L4","L5"],default = "L4")

sec_grid = gr[5,1] = GridLayout()
Label(sec_grid[1,1:2],"Seccio a calcular")
Label(sec_grid[2,1],"Tipus:")
menu_type = Menu(sec_grid[2,2], options = ["L4","L5","C1","C2","0"],default = "L4")
Label(sec_grid[3,1],"Radi:")
input_radi = Textbox(sec_grid[3,2], placeholder = "0.1", validator=Float64)
Label(sec_grid[4,1],"Direccio:")
menu_dir = Menu(sec_grid[4,2], options = ["+","-"],default = "+")

p_grid = gr[6,1] = GridLayout()
Label(p_grid[1,1],"num. punts")
input_num_pE = Textbox(p_grid[1,2],placeholder = "2000", validator=Int64)
Label(p_grid[2,1],"num. punts")
input_num_pI = Textbox(p_grid[2,2],placeholder = "2000", validator=Int64)

button_calc = Button(gr[7,1],label= "Recalcular")

dataE = Float64[]
dataI = Float64[]
pointsE_plt = Observable(Point2f[])
pointsI_plt = Observable(Point2f[])

scatter!(gl,pointsE_plt,color="blue",markersize=5,
	 inspector_label = (_,i,_) -> begin
	 "θ = "*string(dataE[i])
	 end)
scatter!(gl,pointsI_plt,color="red",markersize=5,
	inspector_label = (_,i,_) -> begin
	 "θ = "*string(dataI[i])
	 end)

Box(gr[end+1,1],visible=false)
colsize!(gr,1,Fixed(200))

#GLOBALS
µ = 0.5

num_pointsE = 2000
num_pointsI = 2000

minE_theta = 0
maxE_theta = 2pi
minI_theta = 0
maxI_theta = 2pi

radius = 0.1

#INTERACCIO

c = Channel()

on(µ_input.stored_string) do s
    val = parse(Float64,s)
    if val >= µ_1 && val <= 0.5
	µ = val
    else
	µ = 0.5
	println("Invalid µ")
    end

    #calc_inter()
    put!(c,1)
end

on(µ_slider.value) do v
    µ_input.stored_string[] = string(v)
    µ_input.displayed_string = string(v)
end

on(input_num_pE.stored_string) do s
    val = parse(Int64,s)
    if val <= 0
	global num_pointsE = 1
    else
	global num_pointsE = val
    end
    
end
on(input_num_pI.stored_string) do s
    val = parse(Int64,s)
    if val <= 0
	global num_pointsI = 1
    else
	global num_pointsI = val
    end
    
end

on(angE_min.stored_string) do s
    global minE_theta = deg2rad(parse(Float64,s))
    while minE_theta > maxE_theta
	global maxE_theta += 2pi
    end
end
on(angE_max.stored_string) do s
    global maxE_theta = deg2rad(parse(Float64,s))
    while minE_theta > maxE_theta
	global maxE_theta += 2pi
    end 
end
on(angI_min.stored_string) do s
    global minI_theta = deg2rad(parse(Float64,s))
    while minI_theta > maxI_theta
	global maxI_theta += 2pi
    end
end
on(angI_max.stored_string) do s
    global maxI_theta = deg2rad(parse(Float64,s))
    while minI_theta > maxI_theta
	global maxI_theta += 2pi
    end 
end
on(input_radi.stored_string) do s
    val = parse(Float64,s)
    global radius = max(0.0,val)
end

on(button_calc.clicks) do n
    #calc_inter()
    put!(c,1)
end

#ACCIONS

function eq_sec(center::SVector{4,Float64})
    function eq_Pos(u,t,mu)
	U = u - center
	U[1]*U[1] + U[2]*U[2]
    end
end
function dir_sec(center::SVector{4,Float64},dir)
    comp = (dir == "+") ? (>) : (<)
    function dir_Pos(u,t,mu)
	U = u - center
	comp(2*U[1]*(u[2]+u[3]) + 2*U[2]*(u[4]-u[1]),0)
    end
end

function num2ang(n)
    a = 2^ceil(log2(n))
    return (2*n-a-1)/a
end


function calc_inter()
    empty!(pointsE_plt.val)
    empty!(pointsI_plt.val)
    empty!(dataE)
    empty!(dataI)

    sisE = (menu_E.selection[] == "L5") ? sistema_L5(µ) : sistema_L4(µ)
    sisI = (menu_I.selection[] == "L5") ? sistema_L5(µ) : sistema_L4(µ)

    if menu_type.selection[] == "L4"
	center = SA[0.5-µ,sqrt(3)/2,-sqrt(3)/2,0.5-µ]
    elseif menu_type.selection[] == "L5"
	center = SA[0.5-µ,-sqrt(3)/2,sqrt(3)/2,0.5-µ]
    elseif menu_type.selection[] == "B1"
	center = SA[-µ,0,0,-µ]
    elseif menu_type.selection[] == "B2"
	center = SA[1-µ,0,0,1-µ]
    else
	center = SA[0.0,0.0,0.0,0.0]
    end

    S = section(
	radius^2,
	eq_sec(center),
	dir_sec(center,menu_dir.selection[])
    )
    cbE = cbVect(S,E,1)
    cbI = cbVect(S,I,1)
    
    H0 = hamiltonian(sisE.center,µ)#L4 i L5 son simetrics

    @async for i = 1:num_pointsE
	theta = (maxE_theta-minE_theta)*num2ang(i)+minE_theta
	u0 = sisE.center + sisE.eps*(sisE.Evecs[1]*cos(theta)+sisE.Evecs[2]*sin(theta))
	u0sa = SVector{5,Float64}([u0;0])
	    prob = ODEProblem(orbitWithCollBack,u0sa,(0,50),[µ,H0])
	    sol = solve(prob,Vern9(),
		   abstol = 1e-14,reltol = 1e-14,
		    callback = cbE,
		    save_everystep = false,
		    save_start = false,
		    save_end = false
	    	#isoutofdomain = unstable_f(H0,1e-10)
		)

	    if length(sol.u) > 0
		p = sol2pos(sol.u[1],µ)
		push!(pointsE_plt[],Point2f(p[3],p[4]))
		push!(dataE,rad2deg(theta))
		notify(pointsE_plt)
	    end
	yield()
    end
    @async for i = 1:num_pointsI
	theta = (maxI_theta-minI_theta)*num2ang(i)+minI_theta
	    u0 = sisI.center + sisI.eps*(sisI.Ivecs[1]*cos(theta)+sisI.Ivecs[2]*sin(theta))
	    u0sa = SVector{5,Float64}([u0;0])
	    prob = ODEProblem(orbitWithColl,u0sa,(0,50),[µ,H0])
	    sol = solve(prob,Vern9(),
		   abstol = 1e-14,reltol = 1e-14,
		    callback = cbI,
		    save_everystep = false,
		    save_start = false,
		    save_end = false
	    	#isoutofdomain = unstable_f(H0,1e-10)
		)

	    if length(sol.u) > 0
		p = sol2pos(sol.u[1],µ)
		push!(pointsI_plt[],Point2f(p[3],p[4]))
		push!(dataI,rad2deg(theta))
		notify(pointsI_plt)
	    end
	yield()
    end
    
end


function gen_orbits()
    first = true
    while (true)
	println("hey")
	if isready(c) || first
	    if !first
		take!(c)
	    end
	    first = false

	    empty!(pointsE_plt.val)
	    empty!(pointsI_plt.val)
	    empty!(dataE)
	    empty!(dataI)

	    sisE = (menu_E.selection[] == "L5") ? sistema_L5(µ) : sistema_L4(µ)
	    sisI = (menu_I.selection[] == "L5") ? sistema_L5(µ) : sistema_L4(µ)

	    if menu_type.selection[] == "L4"
		center = SA[0.5-µ,sqrt(3)/2,-sqrt(3)/2,0.5-µ]
	    elseif menu_type.selection[] == "L5"
		center = SA[0.5-µ,-sqrt(3)/2,sqrt(3)/2,0.5-µ]
	    elseif menu_type.selection[] == "B1"
		center = SA[-µ,0,0,-µ]
	    elseif menu_type.selection[] == "B2"
		center = SA[1-µ,0,0,1-µ]
	    else
		center = SA[0.0,0.0,0.0,0.0]
	    end

	    S = section(
		radius^2,
		eq_sec(center),
		dir_sec(center,menu_dir.selection[])
	    )
	    cbE = cbVect(S,E,1)
	    cbI = cbVect(S,I,1)
    
	    H0 = hamiltonian(sisE.center,µ)#L4 i L5 son simetrics
	end
	if length(pointsE_plt.val) < num_pointsE
	    theta = (maxE_theta-minE_theta)*rand()+minE_theta
	    u0 = sisE.center + sisE.eps*(sisE.Evecs[1]*cos(theta)+sisE.Evecs[2]*sin(theta))
	    u0sa = SVector{5,Float64}([u0;0])
	    prob = ODEProblem(orbitWithCollBack,u0sa,(0,50),[µ,H0])
	    sol = solve(prob,Vern9(),
		   abstol = 1e-14,reltol = 1e-14,
		    callback = cbE,
		    save_everystep = false,
		    save_start = false,
		    save_end = false
	    	#isoutofdomain = unstable_f(H0,1e-10)
		)

	    if length(sol.u) > 0
		p = sol2pos(sol.u[1],µ)
		push!(pointsE_plt[],Point2f(p[3],p[4]))
		push!(dataE,rad2deg(theta))
		notify(pointsE_plt)
	    end
	end
	if length(pointsI_plt.val) < num_pointsI
	    theta = (maxI_theta-minI_theta)*rand()+minI_theta
	    u0 = sisI.center + sisI.eps*(sisI.Ivecs[1]*cos(theta)+sisI.Ivecs[2]*sin(theta))
	    u0sa = SVector{5,Float64}([u0;0])
	    prob = ODEProblem(orbitWithColl,u0sa,(0,50),[µ,H0])
	    sol = solve(prob,Vern9(),
		   abstol = 1e-14,reltol = 1e-14,
		    callback = cbI,
		    save_everystep = false,
		    save_start = false,
		    save_end = false
	    	#isoutofdomain = unstable_f(H0,1e-10)
		)

	    if length(sol.u) > 0
		p = sol2pos(sol.u[1],µ)
		push!(pointsI_plt[],Point2f(p[3],p[4]))
		push!(dataI,rad2deg(theta))
		notify(pointsI_plt)
	    end
	end

    end
end

Threads.@spawn gen_orbits()

DataInspector(gl)
f
