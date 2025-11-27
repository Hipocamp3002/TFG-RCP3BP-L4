using DataStructures
using GLMakie
using JLD2
using CSV

#Plot with selectors
#   - mu value
#   - 2 systems
#   - 1 section
#   - Axis selector



system1 = "none"
system2 = "none"

section1_path = "none"
section2_path = "none"

sec1 = 0
sec2 = 0
cur_section = 0

xaxis=1
yaxis=2

last_hover=0.0

curr_mu = 0.0


fig = Figure()

ax = Axis3(fig[2,1])
fig[1,1] = inputGrid_sec = GridLayout(tellwidth=false)
fig[3,1] = inputGrid_axis = GridLayout(tellwidth=false)

mu_selection = []
for mu_dir in readdir("data")
    push!(mu_selection,mu_dir)
end

mu_menu = Menu(fig,options=mu_selection,default=nothing)
sys1_menu = Menu(fig,options=["none"], default="none")
sys2_menu = Menu(fig,options=["none"], default="none")
sec_menu = Menu(fig,options=["none"], default="none")

axis_opt = ["q1","q2","p1","p2", "Q_L4_angle","P_L4_angle","Q_L4_mod","P_L4_mod"]
xaxis_menu = Menu(fig,options=zip(axis_opt,1:8)
		  , default = "q1")
yaxis_menu = Menu(fig,options=zip(axis_opt,1:8)
		  , default = "q2")
sec1_check = Toggle(fig,active = false)
sec1_input = Textbox(fig,placeholder="E",validator = Int64)
sec2_check = Toggle(fig,active = false)
sec2_input = Textbox(fig,placeholder="I",validator = Int64)

input_sec = inputGrid_sec[1,1:4] = [mu_menu,sys1_menu, sys2_menu, sec_menu]
input_axis = inputGrid_axis[1,1:6] = [xaxis_menu,yaxis_menu,
    sec1_check,sec1_input,sec2_check,sec2_input]


pointsE = Observable(Vector{Float64}[])
anglesE = Float64[]
pointsE_plt = lift(pointsE) do P
    points::Vector{Point3f} = [Point3f(project(u,xaxis),project(u,yaxis),project(u,8)) for u in P]
    return points
end
anglesI = Float64[]
pointsI = Observable(Vector{Float64}[])
pointsI_plt = lift(pointsI) do P
    points::Vector{Point3f} = [Point3f(project(u,xaxis),project(u,yaxis),project(u,8)) for u in P]
    return points
end

scatter!(ax,pointsE_plt,color="blue",markersize = 5,
	 inspector_label = (self,i,pos) -> begin
	    global last_hover = anglesE[i]
	    "θ = "*string(anglesE[i])
	 end)
scatter!(ax,pointsI_plt,color="red",markersize = 5,
	 inspector_label = (self,i,pos) -> begin
	    global last_hover = anglesI[i]
	    "θ = "*string(anglesI[i])
	 end)



#TODO: fix doble execution
function updateSectionsMenu()
    index = sec_menu.i_selected.val
    if index != 0
	id = sec_menu.options.val[index]
    else
	id = "none"
    end

    if system1 == "none" && system2 == "none"
	sec_menu.options = ["none"]
    elseif system1 == "none"
	sec2_opts = ["none"]
	for sec in readdir(system2)
	    sec_dir = system2*"/"*sec*"/I"
	    if isdir(sec_dir)
		push!(sec2_opts,sec)
	    end
	end
	sec_menu.options = sec2_opts
    elseif system2 == "none"
	sec1_opts = ["none"]
	for sec in readdir(system1)
	    sec_dir = system1*"/"*sec*"/E"
	    if isdir(sec_dir)
		push!(sec1_opts,sec)
	    end
	end
	sec_menu.options = sec1_opts
    else
	sec1_opts = ["none"]
	for sec in readdir(system1)
	    sec_dir = system1*"/"*sec*"/E"
	    if isdir(sec_dir)
		push!(sec1_opts,sec)
	    end
	end
	sec_menu.options = sec1_opts

	sec2_opts = ["none"]
	for sec in readdir(system2)
	    sec_dir = system2*"/"*sec*"/I"
	    if isdir(sec_dir)
		push!(sec2_opts,sec)
	    end
	end
	
	sec_menu.options = intersect(sec1_opts,sec2_opts)
    end

    for (n,i) in enumerate(sec_menu.options.val)
	if id == i
	    sec_menu.i_selected = n
	    notify(sec_menu.selection)
	    break
	end
    end
end

#drawing
function redraw()
    if section1_path != "none"
	global cur_section = sec1_check.active[] ? sec1 : 0
	update_points!(pointsE,anglesE,section1_path)
	notify(pointsE)
    else
	empty!(pointsE[])
	notify(pointsE)
    end

    if section2_path != "none"
	global cur_section = sec2_check.active[] ? sec2 : 0
	update_points!(pointsI,anglesI,section2_path)
	notify(pointsI)
    else
	empty!(pointsI[])
	notify(pointsI)
    end
end

function project(u,axis)
    if axis <= 4
	return u[axis]
    elseif axis == 5
	angle = atan(u[1]-0.5+curr_mu, u[2]-sqrt(3)/2)
	return angle
    elseif axis == 6
	angle = atan(u[3]+sqrt(3)/2, u[4]-0.5+curr_mu)
	return angle
    elseif axis == 7
	x = u[1]-0.5+curr_mu
	y = u[2]-sqrt(3)/2
	len = sqrt(x*x + y*y)
	return len
    elseif axis == 8
	x = u[3]-0.5+curr_mu
	y = u[4]+sqrt(3)/2
	len = sqrt(x*x + y*y)
	return len
    end
    return 0.0
end

function update_points!(points,angles,path)
    empty!(points.val)
    empty!(angles)
    if path == "none" return end
    for sec_file in readdir(path)
	section = parse(Int64,split(sec_file,".")[1])
	if cur_section != 0 && cur_section != section
	    continue
	end
	data = CSV.File(path*"/"*sec_file)
	for (uStr,t,theta) in data
	    u = parse.(Float64, split(chop(uStr,head=1,tail=1),','))
	    push!(points.val,u)
	    push!(angles,theta)
	end
    end
end


#events
on(mu_menu.selection) do s
    if s == nothing return
    else global curr_mu = parse(Float64,s)end
    mu_dir = "data/"*s
    systems_opt = [("none","none")]
    for sys in readdir(mu_dir)
	push!(systems_opt, (sys, mu_dir*"/"*sys))
    end
    
    #get current id of selection
    index = sys1_menu.i_selected.val
    if index != 0
	id1,_ = sys1_menu.options.val[index]
    else
	id1 = "none"
    end
    index = sys2_menu.i_selected.val
    if index != 0
	id2,_ = sys2_menu.options.val[index]
    else
	id2 = "none"
    end

	
    sys1_menu.options = systems_opt
    for (n,(i,d)) in enumerate(systems_opt)
	if id1 == i
	    sys1_menu.i_selected = n
	    break
	end
    end

    sys2_menu.options = systems_opt
    for (n,(i,d)) in enumerate(systems_opt)
	if id2 == i
	    sys2_menu.i_selected = n
	    break
	end
    end
end

on(sys1_menu.selection) do s
    if s == nothing 
	return
    else
	global system1 = s
    end
    updateSectionsMenu()
end
on(sys2_menu.selection) do s
    if s == nothing 
	return
    else
	global system2 = s
    end
    updateSectionsMenu()
end

on(sec_menu.selection) do s
    if s==nothing || s == "none" return end
    if system1 == "none"
	global section1_path = "none"
    else
	global section1_path = system1*"/"*s*"/E"
    end
    if system2 == "none" 
	global section2_path = "none"
    else
	global section2_path = system2*"/"*s*"/I"
    end
    redraw()
end

on(xaxis_menu.selection) do n
    global xaxis = n
    redraw()
end
on(yaxis_menu.selection) do n
    global yaxis = n
    redraw()
end

on(sec1_input.stored_string) do s
    global sec1 = parse(Int64,s)
    redraw()
end
on(sec2_input.stored_string) do s
    global sec2 = parse(Int64,s)
    redraw()
end

on(events(ax).mousebutton, priority = -1) do event
    if event.button == Mouse.left && event.action == Mouse.press
	plt,i = pick(ax)
	if plt isa Scatter
	    clipboard(string(last_hover))
	end
    end
end

DataInspector(ax)
fig
