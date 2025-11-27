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

section = "none"

section1_path = "none"
section2_path = "none"

xaxis=1
yaxis=2

last_hover=0.0

curr_mu = 0.0


fig = Figure()

ax = Axis(fig[2,1])
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

xaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2","θ","Q_L4_angle","P_L4_angle"],1:7)
		  , default = "q1")
yaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2","θ","Q_L4_angle","P_L4_angle"],1:7)
		  , default = "q2")


input_sec = inputGrid_sec[1,1:4] = [mu_menu,sys1_menu, sys2_menu, sec_menu]
input_axis = inputGrid_axis[1,1:2] = [xaxis_menu,yaxis_menu]

Splt = scatter!(ax,1:10,1:10)
empty!(ax)

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
    empty!(ax)
    if section1_path != "none"
	draw_points(section1_path,"blue")
    end

    if section2_path != "none"
	draw_points(section2_path,"red")
    end
end

function inspector(sec_file,data)
    function inspect(self,i,pos)
	global last_hover = data[i].theta
	sec_file*"\nθ="*string(data[i].theta)
    end
    return inspect
end

function draw_points(path,color)
    for sec_file in readdir(path)
	section = parse(Int64,split(sec_file,".")[1])
	xpos = []
	ypos = []
	data = CSV.File(path*"/"*sec_file)
	for (uStr,t,theta) in data
	    u = parse.(Float64, split(chop(uStr,head=1,tail=1),','))
	    if xaxis <= 4
		push!(xpos,u[xaxis])
	    elseif xaxis == 5
		push!(xpos,theta)
	    elseif xaxis == 6
		angle = atan(u[1]-0.5+curr_mu, u[2]-sqrt(3)/2)
		push!(xpos,angle)
	    elseif xaxis == 7
		angle = atan(u[3]+sqrt(3)/2, u[4]-0.5+curr_mu)
		push!(xpos,angle)
	    end
	    if yaxis <= 4
		push!(ypos,u[yaxis])
	    elseif yaxis == 5
		push!(ypos,theta)
	    elseif yaxis == 6
		angle = atan(u[1]-0.5+curr_mu, u[2]-sqrt(3)/2)
		push!(ypos,angle)
	    elseif yaxis == 7
		angle = atan(u[3]+sqrt(3)/2, u[4]-0.5+curr_mu)
		push!(ypos,angle)
	    end
	end
	scatter!(ax,Point2f.(xpos,ypos),color=color,
	  inspector_label = inspector(sec_file,data))
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

on(events(Splt).mousebutton, priority = -1) do event
    if event.button == Mouse.left && event.action == Mouse.press
	plt,i = pick(Splt)
	if plt isa Scatter
	    clipboard(string(last_hover))
	end
    end
end

DataInspector()
fig
