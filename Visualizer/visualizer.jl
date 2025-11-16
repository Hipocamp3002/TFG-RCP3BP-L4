include("definitions.jl")
using DataStructures
using GLMakie
using JLD2
using CSV

#Plot with selectors
#   - mu value
#   - 2 sections
#   - Axis selector



section1 = "none"
section2 = "none"

section1_path = "none"
section2_path = "none"

xaxis=1
yaxis=2

#buscar en data totes les combinacions
sections = Dict{String,Tuple{String,Array{String}}}("none"=> ("none",[]))
for folder in readdir("data")
    for mu in readdir("data/"*folder)
	path = "data/"*folder*"/"*mu
	if isdir(path*"/E")
	    mu_path = path*"/E"
	    entry = folder*" Estable"
	    if haskey(sections, entry)
		push!(sections[entry][2], mu)
	    else
		push!(sections, entry => ("data/"*folder*"/{mu}/E",[mu]))
	    end

	end
	if isdir(path*"/I")
	    mu_path = path*"/I"
	    entry = folder*" Inestable"
	    if haskey(sections, entry)
		push!(sections[entry][2], mu)
	    else
		push!(sections, entry => ("data/"*folder*"/{mu}/I",[mu]))
	    end
	end
    end
end


fig = Figure()

ax = Axis(fig[2,1])
fig[1,1] = inputGrid_sec = GridLayout(tellwidth=false)
fig[3,1] = inputGrid_axis = GridLayout(tellwidth=false)

sections_options = collect(keys(sections))
sort!(sections_options)
section1_menu = Menu(fig,options=sections_options, default="none")
section2_menu = Menu(fig,options=sections_options, default="none")
mu_menu = Menu(fig,options=["none"])

xaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2","θ"],1:5), default = "q1")
yaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2","θ"],1:5), default = "q2")



input_sec = inputGrid_sec[1,1:3] = [section1_menu, section2_menu, mu_menu]
input_axis = inputGrid_axis[1,1:2] = [xaxis_menu,yaxis_menu]

scatter!(ax,1:10,1:10)
empty!(ax)

function updateMuMenu()
    if section1 == "none" && section2 == "none"
	mu_menu.options = ["none"]
    elseif section1 == "none"
	mus = sections[section2][2]
	mu_menu.options = mus
    elseif section2 == "none"
	mus = sections[section1][2]
	mu_menu.options = mus
    else
	mu1 = sections[section1][2]
	mu2 = sections[section2][2]
	
	mus =  intersect(mu1,mu2)
	mu_menu.options = mus
    end
end

#drawing
function redraw()
    empty!(ax)
    if section1 != "none"
	draw_points(section1_path,"blue")
    end

    if section2 != "none"
	draw_points(section2_path,"red")
    end
end

function draw_points(path,color)
    for sec_file in readdir(path)
	xpos = []
	ypos = []
	data = CSV.File(path*"/"*sec_file)
	for (uStr,t,theta) in data
	    u = parse.(Float64, split(chop(uStr,head=1,tail=1),','))
	    if xaxis <= 4
		push!(xpos,u[xaxis])
	    else xaxis == 5
		push!(xpos,theta)
	    end
	    if yaxis <= 4
		push!(ypos,u[yaxis])
	    else yaxis == 5
		push!(ypos,theta)
	    end
	end
	scatter!(ax,Point2f.(xpos,ypos),color=color,
	  inspector_label = (self, i ,pos) -> sec_file*"\nθ="*string(data[i].theta))
    end
end


#events
on(section1_menu.selection) do s
    global section1 = s
    updateMuMenu()
end
on(section2_menu.selection) do s
    global section2 = s
    updateMuMenu()
end
on(mu_menu.selection) do s
    if s == nothing return end
    if section1 == "none"
	global section1_path = "none"
    else
	path = sections[section1][1]
	first, second = split(path,"{mu}")
	section1_path = first*s*second
    end
    if section2 == "none"
	global section2_path = "none"
    else
	path = sections[section2][1]
	first, second = split(path,"{mu}")
	global section2_path = first*s*second
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

DataInspector()
fig
