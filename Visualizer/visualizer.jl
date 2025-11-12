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

@enum axisSel q1 q2 p1 p2 Q P
xaxis=1
yaxis=2

#buscar en data totes les combinacions
sections_list = ["none"]
sections_path = ["none"]
for folder in readdir("data")
    for mu in readdir("data/"*folder)
	path = "data/"*folder*"/"*mu
	if isdir(path*"/E")
	    push!(sections_list, folder*" "*mu*" Estable")
	    push!(sections_path,path*"/E")
	end
	if isdir(path*"/I")
	    push!(sections_list, folder*" "*mu*" Inestable")
	    push!(sections_path,path*"/I")
	end
    end
end

#TODO: Read data correcly to show

fig = Figure()

ax = Axis(fig[1,1])
fig[2,1] = inputGrid = GridLayout(tellwidth=false)

sections_options = zip(sections_list,sections_path)
section1_menu = Menu(fig,options=sections_options)
section2_menu = Menu(fig,options=sections_options)

xaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2","Q","P"],1:6), default = "q1")
yaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2","Q","P"],1:6), default = "q2")

input = inputGrid[1, 1:4] = [section1_menu, section2_menu, xaxis_menu,yaxis_menu]

scatter!(ax,1:10,1:10)
empty!(ax)

#drawing
function redraw()
    empty!(ax)
    if section1 != "none"
	draw_points(section1,"blue")
	#draw_lines(section1,"blue")
    end

    if section2 != "none"
	draw_points(section2,"red")
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
	    end
	    if yaxis <= 4
		push!(ypos,u[yaxis])
	    end
	end
	scatter!(ax,Point2f.(xpos,ypos),color=color,
	  inspector_label = (self, i ,pos) -> sec_file*"\nθ="*string(data[i].theta))
    end
end


#TODO:Fix line problem
function draw_lines(path,line_color)
    data = load_object(path)

    max_length = 0
    for (k,V) in data
	max_length = max(max_length, length(V))
    end

    for step=1:max_length
	xpos = []
	ypos = []
	for (k,V) in data
	    if length(V) >= step
		pos = V[step].u
		x = 0
		if xaxis <= 4
		    x = pos[xaxis]
		else
		    #TODO: Modulos
		end

		y = 0
		if yaxis <= 4
		    y = pos[yaxis]
		else
		    #TODO: Modulos
		end
		push!(xpos,x)
		push!(ypos,y)
	    elseif !isempty(xpos)
		#hi ha un forat
		if length(xpos) == 1
		    scatter!(ax,Point2f.(xpos,ypos),color=line_color)
		else
		    lines!(ax,Point2f.(xpos,ypos),color=line_color)
		end
		empty!(xpos)
		empty!(ypos)
	    end
	end
	if !isempty(xpos)
	    #juntar con inicio
	    (k,V) = first(data)
	    if length(V) >= step
		pos = V[step].u
		x = 0
		if xaxis <= 4
		    x = pos[xaxis]
		else
		    #TODO: Modulos
		end

		y = 0
		if yaxis <= 4
		    y = pos[yaxis]
		else
		    #TODO: Modulos
		end
		push!(xpos,x)
		push!(ypos,y)
	    end
	    lines!(ax,Point2f.(xpos,ypos),color=line_color)
	    empty!(xpos)
	    empty!(ypos)
	end
    end
end

#events
on(section1_menu.selection) do s
    global section1 = s
    redraw()
end
on(section2_menu.selection) do s
    global section2 = s
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
