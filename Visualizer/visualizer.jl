include("definitions.jl")
using DataStructures
using GLMakie
using JLD2

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
for (path, dirs, files) in walkdir("data")
    for file in files
	mu_string = replace(file,".jld2"=>"")
	mu = tryparse(Float64,mu_string)
	if(mu !== nothing)
	    push!(sections_path,path*"/"*file)
	    path_list = split(path,"/")
	    push!(sections_list, path_list[2]*" "*mu_string)
	end
    end
end
@show sections_list

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
    data = load_object(path)

    xpos = []
    ypos = []
    for (k,V) in data
	for v in V
	    pos = v.u
	    x = 0
	    if xaxis <= 4
		x = pos[xaxis]
	    elseif xaxis == 5
		x = pos[1]^2 + pos[2]^2
	    elseif xaxis == 6
		x = pos[3]^2 + pos[4]^2
	    end

	    y = 0
	    if yaxis <= 4
		y = pos[yaxis]
	    elseif yaxis == 5
		y = pos[1]^2 + pos[2]^2
	    elseif yaxis == 6
		y = pos[3]^2 + pos[4]^2
	    end
	    push!(xpos,x)
	    push!(ypos,y)
	end
    end

    scatter!(ax,Point2f.(xpos,ypos),color=color)
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

fig
