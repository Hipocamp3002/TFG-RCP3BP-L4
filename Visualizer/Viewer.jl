using GLMakie
using GZip
using StringViews
using NativeFileDialog

include("ViewerOrbit.jl")
sections = include("ViewerSections.jl")

fig = Figure()

lg = fig[1,1] = GridLayout()
rg = fig[1,2] = GridLayout(tellheight=false)
data_grid = lg[1,1] = GridLayout(tellwidth=false)
data_button = Button(data_grid[1,1],label="directori de dades")
data_label = Label(data_grid[1,2],"./data")
Box(data_grid[1,3],visible = false)
axq = Axis(lg[2,1],autolimitaspect = 1, title="Posició (q)")
axp = Axis(lg[3,1],autolimitaspect=1, title="Moment (p)")
mu_sl_grid = lg[4,1] = GridLayout()
mu_label = Label(mu_sl_grid[1,1],"0.25",width=50)
mu_slider = Slider(mu_sl_grid[1,2], range = 0.0386:0.0001:0.5,startvalue = 0.25)

sec_grid = rg[1,1] = GridLayout()
Label(sec_grid[1,1],"Secció:")
sec_menu = Menu(sec_grid[1,2],options = ["dq1_ltdq2","dq1_gtdq2","dq2_ltdq1","dq2_gtdq1",
					 "dp1_ltdp2","dp1_gtdp2","dp2_ltdp1","dp2_gtdp1"])

ori_grid = rg[2,1] = GridLayout()
Label(ori_grid[1,1:2],"Origen de la variant")
Label(ori_grid[2,1],"Estable")
Label(ori_grid[3,1],"Inestable")
oriE_menu = Menu(ori_grid[2,2],options = ["L4","L5"])
oriI_menu = Menu(ori_grid[3,2],options = ["L4","L5"])

int_grid = rg[3,1] = GridLayout(tellwidth=false)
Label(int_grid[1,1],"Interseccio a visualitzar",tellwidth=false)
int_sliders = SliderGrid(int_grid[2,1],
			    (label = "Estable", range = 1:16, startvalue = 4,
			 color_active_dimmed = :royalblue,color_active=:blue),
			    (label = "Inestable", range = 1:16, startvalue = 4,
			 color_active_dimmed = :firebrick,color_active=:red)
			 )

time_grid = rg[4,1] = GridLayout(tellwidth=false)
Label(time_grid[1,1],"Duracio orbita",tellwidth=false)
time_input = Textbox(time_grid[1,2],validator=Float64,
		     stored_string = "30",width=200)

Label(rg[5,1],"Mida punts (pixels)",tellwidth = false)
size_slider = Slider(rg[6,1],range = 0.1:0.1:10,startvalue = 5,tellwidth = false)

mu_grid = rg[7,1] = GridLayout(tellwidth=false)
Label(mu_grid[1,1],"µ estable")
Label(mu_grid[2,1],"µ inestable")
Label(mu_grid[3,1],"µ mitjana")
muE_label = Label(mu_grid[1,2],"da")
muI_label = Label(mu_grid[2,2],"ad")
muR_label = Label(mu_grid[3,2],"ad")
Box(mu_grid[1:3,3],visible=false)

point_grid = rg[8,1] = GridLayout(tellwidth=false)
Label(point_grid[1,1],"q_x:")
Label(point_grid[2,1],"q_y:")
Label(point_grid[3,1],"p_x:")
Label(point_grid[4,1],"p_y:")
q1_label = Label(point_grid[1,2],"~")
q2_label = Label(point_grid[2,2],"~")
p1_label = Label(point_grid[3,2],"~")
p2_label = Label(point_grid[4,2],"~")
Box(point_grid[1:4,3],visible=false)

Box(rg[9,1],visible=false)#separator

move_grid = rg[10,1] = GridLayout(tellwidth=false)
step_input = Textbox(move_grid[1,1],validator=Float64, stored_string ="0.01",width=300-128-16)
add_button = Button(move_grid[1,2],label="+",width=64)
sub_button = Button(move_grid[1,3],label="-",width=64)
set_input = Textbox(move_grid[2,1],validator=Float64, stored_string ="0.25",width=300-128-16)
set_button = Button(move_grid[2,2:3],label="=",width=128)


colsize!(fig.layout,2,300)

#Dades

data = "data"

E_inter = int_sliders.sliders[1].value
I_inter = int_sliders.sliders[2].value

E_mu_list = []
I_mu_list = []

cur_mu = 0.25
E_mu = -1
I_mu = -1
real_mu = 0.25

step_mu = @lift(parse(Float64, $(step_input.stored_string)))

time = @lift(max(0.0,parse(Float64, $(time_input.stored_string))))

orbit_point = [0.0,0.0,0.0,0.0]

Eq_points = Observable(Point2f[])
Ep_points = Observable(Point2f[])
Iq_points = Observable(Point2f[])
Ip_points = Observable(Point2f[])
orbitq_points = Observable(Point2f[])
orbitp_points = Observable(Point2f[])
orbit_color = Observable(Float64[])
orbit_iniq = Observable(Point2f[])
orbit_inip = Observable(Point2f[])

planets_q = Observable(Point2f[])
planets_p = Observable(Point2f[])
function set_planets()
    mu = real_mu
    planets_q[] = [Point2f(-mu,0),
	Point2f(1-mu,0),
	Point2f(0.5-mu,sqrt(3)/2),
	Point2f(0.5-mu,-sqrt(3)/2)]
    planets_p[] = [Point2f(0.0,-mu),
	Point2f(0.0,1-mu),
	Point2f(-sqrt(3)/2,0.5-mu),
	Point2f(sqrt(3)/2,0.5-mu)]

    notify(planets_q)
    notify(planets_p)
end
set_planets()

point_size = Observable(5)

scatter!(axq,Eq_points,color="blue",markersize = size_slider.value)
#datashader!(axq,Eq_points,colormap=[(:white,0.0),(:blue,1.0)],binsize=2)
scatter!(axp,Ep_points,color="blue",markersize = size_slider.value)
scatter!(axq,Iq_points,color="red",markersize = size_slider.value)
#datashader!(axq,Iq_points,colormap=[(:white,0.0),(:red,1.0)],binsize=2)
scatter!(axp,Ip_points,color="red",markersize = size_slider.value)
scatter!(axq,planets_q,marker = [:circle,:circle,:star4,:star4], color = "black")
scatter!(axp,planets_p,marker = [:circle,:circle,:star4,:star4], color = "black")
lines!(axq,orbitq_points,color = orbit_color,colormap = [:cyan3, :magenta3])
lines!(axp,orbitp_points,color = orbit_color,colormap = [:cyan3, :magenta3])
scatter!(axq,orbit_iniq,color="green4",marker = :cross, markersize = 12)
scatter!(axp,orbit_inip,color="green4",marker = :cross, markersize = 12)

#Interaccions
on(mu_slider.value) do v
    global cur_mu = v
    mu_label.text = string(cur_mu)[1:min(end,6)]
    set_mu()
end

on(sec_menu.selection) do _ reload() end
on(oriE_menu.selection) do _ reload() end
on(oriI_menu.selection) do _ reload() end
on(E_inter) do _ reload() end
on(I_inter) do _ reload() end

on(add_button.clicks) do _
    global cur_mu += step_mu[]
    if cur_mu > 0.5 global cur_mu = 0.5 end
    mu_label.text = string(cur_mu)[1:min(end,6)]
    set_mu()
end
on(sub_button.clicks) do _
    global cur_mu -= step_mu[]
    if cur_mu < routh global cur_mu = routh end
    mu_label.text = string(cur_mu)[1:min(end,6)]
    set_mu()
end
on(set_button.clicks) do _
    mu = parse(Float64,set_input.stored_string[])
    global cur_mu = clamp(mu,routh,0.5)
    mu_label.text = string(cur_mu)[1:min(end,6)]
    set_mu()
end

on(data_button.clicks) do _
    d = pick_folder()
    global data = d
    data_label.text = d
end

#Q axis interaction
register_interaction!(axq, :my_interaction) do event::MouseEvent, axis
    if event.type === MouseEventTypes.leftclick
	get_point_q()
    end
end

on(events(fig).keyboardbutton) do k
    if k.action == Keyboard.release
	if k.key == Keyboard.right
	    add_button.clicks[] += 1
	elseif k.key == Keyboard.left
	    sub_button.clicks[] += 1
	end
    end
end

#Accions

function find_mu()
    empty!(E_mu_list)
    empty!(I_mu_list)

    if !isdir(data)
	return
    end

    sec = sec_menu.selection[]
    Eori = oriE_menu.selection[]
    Iori = oriI_menu.selection[]

    Eint = string(E_inter[])
    Iint = string(I_inter[])

    for mu_file in readdir(data)
	mu = parse(Float64,mu_file)
	#Estable
	pathE = data*"/"*mu_file*"/"*Eori*"/"*sec*"/E/"*Eint
	if isfile(pathE)
	    push!(E_mu_list,mu)
	end
	#Inestable
	pathI = data*"/"*mu_file*"/"*Iori*"/"*sec*"/I/"*Iint
	if isfile(pathI)
	    push!(I_mu_list,mu)
	end
    end
    sort!(E_mu_list)
    sort!(I_mu_list)
end

function set_mu(force_load = false)
    #empty everything
    empty!(orbitq_points[])
    empty!(orbitp_points[])
    empty!(orbit_iniq[])
    empty!(orbit_inip[])
    empty!(orbit_color[])
    notify(orbitq_points)
    notify(orbitp_points)
    notify(orbit_iniq)
    notify(orbit_inip)
    q1_label.text="~"
    q2_label.text="~"
    p1_label.text="~"
    p2_label.text="~"

    new_E_mu = find_nearest_sorted(E_mu_list,cur_mu)
    new_I_mu = find_nearest_sorted(I_mu_list,cur_mu)
    if new_E_mu == nothing
	empty!(Eq_points[])
	empty!(Ep_points[])
	muE_label.text = "Cap"
	notify(Eq_points)
	notify(Ep_points)
    elseif new_E_mu != E_mu
	global E_mu = new_E_mu
	load_points(Eq_points,Ep_points,E)
	muE_label.text = string(new_E_mu)
    elseif force_load
	load_points(Eq_points,Ep_points,E)
	muE_label.text = string(new_E_mu)
    end
    if new_I_mu == nothing
	empty!(Iq_points[])
	empty!(Ip_points[])
	muI_label.text[] = "Cap"
	notify(Iq_points)
	notify(Ip_points)
    elseif new_I_mu != I_mu
	global I_mu = new_I_mu
	load_points(Iq_points,Ip_points,I)
	muI_label.text = string(new_I_mu)
    elseif force_load
	load_points(Iq_points,Ip_points,I)
	muI_label.text = string(new_I_mu)
    end

    if new_E_mu == nothing && new_I_mu == nothing
	global real_mu = cur_mu
    elseif new_E_mu == nothing
	global real_mu = I_mu
    elseif new_I_mu == nothing
	global real_mu = E_mu
    else
	global real_mu = (I_mu+E_mu)/2
    end
    muR_label.text = string(real_mu)
    set_orbit_system(real_mu)
    set_planets()
end

function find_nearest_sorted(a,x)
    m = 1
    M = length(a)
    if M == 0
	return nothing
    end
    while M-m > 1
	b = div(M+m,2)
	if a[b] > x
	    M = b
	elseif a[b] < x
	    m = b
	else
	    return a[b]
	end
    end

    dm = abs(a[m]-x)
    dM = abs(a[M]-x)
    if dm > dM
	return a[M]
    end
    return a[m]
end

function sec_path(type::orbit_type)
    mu_file = string(type == E ? E_mu : I_mu)
    ori = type == E ? oriE_menu.selection[] : oriI_menu.selection[]
    sec = sec_menu.selection[]
    type_file = type == E ? "/E/" : "/I/"
    inter = string(type == E ? E_inter[] : I_inter[])

    return data*"/"*mu_file*"/"*ori*"/"*sec*type_file*inter
end

function load_points(q_points,p_points,type::orbit_type)
    empty!(q_points[])
    empty!(p_points[])
    path = sec_path(type)
    if !isfile(path)
	notify(q_points)
	notify(p_points)
	return
    end

    #TODO: preallocate
    empty!(q_points[])
    empty!(p_points[])
    u = Vector{Float64}(undef,4)
    GZip.open(path,"r") do f
	while !eof(f)
	    read!(f,u)
	    push!(q_points[],Point2f(u[1],u[2]))
	    push!(p_points[],Point2f(u[3],u[4]))
	end
    end

    notify(q_points)
    notify(p_points)
end

function reload()
    find_mu()
    set_mu(true)
end

#Orbit generation
function get_point_q()
    sec = sec_menu.selection[]
    if !haskey(sections,sec)
	return
    end
    if sections[sec].q_solve == nothing
	return
    end

    q1,q2 = mouseposition(axq.scene)
    p1,p2 = sections[sec].q_solve((q1,q2),real_mu,H_0)
    global orbit_point =  [q1,q2,p1,p2]

    q1_label.text = string(q1)
    q2_label.text = string(q2)
    p1_label.text = string(p1)
    p2_label.text = string(p2)

    get_orbits()
end

function get_orbits()
    path, mid = calc_orbit(orbit_point,real_mu,time[])
    p = path[mid]
    len = length(path)
    orbit_iniq[] = [Point2f(p[1],p[2])]
    orbit_inip[] = [Point2f(p[3],p[4])]
    orbit_color[] = [i < mid ? 0 : (i == mid ? 0.5 : 1) for i in 1:len]
    orbitq_points[] = [Point2f(u[1],u[2]) for u in path]
    orbitp_points[] = [Point2f(u[3],u[4]) for u in path]
    

end

#init functions
reload()
display(fig) |> wait
