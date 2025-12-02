using GLMakie
using OrdinaryDiffEq
include("systems.jl")

fig = Figure()

ax = Axis(fig[2,1], aspect = AxisAspect(1), autolimitaspect = 1)
fig[1,1] = inputGrid_orbit = GridLayout(tellwidth=false)
fig[3,1] = inputGrid_axis = GridLayout(tellwidth=false)

angle_input = Textbox(fig,placeholder = "Enter angle",validator = Float64)
mu_input = Textbox(fig,placeholder = "Enter mu",validator = Float64)
time_input = Textbox(fig,placeholder = "Enter time",validator = Float64)
calc_button = Button(fig,label="calculate")


xaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2"],1:4), default = "q1")
yaxis_menu = Menu(fig,options=zip(["q1","q2","p1","p2"],1:4), default = "q2")

input_orbit = inputGrid_orbit[1,1:4] = [angle_input,mu_input,calc_button,time_input]
input_axis = inputGrid_axis[1,1:2] = [xaxis_menu,yaxis_menu]

mu = 0.0
angle = 0.0
time = 100.0
x_axis = 1
y_axis = 2

path = Observable(Vector{Vector{Float64}}(undef,0))
path_points = lift(path) do P
    points = Point2f[]
    for u in P
	push!(points,Point2f(u[x_axis],u[y_axis]))
    end
return points
end

bodies = Observable([[-mu,0,0,0],[1-mu,0,0,0]])
bodies_points = lift(bodies) do u
    Point2f.([u[1][x_axis],u[2][x_axis]],[u[1][y_axis],u[2][y_axis]])
end

lagrange = Observable([0.5-mu sqrt(3)/2 -sqrt(3)/2 0.5-mu ;
		      0.5-mu -sqrt(3)/2 sqrt(3)/2 0.5-mu])
lagrange_points = lift(lagrange) do L
    Point2f.(L[:,x_axis],L[:,y_axis])
end

on(mu_input.stored_string) do s
    global mu = parse(Float64,s)
    bodies[][1] = [-mu,0,0,0]
    bodies[][2] = [1-mu,0,0,0]
    lagrange[][1,:] = [0.5-mu, sqrt(3)/2, -sqrt(3)/2, 0.5-mu]
    lagrange[][2,:] = [0.5-mu, -sqrt(3)/2, sqrt(3)/2, 0.5-mu]
    empty!(path[])
    notify(path)
    notify(bodies)
    notify(lagrange)
end

on(angle_input.stored_string) do s
    global angle = parse(Float64,s)
    empty!(path[])
    calc_orbit()
    notify(path)
end

on(time_input.stored_string) do s
    global time = parse(Float64,s)
end

on(xaxis_menu.selection) do v
    global x_axis = v
    notify(bodies)
    notify(lagrange)
    notify(path)
end

on(yaxis_menu.selection) do v
    global y_axis = v
    notify(bodies)
    notify(lagrange)
    notify(path)
end


function max_hamiltonian()
    H_ini = hamiltonian(path[][1],mu)
    max_H = 0.0
    for u in path[]
	cur_H = hamiltonian(u,mu)
	max_H = max(max_H, abs(cur_H - H_ini))
    end
    return max_H
end

function unstable_f(H_ini,err)
    return (dt,u,p,t) -> err < abs(H_ini - hamiltonian(u,mu))
end

function calc_orbit()
    s = sistema_L4(mu)
    
    u0 = s.center + s.eps*(cos(angle)*s.Ivecs[1] + sin(angle)*s.Ivecs[2])

    u0sa = SVector{4,Float64}(u0)
    prob = ODEProblem(orbit,u0sa,(0.0,time),(mu))
    sol = solve(prob,Vern9(),
		abstol = 1e-14,reltol = 1e-14,
		unstable_check = unstable_f(hamiltonian(u0,mu),1e-10))
    append!(path[],[collect(x) for x in sol.u])
end

on(calc_button.clicks) do n
    empty!(path[])
    calc_orbit()
    notify(path)
end

scatter!(ax,bodies_points)
scatter!(ax,lagrange_points)
lines!(ax,path_points)

fig
