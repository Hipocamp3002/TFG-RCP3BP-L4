using GLMakie
using OrdinaryDiffEq
include("systems.jl")
include("LeviCivita.jl")

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

path = Observable(Vector{SVector{4,Float64}}(undef,0))
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



function orbit_with_coll(U,p,t)::SVector{5,Float64}
    mu,E = p
    u = SVector{4}(U[1:4])
    coll = U[5]
    if coll == 2.0
	res = LCorbitP2(u,(mu,E),t) 
	return [res ;0.0]
    elseif coll == 1.0
	res = LCorbitP1(u,(mu,E),t)
	return [res ;0.0]
    else
	res = orbit(u,mu,t)
	return [res ; 0.0]
    end
    SA[0.0,0.0,0.0,0.0,0.0]
end

min_coll = 0
function collCondition(u,_,integrator)
    c = integrator.u[5]
    mu = integrator.p[1]
    if c == 0
	q1 = u[1]
	q2 = u[2]
	coll1 = (q1+mu)^2 + q2^2 - 0.01
	coll2 = (q1+mu-1)^2 + q2^2 - 0.01
	if abs(coll1) < abs(coll2)
	    global min_coll = 1
	    return coll1
	else
	    global min_coll = 2
	    return coll2
	end
    elseif c == 1
	u1 = u[1]
	u2 = u[2]
	coll1 = (u1^2 + u2^2)^2 - 0.01
	global min_coll = 1
	return coll1
    elseif c == 2
	u1 = u[1]
	u2 = u[2]
	coll2 = (u1^2 + u2^2)^2 - 0.01
	global min_coll = 2
	return coll2
    end
end

function collAffect(integrator)
    U = integrator.u
    u = SA[U[1],U[2],U[3],U[4]]
    c = integrator.u[5]
    mu = integrator.p[1]
    if min_coll == 1
	if c == 0
	    q1 = u[1] + mu
	    q2 = u[2]
	    p1 = u[3]
	    p2 = u[4] + mu
	    if q1*(p1+q2) + q2*(p2-q1) < 0
		integrator.u = SVector{5}([qp2uU(SA[q1,q2,p1,p2]);1])
	    end
	elseif c == 1
	    u1,u2,U1,U2 = u
	    if (U1*u1 + U2*u2)/2 > 0
		u_new = uU2qp(u) + [-mu,0.0,0.0,-mu]
		integrator.u = SVector{5}([u_new ; 0])
	    end

	end
    elseif min_coll == 2
	if c == 0
	    q1 = u[1] - 1 + mu
	    q2 = u[2]
	    p1 = u[3]
	    p2 = u[4] - 1 + mu
	    if q1*(p1+q2) + q2*(p2-q1) < 0
		integrator.u = SVector{5}([qp2uU(SA[q1,q2,p1,p2]);2])
	    end
	elseif c == 2
	    u1,u2,U1,U2 = u
	    if (U1*u1 + U2*u2)/2 > 0
		u_new = uU2qp(u) + [1-mu,0,0,1-mu]
		integrator.u = SVector{5}([u_new ; 0])
	    end
	end
    end
end

cb = ContinuousCallback(collCondition,collAffect)

function sol2pos(U)
    u = U[1:4]
    c = U[5]
    if c == 1
	return uU2qp(u) - [mu,0,0,mu]
    elseif c== 2
	return uU2qp(u) + [1-mu,0,0,1-mu]
    end
    return u
end

function calc_orbit()
    s = sistema_L4(mu)
    u0 = s.center + s.eps*(cos(angle)*s.Ivecs[1] + sin(angle)*s.Ivecs[2])
    H0 = hamiltonian(s.center,mu)
    u0sa = SVector{4,Float64}(u0)
    
    prob = ODEProblem(orbit_with_coll,SVector{5}([u0sa;0]),(0.0,time),[mu,H0])
	sol = solve(prob,Vern9(),
	    	abstol = 1e-14,reltol = 1e-14,
		dense = false,
		callback = cb,
	    	#isoutofdomain = unstable_f(H0,1e-10)
	     )
    append!(path[],sol2pos.(sol.u))
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
