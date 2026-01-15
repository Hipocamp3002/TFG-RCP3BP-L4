include("RCP3BP.jl")
using PlotlyJS

mu = 0.4
s = Sistem(mu)

theta = 5.2


cond(u,t,integrator) = u[3]+u[2]
function affect!(integrator)
    u = integrator.u
    if(u[4] - u[1] > 0.0)
	res = savevalues!(integrator, true)
    end
    #@show res
end
cb = ContinuousCallback(cond,affect!,
			save_positions = (false,false))

tspan = (0.0,400)
prob = problem(s,theta,tspan)
integrator = init(prob,Feagin14(),
		  dtmax = 0.01, 
		  callback = cb,
		  save_everystep = false,
		  save_start = false,
		  save_end = false)

H_ini =  H0(s,theta) 
for (u,i) in tuples(integrator)
    H = hamiltonian((u[3],u[4]),(u[1],u[2]),s.mu)
    if(1e-5 < abs(H_ini - H))
	break
    end
end


sol = integrator.sol


Qplot = plot([
    scatter(x=sol[1,:],y=sol[2,:],name="orbit",mode="markers"),
    scatter(x=[s.L4[1],s.L4[1]],y=[s.L4[2],-s.L4[2]],name="Lagrange",mode="markers"),
    scatter(x=[-mu,1-mu],y=[0,0],name="bodies",mode="markers")
    ],Layout(title="Q"))

Pplot = plot([
    scatter(x=sol[3,:],y=sol[4,:],name="orbit",mode="markers"),
    scatter(x=[s.L4[3],-s.L4[3]],y=[s.L4[4],s.L4[4]],name="Lagrange",mode="markers")
    ],Layout(title="P"))

H = []
for u in sol.u
    push!(H,hamiltonian((u[3],u[4]),(u[1],u[2]),(mu)))
end
H = H .- H[1]

Hplot = plot(scatter(x=1:length(H),y=H,mode="lines"),Layout(title="H diff"))

p = [Qplot Pplot;Hplot]
relayout!(p,showlegend=false)
p
