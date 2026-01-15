include("RCP3BP.jl")
using PlotlyJS

mu = 0.1
s = Sistem(mu)

tspan = (0.0,50.0)

num_orbits = 500
rangeI = collect(LinRange(0,2π,num_orbits))
rangeE = collect(LinRange(0,2π,num_orbits))
#range = [5.3]

q1I = []
q2I = []
p1I = []
p2I = []
angI = []

q1E = []
q2E = []
p1E = []
p2E = []
angE = []


Htest = getHamiltonianTest(s, error=1e-10)

#eqCond(u) = u[1] - 0.5 + mu
#dirCond(u) = (u[2] + u[3]) > 0.0

eqCond(u) = u[1]
dirCond(u) = (u[2] + u[3]) > 0.0


for (i,ang) in pairs(rangeI)
    prob = problemI(s,ang,tspan)
    cb = seccio_callback(eqCond,dirCond)
    
    integrator = init(prob,Feagin14(),
		  dtmax = 0.01, 
		  callback = cb,
		  save_everystep = false,
		  save_start = false,
		  save_end = false)

    for (u,i) in tuples(integrator)
	if (Htest(u))
	    break
	end
    end

    sol = integrator.sol

    l = length(sol.u)
    if l > 0
	append!(q1I,sol[1,:])
	append!(q2I,sol[2,:])
	append!(p1I,sol[3,:])
	append!(p2I,sol[4,:])
	append!(angI,fill(ang,l))
    end
end

for (i,ang) in pairs(rangeE)
    prob = problemE(s,ang,tspan)
    cb = seccio_callback(eqCond,dirCond)
    
    integrator = init(prob,Feagin14(),
		  dtmax = 0.01, 
		  callback = cb,
		  save_everystep = false,
		  save_start = false,
		  save_end = false)

    for (u,i) in tuples(integrator)
	if (Htest(u))
	    break
	end
    end

    sol = integrator.sol
    l = length(sol.u)
   if l > 0
	append!(q1E,sol[1,:])
	append!(q2E,sol[2,:])
	append!(p1E,sol[3,:])
	append!(p2E,sol[4,:])
	append!(angE,fill(ang,l))
    end
end

q1q2 = plot([scatter(x=q1I,y=q2I,name="Inestable",text = angI,mode="markers"),
    scatter(x=q1E,y=q2E,name="Estable",text = angE,mode="markers")],Layout(title="q1q2"))
q1p1 = plot([scatter(x=q1I,y=p1I,name="Inestable",text = angI,mode="markers"),
    scatter(x=q1E,y=p1E,name="Estable",text = angE,mode="markers")],Layout(title="q1p1"))
q1p2 = plot([scatter(x=q1I,y=p2I,name="Inestable",text = angI,mode="markers"),
    scatter(x=q1E,y=p2E,name="Estable",text = angE,mode="markers")],Layout(title="q1p2"))
q2p1 = plot([scatter(x=q2I,y=p1I,name="Inestable",text = angI,mode="markers"),
    scatter(x=q2E,y=p1E,name="Estable",text = angE,mode="markers")],Layout(title="q2p1"))
q2p2 = plot([scatter(x=q2I,y=p2I,name="Inestable",text = angI,mode="markers"),
    scatter(x=q2E,y=p2E,name="Estable",text = angE,mode="markers")],Layout(title="q2p2"))
p1p2 = plot([scatter(x=p1I,y=p2I,name="Inestable",text = angI,mode="markers"),
    scatter(x=p1E,y=p2E,name="Estable",text = angE,mode="markers")],Layout(title="p1p2"))

QI =sqrt.(q1I.^2 + q2I.^2)
QE =sqrt.(q1E.^2 + q2E.^2)

PI =sqrt.(p1I.^2 + p2I.^2)
PE =sqrt.(p1E.^2 + p2E.^2)

#QP = plot([scatter(x=QI,y=PI,name="Inestable",text = angI,mode="markers"),
#    scatter(x=QE,y=PE,name="Estable",text = angE,mode="markers")],Layout(title="QP"))


plt = [q1q2 q1p1 q1p2; q2p1 q2p2 p1p2]
relayout!(plt,showlegend=false)
plt
