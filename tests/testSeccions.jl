using PlotlyJS
include("RCP3BP.jl")

mu = 0.4
s = Sistem(mu)

tspan = (0.0,32.0)

num_orbits = 2000
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

for (i,ang) in pairs(rangeI)
    #solve sistems
    orbit_terminateI = end_callback(s,ang)
    pI = problemI(s,ang,tspan)

    solI = solve(pI,Feagin14(),
    	     dtmax=0.01,
    	     save_end=false,
    	     callback = orbit_terminateI )
    
    
    dq1_prev = solI[2,1] + solI[3,1]
    dq2_prev = solI[4,1] - solI[1,1]
    for u in solI.u[2:end]
	dq1 = u[2] + u[3]
	dq2 = u[4] - u[1]

	if dq1*dq1_prev <= 0 && dq2_prev >= 0.0 && dq2 > 0.0
	    push!(q1I,u[1])
	    push!(q2I,u[2])
	    push!(p1I,u[3])
	    push!(p2I,u[4])
	    push!(angI,ang)
    
	end
	dq1_prev = dq1
	dq2_prev = dq2
    end
end

for (i,ang) in pairs(rangeE)
    #solve sistems
    orbit_terminateE = end_callback(s,ang,esEstable=true)
    pE = problemE(s,ang,tspan)

    solE = solve(pE,Feagin14(),
    	     dtmax=0.01,
    	     save_end=false,
    	     callback = orbit_terminateE )

    dq1_prev = solE[2,1] + solE[3,1]
    dq2_prev = solE[4,1] - solE[1,1]
    for u in solE.u[2:end]
	dq1 = u[2] + u[3]
	dq2 = u[4] - u[1]

	#if dq1*dq1_prev <= 0 && dq2_prev <= 0.0 && dq2 < 0.0
	if dq1*dq1_prev <= 0 && dq2_prev >= 0.0 && dq2 > 0.0
	    push!(q1E,u[1])
	    push!(q2E,u[2])
	    push!(p1E,u[3])
	    push!(p2E,u[4])
	    push!(angE,ang)

    
	end
	dq1_prev = dq1
	dq2_prev = dq2
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
