#using DifferentialEquations
#using DiffEqCallbacks
using PlotlyJS
#using LinearAlgebra
include("RCP3BP.jl")

mu = 0.4
thetaE = 5.673779374376215
thetaI = 5.8418913558392145
tspan = (0.0,100.0)

s = Sistem(mu)


orbit_terminateI = end_callback(s,thetaI)
orbit_terminateE = end_callback(s,thetaE,esEstable=true)
pI = problemI(s,thetaI,tspan)
pE = problemE(s,thetaE,tspan)

solI = solve(pI,Feagin14(),
	     dtmax=0.01,
	     save_end=false,
	     callback = orbit_terminateI )

solE = solve(pE,Feagin14(),
	     dtmax=0.01,
	     save_end=false,
	     callback = orbit_terminateE )


Qplot = plot([scatter(x=solI[1,:], y=solI[2,:],name="Inestable",mode="lines"),
	scatter(x=solE[1,:], y=solE[2,:],name="Estable",mode="lines")])
Pplot = plot([scatter(x=solI[3,:], y=solI[4,:],name="Inestable",mode="lines"),
	scatter(x=solE[3,:], y=solE[4,:],name="Estable",mode="lines")])

plt = [Qplot;Pplot]
relayout!(plt,showlegend=false)
plt
