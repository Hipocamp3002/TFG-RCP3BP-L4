using DifferentialEquations
#using Plots
using PlotlyJS
using LinearAlgebra
include("RCP3BP.jl")

mu = 0.4

s = Sistem(mu)

#theta = 5.1001550005 #Viatja a L5
#theta = 5.178331 #torna a L4
theta = 3.0

tspan = (0.0,102)
prob = problem(s,theta,tspan)
sol = solve(prob,Feagin14(),dtmax=0.01)


Qplot = plot([
    scatter(x=sol[1,:],y=sol[2,:],name="orbit",mode="lines"),
    scatter(x=[s.L4[1],s.L4[1]],y=[s.L4[2],-s.L4[2]],name="Lagrange",mode="markers"),
    scatter(x=[-mu,1-mu],y=[0,0],name="bodies",mode="markers")
    ],Layout(title="Q"))

Pplot = plot([
    scatter(x=sol[3,:],y=sol[4,:],name="orbit",mode="lines"),
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
