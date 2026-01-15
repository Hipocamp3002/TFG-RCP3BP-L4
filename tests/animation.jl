using Plots
include("RCP3BP.jl")

@gif for mu=0.1:0.01:0.5
    s = Sistem(mu)

    tspan = (0.0,32.0)

    num_orbits = 500
    range = collect(LinRange(0,2π,num_orbits))
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

    for (i,ang) in pairs(range)
    #solve sistems
	orbit_terminateI = end_callback(s,ang)
	pI = problemI(s,ang,tspan)

	pE = problemE(s,ang,tspan)#L5
	orbit_terminateE = end_callback(pE,mu)

	solI = solve(pI,Feagin14(),
		dtmax=0.01,
		save_end=false,
		callback = orbit_terminateI )
    
	solE = solve(pE,Feagin14(),
		dtmax=0.01,
		save_end=false,
		callback = orbit_terminateE )

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

	dq1_prev = solE[2,1] + solE[3,1]
	dq2_prev = solE[4,1] - solE[1,1]
	for u in solE.u[2:end]
	    dq1 = u[2] + u[3]
	    dq2 = u[4] - u[1]

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


    QI =sqrt.(q1I.^2 + q2I.^2)
    QE =sqrt.(q1E.^2 + q2E.^2)

    PI =sqrt.(p1I.^2 + p2I.^2)
    PE =sqrt.(p1E.^2 + p2E.^2)

    plot([p1I,p1E],[p2I,p2E], seriestype=:scatter, xlims=(-1.5,1.5),ylims=(-1,3))
end
