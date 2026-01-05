struct Section
    q_solve ::Union{Function,Nothing}
    p_solve ::Union{Function,Nothing}
end

dq1_ltdq2 = Section(
    (q,mu,H0) -> begin
	q1,q2 = q
	p1 = -q2
	U = mu/sqrt((q1-1+mu)^2+q2^2) + (1-mu)/sqrt((q1+mu)^2+q2^2)
        c = -U - 0.5*q2^2 - H0
        p2 = q1 - sqrt(q1^2-2*c)
	return (p1,p2)
    end,
    nothing
)
dq1_gtdq2 = Section(
    (q,mu,H0) -> begin
	q1,q2 = q
	p1 = -q2
	U = mu/sqrt((q1-1+mu)^2+q2^2) + (1-mu)/sqrt((q1+mu)^2+q2^2)
        c = -U - 0.5*q2^2 - H0
        p2 = q1 + sqrt(q1^2-2*c)
	return (p1,p2)
    end,
    nothing
)
dq2_ltdq1 = Section(
    (q,mu,H0) -> begin
	q1,q2 = q
	p2 = q1
	U = mu/sqrt((q1-1+mu)^2+q2^2) + (1-mu)/sqrt((q1+mu)^2+q2^2)
        c = -U - 0.5*q1^2 - H0
        p1 = -q2 - sqrt(q2^2-2*c)
	return (p1,p2)
    end,
    nothing
)
dq2_gtdq1 = Section(
    (q,mu,H0) -> begin
	q1,q2 = q
	p2 = q1
	U = mu/sqrt((q1-1+mu)^2+q2^2) + (1-mu)/sqrt((q1+mu)^2+q2^2)
        c = -U - 0.5*q1^2 - H0
        p1 = -q2 + sqrt(q2^2-2*c)
	return (p1,p2)
    end,
    nothing
)
return Dict("dq1_ltdq2" => dq1_ltdq2,
	    "dq1_gtdq2" => dq1_gtdq2,
	    "dq2_ltdq1" => dq2_ltdq1,
	    "dq2_gtdq1" => dq2_gtdq1)

