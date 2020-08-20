include("ISS_mdl.jl")

struct FreeFlyerParameters<:ModelParameters
	id_r::Array{Integer,1}
	id_v::Array{Integer,1}
	id_q::Array{Integer,1}
	id_w::Array{Integer,1}
	id_F::Array{Integer,1}
	id_M::Array{Integer,1}
	F_nrm_max::Float64
	M_nrm_max::Float64
	r_nrm_max::Float64
	v_nrm_max::Float64
	w_nrm_max::Float64
	mass::Float64
	inertia::Array{Float64,2}
	obsN::Integer
	obsiH::Array{Float64,3}
	obsC::Array{Float64,2}
	kozN::Integer
	koz::Array{PolygonalObstacle,1}
end

function dynamics(t::Float64,x,u,t_grid,pars::T) where {T<:ModelParameters}
	# interpolate controls
	ut = interp_vec(t,u,t_grid)
	# return using regular dynamics function
	return dynamics(t,x,ut,pars)
end # dynamics for integration

function dynamics(t::Float64,x,u,pars::T) where {T<:ModelParameters}
	v = x[pars.id_v]
	q = x[pars.id_q]
	w = x[pars.id_w]
	F = u[pars.id_F]
	M = u[pars.id_M]

	m = pars.mass
	J = pars.inertia

	dr = v
	dv = F/m
	dq = 0.5 * quat_mult(q,w)
	dw = J\( M - skew(w) * J * w )
	dx = [ dr; dv; dq; dw ]

	return dx
end # dynamics for evaluation

function linearize(t::Float64,x,u,pars::T) where {T<:ModelParameters}
	v = x[pars.id_v]
	q = x[pars.id_q]
	w = x[pars.id_w]

	m = pars.mass
	J = pars.inertia

	dfq_dq = 0.5 * quat_skew_star(w)
	dfq_dw = 0.5 * quat_skew(q)
	dfw_dw = -J\( skew(w)*J - skew(J*w) )

	A = zeros(13,13)
	A[pars.id_r,pars.id_v] = I(3)
	A[pars.id_q,pars.id_q] = dfq_dq
	A[pars.id_q,pars.id_w] = dfq_dw[:,1:3]
	A[pars.id_w,pars.id_w] = dfw_dw

	B = zeros(13,6)
	B[pars.id_v,pars.id_F] = (1.0/m)*I(3)
	B[pars.id_w,pars.id_M] = J\I(3)

	return A, B
end

function opt_cost(x,u,t,N::Integer)
	J = 0.0;
	id_F = 1:3
	id_M = 4:6
	for k = 1:N-1
	    Fk  = u[id_F,k]
	    Fkp = u[id_F,k+1]
		Mk  = u[id_M,k]
		Mkp = u[id_M,k+1]
	    J += 0.5 * ( dot(Fk,Fk) + dot(Fkp,Fkp) + dot(Mk,Mk) + dot(Mkp,Mkp) )
	end
	return 1e-1*J
end

function mdl_cvx_constraints!(socp,xk,uk,pars::T) where T<:ModelParameters
	F_nrm_max = pars.F_nrm_max
	M_nrm_max = pars.M_nrm_max
	r_nrm_max = pars.r_nrm_max
	v_nrm_max = pars.v_nrm_max
	w_nrm_max = pars.w_nrm_max

	rk = xk[pars.id_r]
	vk = xk[pars.id_v]
	wk = xk[pars.id_w]
	Fk = uk[pars.id_F]
	Mk = uk[pars.id_M]
	# add norm constraints
	socp.constraints += norm(rk,Inf) - r_nrm_max <= 0.0
	socp.constraints += norm(vk) - v_nrm_max <= 0.0
	socp.constraints += norm(wk) - w_nrm_max <= 0.0
	socp.constraints += norm(Fk) - F_nrm_max <= 0.0
	socp.constraints += norm(Mk) - M_nrm_max <= 0.0
	return nothing
end

function mdl_ncvx_constraints!(socp,xk,uk,xbk,ubk,vbk,pars::T) where T<:ModelParameters
	# loop through nonconvex constraints and add approximations
	for i = 1:pars.obsN
		_,A,b = obstacle_constraint(xbk,pars,i)
		socp.constraints += dot(A,xk)+b - vbk <= 0.0
	end

	for i = 1:pars.kozN
		_,A,b = koz_constraint(xbk,pars,i)
		socp.constraints += dot(A,xk)+b - vbk <= 0.0
	end

	return nothing
end

function obstacle_constraint(xk,pars::T,id::Integer) where T<:ModelParameters
	id_r = pars.id_r
	rk 	 = xk[id_r]
	# obstacle center & dimensions
	iH = pars.obsiH[:,:,id]
	c  = pars.obsC[:,id]
	temp = iH*(rk-c)

	# compute constraint value s.t. f(x)<=0
	f = 1 - norm(temp)
	# compute constraint derivative s.t. f(x)<=0 approx A*x+b<=0
	A = zeros(13)
	if norm(temp)>eps()
		A[id_r] = -transpose(rk-c)*(transpose(iH)*iH)/norm(temp)
	end
	b = f - dot(A,xk)

	return f,A,b
end

function koz_constraint(xk,pars::T,id::Integer) where T<:ModelParameters
	id_r = pars.id_r
	rk 	 = xk[id_r]
	obs  = pars.koz[id]
	# compute signed distance & closest point to r on obstacle
	sd, close_pt = signed_distance(rk,obs)
	# compute constraint value s.t. f(x) <= 0
	f = -sd
	# compute constraint derivative s.t. f(x)<=0 approx A*x+b<=0
	A = zeros(13)
	# if abs(sd)>eps()
		A[id_r] = -(rk-close_pt)/sd
	# end
	b = f - dot(A,xk)

	return f,A,b
end
