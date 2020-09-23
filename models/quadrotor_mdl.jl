struct QuadrotorParameters<:ModelParameters
	id_r::Array{Integer,1}
	id_v::Array{Integer,1}
	id_w::Array{Integer,1}
	id_G::Integer
	u_nrm_max::Float64
	u_nrm_min::Float64
	tilt_max::Float64
	obsN::Integer
	obsiH::Array{Float64,3}
	obsC::Array{Float64,2}
end

function dynamics(t::Float64,x,u,t_grid,pars::T) where {T<:ModelParameters}
	# interpolate controls
	ut = interp_vec(t,u,t_grid)
	# return using regular dynamics function
	return dynamics(t,x,ut,pars)
end # dynamics for integration

function dynamics(t::Float64,x,u,pars::T) where {T<:ModelParameters}
	r = x[pars.id_r]
	v = x[pars.id_v]
	w = u[pars.id_w]

	g, = gravity(r)

	dr = v
	dv = w + g
	dx = [ dr; dv ]

	return dx
end # dynamics for evaluation

function gravity(r)
	g = [0.0;0.0;-9.81]
	dg_dr = zeros(3,3)
	return g, dg_dr
end

function linearize(t::Float64,x,u,pars::T) where {T<:ModelParameters}
	r = x[pars.id_r]
	v = x[pars.id_v]
	w = u[pars.id_w]

	_,dg_dr = gravity(r)

	A = zeros(6,6)
	A[pars.id_r,pars.id_v] = I(3)
	A[pars.id_v,pars.id_r] = dg_dr

	B = zeros(6,4)
	B[pars.id_v,pars.id_w] = I(3)

	return A, B
end

function init_solution(bnds::ScvxBnds,N::Integer)

	# initial state guess using straight line interpolation
	x0_min = bnds.init.x_min
	x0_max = bnds.init.x_max
	xf_min = bnds.trgt.x_min
	xf_max = bnds.trgt.x_max
	x0 = zeros(6,1)
	xf = zeros(6,1)
	for k = 1:6
		x0[k] = 0.5*(x0_min[k] + x0_max[k])
	    xf[k] = 0.5*(xf_min[k] + xf_max[k])
	end
	# initial state guess using straight line interpolation
	x = init_straightline(x0,xf,N)
	# initial control guess using straight line interpolation
	u = init_straightline(bnds.path.u_min,bnds.path.u_min,N)
	g, = gravity([0.0;0.0;1.0]) # input here is arbitrary
	u[:,1] = [-g; norm(g)]
	u[:,N] = [-g; norm(g)]
	# initial time guess halfway between bounds
	p = 0.5 * (bnds.trgt.t_min + bnds.trgt.t_max)
	return x,u,p
end

function opt_cost(x,u,t,N::Integer)
	J = 0.0;
	id_G = 4
	for k = 1:N-1
	    Gk  = u[id_G,k];
	    Gkp = u[id_G,k+1];
	    J += 0.5 * ( Gk*Gk + Gkp*Gkp );
	end
	return 1e0*J
end

function mdl_cvx_constraints!(socp,xk,uk,pars::T) where T<:ModelParameters
	id_w = pars.id_w
	id_G = pars.id_G
	u_nrm_max = pars.u_nrm_max
	u_nrm_min = pars.u_nrm_min
	tilt_max  = pars.tilt_max

	wk = uk[id_w]
	Gk = uk[id_G]
	# add relaxed norm constraints
	socp.constraints += u_nrm_min - Gk <= 0.0
	socp.constraints += Gk - u_nrm_max <= 0.0
	socp.constraints += norm(wk) - Gk <= 0.0
	# add relaxed thrust pointing constraint
	socp.constraints += Gk*cos(tilt_max) - wk[3] <= 0.0
	return nothing
end

function mdl_ncvx_constraints!(socp,xk,uk,xbk,ubk,vbk,pars::T) where T<:ModelParameters
	# loop through nonconvex constraints and add approximations
	for i = 1:pars.obsN
		_,A,b = obstacle_constraint(xbk,pars,i)
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
	A = zeros(6)
	if norm(temp)>eps()
		A[id_r] = -transpose(rk-c)*(transpose(iH)*iH)/norm(temp)
	end
	b = f - dot(A,xk)

	return f,A,b
end
