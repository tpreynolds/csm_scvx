struct PointMassParameters<:ModelParameters
	m::Real
	cd::Real
	Sd::Real
	ρ::Real
	id_r::Array{Integer,1}
	id_v::Array{Integer,1}
end

function dynamics(t::Float64,x,u,t_grid,pars::T) where {T<:ModelParameters}
	r = x[pars.id_r]
	v = x[pars.id_v]

	ut = interp_vec(t,u,t_grid)
	fD, = drag(v,pars)
	g, = gravity(r)

	dr = v
	dv = (1/pars.m)*(ut+fD) + g
	dx = [ dr; dv ]

	return dx
end # dynamics for integration

function dynamics(t::Float64,x,u,pars::T) where {T<:ModelParameters}
	r = x[pars.id_r]
	v = x[pars.id_v]

	fD, = drag(v,pars)
	g, = gravity(r)

	dr = v;
	dv = (1/pars.m)*(u+fD) + g
	dx = [ dr; dv ]

	return dx
end # dynamics for evaluation

function linearize(t::Float64,x,u,pars::T) where {T<:ModelParameters}
	r = x[pars.id_r]
	v = x[pars.id_v]

	_,dg_dr = gravity(r)
	_,dfD_dv = drag(v,pars)

	A = zeros(6,6)
	A[pars.id_r,pars.id_v] = I(3)
	A[pars.id_v,pars.id_r] = dg_dr
	A[pars.id_v,pars.id_v] = dfD_dv

	B = zeros(6,3)
	B[pars.id_v,:] = (1/pars.m)*I(3)

	return A, B
end

function gravity(r)
	g = [0.0;0.0;-1.0]
	dg_dr = zeros(3,3)
	return g, dg_dr
end

function drag(v,pars::T) where {T<:ModelParameters}
	q = -0.5 * pars.ρ * pars.Sd * pars.cd;
	speed = norm(v);

	fD = q * speed .* v;

	dfD_dv = zeros(3,3)
	if speed>1e-12
	    dfD_dv = q * (speed * I(3) + (v*transpose(v))/speed )
	end

	return fD, dfD_dv
end
