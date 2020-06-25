# include("types.jl")

struct ModelParameters
	m::Real
	cd::Real
	Sd::Real
	ρ::Real
	id_r::Array{Integer,1}
	id_v::Array{Integer,1}
end
function ModelParameters(m::Real,cd::Real,Sd::Real,ρ::Real)
	id_r = 1:3
	id_v = 4:6
	ModelParameters(m,cd,Sd,ρ,id_r,id_v)
end

mutable struct Model
	nx::Integer
	nu::Integer
	pars::ModelParameters
end

function dynamics(t::Float64,x,u,t_grid,pars::ModelParameters) 
	r = x[pars.id_r]
	v = x[pars.id_v]

	ut = interp_vec(t,u,t_grid);
	fD = drag(v,pars);

	dr = v;
	dv = (1/pars.m)*(ut+fD) + gravity(r);
	dx = [ dr; dv ];

	return dx
end # dynamics for integration

function dynamics(t::Float64,x,u,pars::ModelParameters) 
	r = x[pars.id_r]
	v = x[pars.id_v]

	fD = drag(v,pars);

	dr = v;
	dv = (1/pars.m)*(u+fD) + gravity(r);
	dx = [ dr; dv ];

	return dx
end # dynamics for evaluation

function gravity(r)
	g = [0.0;0.0;-1.0];
end

function drag(v,pars::ModelParameters)
	q = -0.5 * pars.ρ * pars.Sd * pars.cd;
	speed = norm(v);

	fD = q * speed .* v;
end
