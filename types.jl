abstract type Model end
abstract type ModelParameters end

struct ScvxParameters{T<:ModelParameters}
	N::Integer
	Nsub::Integer
	nx::Integer
	nu::Integer
	iter_max::Integer
	wvc::Float64
	wtr::Float64
	wtrp::Float64
	cvrg_tol::Float64
	feas_tol::Float64
	mdl_pars::T
end

struct ScvxInitBnds
	t_min::Real
	t_max::Real
	x_min
	x_max
	u_min
	u_max
end

struct ScvxTrgtBnds
	t_min::Real
	t_max::Real
	x_min
	x_max
	u_min
	u_max
end

struct ScvxPathBnds
	x_min
	x_max
	u_min
	u_max
end

struct ScvxBnds
	init::ScvxInitBnds
	trgt::ScvxTrgtBnds
	path::ScvxPathBnds
end

mutable struct ScvxSolution
	state::Array{Float64,2}
	control::Array{Float64,2}
	tf::Float64

	Ad::Array{Float64,3}
	Bdm::Array{Float64,3}
	Bdp::Array{Float64,3}
	Sd::Array{Float64,2}
	Rd::Array{Float64,2}
	defects::Array{Float64,1}

	feas::Bool
end
function ScvxSolution(x::Array{Float64,2},u::Array{Float64,2},t::Float64,nx::Integer,nu::Integer,N::Integer)
	ScvxSolution(x,u,t,
		Array{Float64}(undef,(nx,nx,N-1)),
		Array{Float64}(undef,(nx,nu,N-1)),
		Array{Float64}(undef,(nx,nu,N-1)),
		Array{Float64}(undef,(nx,N-1)),
		Array{Float64}(undef,(nx,N-1)),
		Array{Float64}(undef,N-1),
		false)
end

mutable struct ScvxProblem
	bnds::ScvxBnds
	pars::ScvxParameters
	prv_sol::ScvxSolution
	new_sol::ScvxSolution
	prv_J::Float64
end
function ScvxProblem(bnds::ScvxBnds,pars::ScvxParameters,sol::ScvxSolution)
	ScvxProblem(bnds,pars,sol,sol,0.0)
end
