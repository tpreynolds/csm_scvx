module Scvx

include("utils.jl")

using Convex
using ECOS
using Printf
using LinearAlgebra
import Base.+

using .Utils

export
	ScvxParameters,
	ScvxInitBnds,
	ScvxTrgtBnds,
	ScvxPathBnds,
	ScvxBnds,
	ScvxSolution,
	initialize,
	disp_sol,
	PointMassParameters,
	disp_Ad

################################################################################
################################################################################
################################################################################

abstract type Model end
abstract type ModelParameters end

struct PointMassParameters<:ModelParameters
	m::Real
	cd::Real
	Sd::Real
	ρ::Real
	id_r::Array{Integer,1}
	id_v::Array{Integer,1}
end
function PointMassParameters(m::Real,cd::Real,Sd::Real,ρ::Real)
	id_r = 1:3
	id_v = 4:6
	PointMassParameters(m,cd,Sd,ρ,id_r,id_v)
end

# mutable struct PointMass{T<:ModelParameters}<:Model
# 	nx::Integer
# 	nu::Integer
# 	pars::T
# end

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

################################################################################
################################################################################
################################################################################

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
# function ScvxParameters()

# end

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

function initialize(bnds::ScvxBnds,pars::ScvxParameters)::ScvxProblem

	# constants
	nx = pars.nx
	N = pars.N

	# initial state guess using straight line interpolation
	x0_min = bnds.init.x_min
	x0_max = bnds.init.x_max
	xf_min = bnds.trgt.x_min
	xf_max = bnds.trgt.x_max

	x0 = zeros(nx,1)
	xf = zeros(nx,1)
	for k = 1:nx
		x0[k] = 0.5*(x0_min[k] + x0_max[k])
	    xf[k] = 0.5*(xf_min[k] + xf_max[k])
	end
	x = init_straightline(x0,xf,N)

	# initial control guess using straight line interpolation
	u = init_straightline(bnds.path.u_min,bnds.path.u_min,N)

	# initial time guess halfway between bounds
	tf = 0.5 * (bnds.trgt.t_min + bnds.trgt.t_max)

	# create initial solution struct
	init_sol = ScvxSolution(x,u,tf,nx,pars.nu,N)

	# convexify along initial guess
	convexify!(init_sol,pars)

	return ScvxProblem(bnds,pars,init_sol)
end

function init_straightline(v0,vf,N)::Array{Float64,2}
	nv = size(v0,1)
	if (size(vf,1)!=nv)
		error("Inputs v0 and vf have incompatible sizes")
	end
	v = zeros(nv,N)
	for k = 1:nv
	    v[k,:] = range(v0[k],stop=vf[k],length=N)
	end
	return v
end

# function solve()

# end

# function solve_socp()

# solver = () -> ECOS.Optimizer(verbose=0);

# x = Variable(4);
# c = [1;2;3;4];
# A = I(4);
# b = [10;10;10;10];
# p = minimize(dot(c,x))
# p.constraints += A * x <= b;
# p.constraints += [ x>=1; x<=10; x[2]<=5; x[1]+x[4]-x[2]<=10];
# solve!(p,solver);

# println(p.optval);
# println(x.value);
# println(evaluate(x[1]+x[4]-x[2]));

# end

# function check_convergence()

# end

function disp_sol(sol::ScvxSolution,pars::ScvxParameters)
	t = range(0.,stop=sol.tf,length=pars.N)
	println("        Scvx solution:       ")
	println("=============================")
	println("Time   :   State   :  Control")
	for k = 1:pars.N
	    @printf "%05.2f  :  " t[k]
	    for i = 1:pars.nx
	        @printf "%+2.2e " sol.state[i,k]
	    end
	    @printf "  :  "
	    for i = 1:pars.nu
	    	@printf "%+2.2e " sol.control[i,k]
	    end
	    @printf "\n"
	end
end

function disp_Ad(sol::ScvxSolution,pars::ScvxParameters,rng::Integer=0)
	nx = pars.nx
	if rng==0
	    rng = 1:pars.N-1
	end
	println("Discrete A Matrices")
	println("===================")
	for k ∈ rng
		@printf "Interval %02d\n" k
		for i = 1:nx
			for j=1:nx
				if j < nx
					@printf "%5.4f, " sol.Ad[i,j,k]
				else
					@printf "%5.4f\n" sol.Ad[i,j,k]
				end
			end
		end
		print("\n")
	end
end

################################################################################
################################################################################
################################################################################


function +(i::Int64,rng::UnitRange{Int64})
	return UnitRange{Int64}(rng.start+i,rng.stop+i)
end

mutable struct cvxfy_params
	nx::Integer
	id_x::UnitRange{Integer}
	id_A::UnitRange{Integer}
	id_Bm::UnitRange{Integer}
	id_Bp::UnitRange{Integer}
	id_S::UnitRange{Integer}
	id_R::UnitRange{Integer}
	u::Array{Float64,2}
	tf::Float64
	tauspan::Array{Float64,1}
	mdl_pars
end
function cvxfy_params(nx::Integer,nu::Integer,tf::Float64,pars::T) where T<:ModelParameters
	id_x  = (1:nx)
	id_A  = id_x[end] + (1:nx*nx)
	id_Bm = id_A[end] + (1:nx*nu)
	id_Bp = id_Bm[end] + (1:nx*nu)
	id_S  = id_Bp[end] + (1:nx)
	id_R  = id_S[end] + (1:nx)

	return cvxfy_params(nx,id_x,id_A,id_Bm,id_Bp,id_S,id_R,
					Array{Float64}(undef,nu,2),
					tf,
					Array{Float64}(undef,2),
					pars)
end

function convexify!(sol::ScvxSolution,pars::ScvxParameters)
	nx = pars.nx
	nu = pars.nu
	N = pars.N
	Nsub = pars.Nsub
	tau = LinRange(0.0,1.0,N)

	# initialize data
	cvxfy_pars = cvxfy_params(nx,nu,sol.tf,pars.mdl_pars)
	feas = true

	V0 = zeros(nx+nx*nx+2*nx*nu+2*nx)
	V0[cvxfy_pars.id_A] = vec(I(nx))

	for k = 1:N-1
		# set initial condition
		# V = V0
		V0[cvxfy_pars.id_x] = sol.state[:,k]
		# set remaining parameters
		tauspan = LinRange(tau[k],tau[k+1],Nsub)
		cvxfy_pars.tauspan = [ tauspan[1]; tauspan[end] ]
		cvxfy_pars.u = sol.control[:,k:k+1]
		# integrate
		f(t,v) = deriv(t,v,cvxfy_pars)
		V = rk4(V0,f,tauspan)
		# get results
		xV  = V[cvxfy_pars.id_x]
		AV  = V[cvxfy_pars.id_A]
		BmV = V[cvxfy_pars.id_Bm]
		BpV = V[cvxfy_pars.id_Bp]
		SV  = V[cvxfy_pars.id_S]
		RV  = V[cvxfy_pars.id_R]
		# reshape and get final results
		Adk  = reshape(AV,(nx,nx))
		Bdmk = Adk * reshape(BmV,(nx,nu))
		Bdpk = Adk * reshape(BpV,(nx,nu))
		Sdk  = Adk * SV
		Rdk  = Adk * RV
		# fill up matrices
		sol.Ad[:,:,k]  = Adk
		sol.Bdm[:,:,k] = Bdmk
		sol.Bdp[:,:,k] = Bdpk
		sol.Sd[:,k]    = Sdk
		sol.Rd[:,k]    = Rdk
		# compute defect
		xkp = sol.state[:,k+1]
		sol.defects[k] = norm(xkp-xV)
		if sol.defects[k]>pars.feas_tol&&feas
		    feas = false
		end
	end
	sol.feas = feas
	return sol
end

function deriv(tau,V,cvxfy_pars)
	nx = cvxfy_pars.nx
	tauspan = cvxfy_pars.tauspan
	x_tau = V[1:nx]
	u_tau = interp_vec(tau,cvxfy_pars.u,tauspan)
	tf = cvxfy_pars.tf
	PHI = reshape(V[nx+(1:nx*nx)],(nx,nx))
	lm = (tauspan[2]-tau)/(tauspan[2]-tauspan[1])
	lp = (tau-tauspan[1])/(tauspan[2]-tauspan[1])

	f = dynamics(tau,x_tau,u_tau,cvxfy_pars.mdl_pars)
	A,B = linearize(tau,x_tau,u_tau,cvxfy_pars.mdl_pars)

	Bm = lm * B
	Bp = lp * B
	R  = - A * x_tau - B * u_tau
	iPHI = PHI\I(nx)

	dx = tf * f
	dPHI = (tf * A) * PHI
	dBm = iPHI * (tf * Bm)
	dBp = iPHI * (tf * Bp)
	dS = iPHI * f
	dR = iPHI * (tf * R)

	dV = [ dx; vec(dPHI); vec(dBm); vec(dBp); vec(dS); vec(dR) ]
end

################################################################################
################################################################################
################################################################################

end # module
