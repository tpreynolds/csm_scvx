# module Scvx

include("types.jl")
include("utils.jl")
include("convexify.jl")

using Convex
using ECOS
using Printf
using LinearAlgebra

using .Utils
# using .Convexify

# export
# 	initialize,
# 	disp_sol,
# 	disp_Ad

################################################################################
################################################################################
################################################################################


################################################################################
################################################################################
################################################################################

function initialize(bnds::ScvxBnds,pars::ScvxParameters{<:ModelParameters})::ScvxProblem

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

# end # module
