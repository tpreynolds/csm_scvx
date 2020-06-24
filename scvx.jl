module Scvx

using Convex
using ECOS

export
	ScvxParameters,
	ScvxInitBnds,
	ScvxTrgtBnds,
	ScvxPathBnds,
	ScvxBnds

struct ScvxParameters
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

struct ScvxSolution
	state::Array{Float64,2}
	control::Array{Float64,2}
end

# function initialize()

# end

# function solve()

# end

# function convexify()

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



end # module