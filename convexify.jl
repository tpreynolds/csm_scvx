module Convexify

include("types.jl")
include("model.jl")

export 
	convexify!

function convexify!(sol::ScvxSolution,pars::ScvxParameters)
	x = sol.state[:,1]
	u = sol.control[:,1]
	t = 0.

	model_pars = pars.model_pars

	dx = dynamics(t,x,u,model_pars)
	println(dx)
end

end # module