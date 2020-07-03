include("types.jl")
include("utils.jl")
include("convexify.jl")

using Convex
using ECOS
using Printf
using LinearAlgebra

using .MyUtils

function scvx_initialize(bnds::ScvxBnds,
						 pars::ScvxParameters{<:ModelParameters},
						 tr::Float64)::ScvxProblem
	 print("Initializing")
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
	print(".")
	# create initial solution struct
	init_sol = ScvxSolution(x,u,tf,nx,pars.nu,N)
	print(".")
	# convexify along initial guess
	convexify!(init_sol,pars)
	print(".")

	# compute cost along initial guess
	J = opt_cost(init_sol.state,init_sol.control,init_sol.tf,N)
	J += pars.wvc * sum(init_sol.defects)

	println("done.")
	return ScvxProblem(bnds,pars,init_sol,J,tr)
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

function scvx_set_scales(prob::ScvxProblem)::ScvxScale
	# sizes
	nx = prob.pars.nx
	nu = prob.pars.nu

	# min/max state, control and final time
	x_min = prob.bnds.path.x_min
	x_max = prob.bnds.path.x_max
	u_min = prob.bnds.path.u_min
	u_max = prob.bnds.path.u_max
	t_min = prob.bnds.init.t_min
	t_max = prob.bnds.trgt.t_max

	# choose scaled variable intervals
	intrvl_x  = [0;1]
	wdth_x    = intrvl_x[2]-intrvl_x[1]
	intrvl_u  = [0;1]
	wdth_u    = intrvl_u[2]-intrvl_u[1]
	intrvl_t  = [0;1]
	wdth_t    = intrvl_t[2]-intrvl_t[1]

	# state terms
	Sx  = zeros(nx,nx)
	iSx = zeros(nx,nx)
	cx  = zeros(nx)
	for k = 1:nx
	    Sx[k,k]  = (x_max[k]-x_min[k])/wdth_x
	    iSx[k,k] = 1.0/Sx[k,k]
	    cx[k]    = x_min[k]-Sx[k,k]*intrvl_x[1]
	end

	# control terms
	Su = zeros(nu,nu)
	iSu = zeros(nu,nu)
	cu = zeros(nu)
	for k = 1:nu
		Su[k,k]  = (u_max[k]-u_min[k])/wdth_u
		iSu[k,k] = 1.0/Su[k,k]
		cu[k] 	 = u_min[k]-Sx[k,k]*intrvl_u[1]
	end

	# temporal terms
	Sp 	= (t_max-t_min)/wdth_t
	iSp = 1.0/Sp
	cp 	= t_min - Sp * intrvl_t[1]

	# create and return ScvxScale
	return ScvxScale(Sx,cx,Su,cu,Sp,cp,iSx,iSu,iSp)
end

function scvx_solve!(prob::ScvxProblem)::Integer
	# print intro
	println("Solving...")
	cvrg = false

	# set scaling matrices
	scale = scvx_set_scales(prob);

	for iter = 1:prob.pars.iter_max

	    # solve an SOCP
	    varxu = solve_socp!(prob,scale)

	    # convexify along new solution
	    convexify!(prob.new_sol,prob.pars)

	    # perform update step
		reject,change = scvx_update!(prob)

	    # check convergence and exit if done
		cvrg = scvx_convergence(varxu,prob.pars.cvrg_tol)

	    # print iterate information
	    scvx_print_status(prob,iter,varxu,reject,change);

		if cvrg
			break
		end
	end
	# return exitcode
	return scvx_exitcode(prob,cvrg)
end

function solve_socp!(prob::ScvxProblem,scale::ScvxScale)::Float64
	# get problem data
	nx  = prob.pars.nx
	nu  = prob.pars.nu
	N   = prob.pars.N
	wvc = prob.pars.wvc
	tr  = prob.tr

	# scaling matrices
	Sx  = scale.Sx
	cx  = scale.cx
	Su  = scale.Su
	cu  = scale.cu
	Sp  = scale.Sp
	cp  = scale.cp
	iSx = scale.iSx
	iSu = scale.iSu
	iSp = scale.iSp

	# reference solution
	x_ref = prob.prv_sol.state
	u_ref = prob.prv_sol.control
	p_ref = prob.prv_sol.tf
	Ad  = prob.prv_sol.Ad
	Bdm = prob.prv_sol.Bdm
	Bdp = prob.prv_sol.Bdp
	Sd  = prob.prv_sol.Sd
	Rd  = prob.prv_sol.Rd

	# variables
	xb  = Variable(nx,N)
	ub  = Variable(nu,N)
	pb  = Variable(1)
	vc  = Variable(nx,N-1)

	# cost function
	socp = minimize( opt_cost(xb,ub,pb,N) + wvc * norm_1(vec(vc)) )

	# constraints
	socp.constraints += (Sx*xb[:,1]+cx) <= prob.bnds.init.x_max
	socp.constraints += (Sx*xb[:,1]+cx) >= prob.bnds.init.x_min

	socp.constraints += (Sx*xb[:,N]+cx) <= prob.bnds.trgt.x_max
	socp.constraints += (Sx*xb[:,N]+cx) >= prob.bnds.trgt.x_min

	socp.constraints += (Sp*pb+cp) <= prob.bnds.trgt.t_max
	socp.constraints += (Sp*pb+cp) >= prob.bnds.trgt.t_min

	for k = 1:N
		xbk = xb[:,k]
		ubk = ub[:,k]
		xk  = (Sx*xbk+cx)
		uk  = (Su*ubk+cu)
		dxk = xbk - iSx * (x_ref[:,k]-cx)
		duk = ubk - iSu * (u_ref[:,k]-cu)
		# dynamics
		if k<N
			xbkp = xb[:,k+1]
			ubkp = ub[:,k+1]
			socp.constraints += (Sx*xbkp+cx) == Ad[:,:,k]*xk + Bdm[:,:,k]*uk + Bdp[:,:,k]*(Su*ubkp+cu) + Sd[:,k]*(Sp*pb+cp) + Rd[:,k] + vc[:,k]
		end

		# path constraints
		socp.constraints += xk <= prob.bnds.path.x_max
		socp.constraints += xk >= prob.bnds.path.x_min
		socp.constraints += uk <= prob.bnds.path.u_max
		socp.constraints += uk >= prob.bnds.path.u_min

		# trust region
		socp.constraints += dot(dxk,dxk) + dot(duk,duk) <= tr
	end

	# solve the problem
	solver = () -> ECOS.Optimizer(verbose=0)
	solve!(socp,solver)

	# save solution data
	prob.new_sol.state = Sx*(xb.value).+cx
	prob.new_sol.control = Su*(ub.value).+cu
	prob.new_sol.tf = Sp*(pb.value)+cp
	prob.new_sol.vc = norm(vec(vc.value),1)

	# compute max scaled change in state/control
	varxu = 0.0
	for k = 1:N
   		tempx = norm((xb.value)[:,k] - iSx*(x_ref[:,k]-cx),Inf);
   		tempu = norm((ub.value)[:,k] - iSu*(u_ref[:,k]-cu),Inf);
   		varxu = maximum([varxu;tempx;tempu]);
	end

	prob.new_sol.flag = Int(socp.status)
	return varxu
end # solve_socp

function scvx_update!(prob::ScvxProblem)
	cost = opt_cost(prob.new_sol.state,
					prob.new_sol.control,
					prob.new_sol.tf,
					prob.pars.N)
	new_L = cost + prob.pars.wvc * prob.new_sol.vc
	new_J = cost + prob.pars.wvc * sum(prob.new_sol.defects)

	dL = prob.prv_J - new_L
	dJ = prob.prv_J - new_J
	tr = prob.tr
	if dL>1e-12
		# compute performance metric
		ρ = dJ/dL
		#  update trust region radius
		if ( ρ < prob.pars.ρ_1 )
    		reject = true
    		tr = tr / prob.pars.α
    		change = 'S'
		else
    		reject = false;
    		prob.prv_J = new_J;
    		if ( ρ < prob.pars.ρ_1 )
        		# shrink
        		tr = tr / prob.pars.α
        		change = 'S';
    		elseif ( (ρ >= prob.pars.ρ_1) && (ρ < prob.pars.ρ_2) )
				# keep
        		change = 'K';
    		else
        		# grow
        		tr = tr * prob.pars.β;
        		change = 'G';
    		end #if
		end #if
	end #if

	# saturate trust region to lower/upper bound
	if (tr < prob.pars.tr_lb)
	    tr = prob.pars.tr_lb;
	    change = 'L';
	elseif (tr > prob.pars.tr_ub)
	    tr = prob.pars.tr_ub;
	    change = 'U';
	end

	# replace updated trust region
	prob.tr = tr;

	# if the step was not rejected, update the solution
	if !reject
		set_prv_solution!(prob,new_J)
	end # if
	return reject, change
end #scvx_update!

function set_prv_solution!(prob::ScvxProblem,new_J::Float64)
	prob.prv_sol = prob.new_sol
	prob.prv_J = new_J
	return nothing
end

function scvx_print_status(prob::ScvxProblem,iter::Integer,diff::Float64,
							reject::Bool,change::Char)
	vc = prob.new_sol.vc
	tr = prob.tr
	feas = prob.new_sol.feas
	@printf "Iter: %02d | " iter
	@printf "VC: %2.2e | " vc
	@printf "TR: %2.2e | " tr
	@printf "J: %2.2e | " prob.prv_J
	@printf "max diff: %2.2e | " diff
	@printf "feas: %d | " feas
	@printf "update: %d" reject
	@printf "%c\n" change
	return nothing
end

function scvx_convergence(varxu::Float64,tol::Float64)
	if varxu<tol
		return true
	else
		return false
	end # if
end

function scvx_exitcode(prob::ScvxProblem,converged::Bool)::Integer
	# Possible exitcodes are:
	#	0 : converged and feasible
	# 	1 : reached max iterations and IS feasible
	# 	2 : converged and IS NOT feasible
	# 	3 : reached max iterations and IS NOT feasible

	if prob.new_sol.feas
		if converged
			return 0
		else
			return 1
		end
	else
		if converged
			return 2
		else
			return 3
		end
	end
end

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
	return nothing
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
	return nothing
end
