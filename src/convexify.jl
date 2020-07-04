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
	return nothing
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
