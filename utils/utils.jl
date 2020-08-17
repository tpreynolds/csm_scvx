module MyUtils

import Base.+

export
	rk4,
	interp_vec,
	rad2deg_arr!,
	skew,
	quat_mult,
	quat_skew,
	quat_skew_star,
	quat_2_dcm,
	quat_2_rpy

function +(i::Int64,rng::UnitRange{Int64})
	return UnitRange{Int64}(rng.start+i,rng.stop+i)
end

function rk4(f::Function,tspan,x0::Array{Float64,1})

	N = length(tspan)
	n = length(x0)
	X = zeros(n,N)

	X[:,1] = x0

	for k = 1:N-1
	  tk = tspan[k]
  	  tkp1 = tspan[k+1]
	  h = tkp1-tk

  	  xk = X[:,k]
	  k1 = f(tk,xk)
	  k2 = f(tk+h/2,xk+h/2*k1)
	  k3 = f(tk+h/2,xk+h/2*k2)
	  k4 = f(tk+h,xk+h*k3)

	  X[:,k+1] = xk + h/6*(k1+2*k2+2*k3+k4)
	end

	return X
end

function rk4(x::Array{Float64,1},f::Function,tspan)
	N = length(tspan)
	n = length(x)

	for k = 1:N-1
	  tk = tspan[k]
  	  tkp1 = tspan[k+1]
	  h = tkp1-tk

	  k1 = f(tk,x)
	  k2 = f(tk+h/2,x+h/2*k1)
	  k3 = f(tk+h/2,x+h/2*k2)
	  k4 = f(tk+h,x+h*k3)

	  x += h/6*(k1+2*k2+2*k3+k4)
	end
	return x
end

function skew(v::Vector{Float64})
	return [0. -v[3] v[2]; v[3] 0. -v[1]; -v[2] v[1] 0. ]
end

function interp_vec(t::Float64, v, t_grid::Union{Array{Float64,1},LinRange{Float64}})
	if length(t_grid)!= size(v,2)
	    error("incompatible sizes between time grid and input vector")
	end
	k = get_interval(t,t_grid);
	linterp = (t_grid[k+1]-t)/(t_grid[k+1]-t_grid[k]);
	v_t = linterp * v[:,k] + (1-linterp) * v[:,k+1];
	return v_t
end

function get_interval(t::Float64,t_grid::Union{Array{Float64,1},LinRange{Float64}})
	if t<t_grid[1]
		error("time %f is less than the grid lower value of %f",t,t_grid[1])
	end
	if t>t_grid[end]
	    error("time %f is greater than the grid upper value of %f",t,t_grid[end])
	end
	k = 0
	while t>t_grid[k+1]
	    k += 1
	end
	if k==0
	    interval = 1
	else
		interval = k
	end

 	return interval
end

function rad2deg_arr!(v::Array{Float64,2})
	for ii = 1:size(v,1)
		for jj = 1:size(v,2)
			v[ii,jj] = rad2deg(v[ii,jj])
		end
	end
end

### QUATERNION UTILITIES
### Using Hamiltonian scalar last convention

function quat_mult(q::Vector{Float64},p::Vector{Float64})
	if length(q)==3
		q = [ q; 0 ]
	elseif length(q)!=4
		error("Quaternion input 1 must be of length 3 or 4")
	end
	if length(p)==3
		p = [ p; 0 ]
	elseif length(q)!=4
		error("Quaternion input 2 must be of length 3 or 4")
	end

	return [ q[4]*p[1]-q[3]*p[2]+q[2]*p[3]+q[1]*p[4];
			 q[3]*p[1]+q[4]*p[2]-q[1]*p[3]+q[2]*p[4];
			-q[2]*p[1]+q[1]*p[2]+q[4]*p[3]+q[3]*p[4];
			-q[1]*p[1]-q[2]*p[2]-q[3]*p[3]+q[4]*p[4] ]
end

function quat_skew(q::Vector{Float64})
	if length(q)==3
		q = [ q; 0 ]
	elseif length(q)!=4
		error("Quaternion input must be of length 3 or 4")
	end

	return [ q[4] -q[3]  q[2] q[1];
			 q[3]  q[4] -q[1] q[2];
			-q[2]  q[1]  q[4] q[3];
			-q[1] -q[2] -q[3] q[4] ]
end

function quat_skew_star(q::Vector{Float64})
	if length(q)==3
		q = [ q; 0 ]
	elseif length(q)!=4
		error("Quaternion input must be of length 3 or 4")
	end

	return [ q[4]  q[3] -q[2] q[1];
			-q[3]  q[4]  q[1] q[2];
			 q[2] -q[1]  q[4] q[3];
			-q[1] -q[2] -q[3] q[4] ]
end

function quat_2_dcm(q::Vector{Float64})
	if length(q)!=4
		error("Quaternion input must be of length 4")
	end
	q1 = q[1]
	q2 = q[2]
	q3 = q[3]
	q4 = q[4]

	return [ 2*q1^2+2*q4^2-1.0		2*(q1*q2+q3*q4) 	2*(q1*q3-q2*q4);
	 	  	 2*(q1*q2-q3*q4) 		2*q2^2+2*q4^2-1.0 	2*(q2*q3+q1*q4);
		  	 2*(q1*q3+q2*q4) 		2*(q2*q3-q1*q4) 	2*q3^2+2*q4^2-1.0 ]
end

function quat_2_rpy(q::Array{Float64,2})
	if size(q,1)!=4
		q = transpose(q)
	end
	if size(q,1)!=4
		error("Quaternion input must have four rows")
	end
	N = size(q,2)
	rpy = zeros(3,N)
	for k = 1:size(q,2)
		Rk = quat_2_dcm(q[:,k])

		rpy[1,k] =  atan( Rk[2,3], Rk[3,3] ) # roll
		rpy[2,k] = -asin( Rk[1,3] ) 		 # pitch
		rpy[3,k] =  atan( Rk[1,2], Rk[1,1] ) # yaw
	end
	return rpy
end

end # module
