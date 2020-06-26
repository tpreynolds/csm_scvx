include("scvx.jl")
include("model.jl")

# using .Scvx

println("Point Mass Example")

# parameters
m = 1.0
Cd = 0.5
Sd = 0.05
ρ = 1.25
id_r = 1:3
id_v = 4:6

pars = PointMassParameters(m,Cd,Sd,ρ,id_r,id_v)

# initial conditions
t0_min = 0.0
t0_max = 0.0
x0_min = [ 5.0; 6.0; 24.0; -4.0; -2.0; 0.0 ]
x0_max = [ 5.0; 6.0; 24.0; -4.0; -2.0; 0.0 ]
u0_min = [ NaN; NaN; NaN ]
u0_max = [ NaN; NaN; NaN ]

# final conditions
tf_min = 4.0
tf_max = 12.0
xf_min = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ]
xf_max = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ]
uf_min = [ NaN; NaN; NaN ]
uf_max = [ NaN; NaN; NaN ]

# path constraints
x_min = [ -10.0; -10.0; -1.0; -5.0; -5.0; -5.0 ]
x_max = [ 10.0;  10.0; 25.0;  5.0;  5.0;  5.0 ]
u_min = [ -2; -2; -2 ]
u_max = [ 2; 2; 2 ]

# sizes and tolerances
N = 20
Nsub = 10
nx = 6
nu = 3
iter_max = 20
wvc = 1e2
wtr = 1e-2
wtrp = 1e-3
cvrg_tol = 1e-2
feas_tol = 1e-2

init = ScvxInitBnds(t0_min,t0_max,x0_min,x0_max,u0_min,u0_max)
trgt = ScvxTrgtBnds(tf_min,tf_max,xf_min,xf_max,uf_min,uf_max)
path = ScvxPathBnds(x_min,x_max,u_min,u_max)
bnds = ScvxBnds(init,trgt,path)
ctrl = ScvxParameters(N,Nsub,nx,nu,iter_max,wvc,wtr,wtrp,cvrg_tol,feas_tol,pars)

prob = initialize(bnds,ctrl)

# output = solve(prob)
