include("src/scvx.jl")
include("models/quadrotor_mdl.jl")
include("utils/csm_plots.jl")

# parameters
id_r = 1:3
id_v = 4:6
id_w = 1:3
id_G = 4
u_nrm_max = 23.2
u_nrm_min = 0.6
tilt_max = deg2rad(60)
obsN = 2
obsiH1 = diagm([2;2;0])
obsiH2 = diagm([2;2;0])
obsc1 = [1.0;2.0;0.0]
obsc2 = [2.0;5.0;0.0]
obsiH = cat(obsiH1,obsiH2;dims=3)
obsC  = cat(obsc1,obsc2;dims=2)

pars = QuadrotorParameters(id_r,id_v,id_w,id_G,u_nrm_max,u_nrm_min,tilt_max,obsN,obsiH,obsC)

# initial conditions
t0_min = 0.0
t0_max = 0.0
x0_min = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ]
x0_max = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ]
u0_min = [ 0.0; 0.0; 9.81; NaN ]
u0_max = [ 0.0; 0.0; 9.81; NaN ]

# final conditions
tf_min = 2.5
tf_max = 2.5
xf_min = [ 2.5; 6.0; 0.0; 0.0; 0.0; 0.0 ]
xf_max = [ 2.5; 6.0; 0.0; 0.0; 0.0; 0.0 ]
uf_min = [ 0.0; 0.0; 9.81; NaN ]
uf_max = [ 0.0; 0.0; 9.81; NaN ]

# linear path constraints
x_min = [ -10.0; -10.0; -10.0; -8.0; -8.0; -8.0 ]
x_max = [ 10.0; 10.0; 10.0; 8.0; 8.0; 8.0 ]
u_min = [ -u_nrm_max; -u_nrm_max; 0; 0 ]
u_max = [ u_nrm_max; u_nrm_max; u_nrm_max; u_nrm_max ]

# sizes and tolerances
N = 30
Nsub = 10
nx = 6
nu = 4
iter_max = 15
wvc = 1e2
ρ_0 = 0.0
ρ_1 = 0.1
ρ_2 = 0.7
α   = 2.0
β   = 2.0
tr  = 1.0
tr_lb = 0.001
tr_ub = 10.
cvrg_tol = 1e-3
feas_tol = 1e-2

init = ScvxInitBnds(t0_min,t0_max,x0_min,x0_max,u0_min,u0_max)
trgt = ScvxTrgtBnds(tf_min,tf_max,xf_min,xf_max,uf_min,uf_max)
path = ScvxPathBnds(x_min,x_max,u_min,u_max)
bnds = ScvxBnds(init,trgt,path)
ctrl = ScvxParameters(N,Nsub,nx,nu,iter_max,wvc,
                        ρ_0,ρ_1,ρ_2,α,β,tr_lb,tr_ub,
                        cvrg_tol,feas_tol,pars)

# initialize an scvx problem
prob = scvx_initialize(bnds,ctrl,tr)

# solve the scvx problem
flag = scvx_solve!(prob)

# plot the results
csm_plots_quad(prob)
