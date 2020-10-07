include("src/scvx.jl")
include("models/freeflyer_mdl.jl")
include("utils/csm_plots.jl")

# parameters
id_r = 1:3
id_v = 4:6
id_q = 7:10
id_w = 11:13
id_F = 1:3
id_M = 4:6
r_nrm_max = 12
v_nrm_max = 0.4
w_nrm_max = 1
F_nrm_max = 72e-3
M_nrm_max = 2e-3
mass = 7.2
inertia = diagm([0.1083,0.1083,0.1083])
radius = sqrt(3) * (0.05/2)
obsN = 3
obsiH1 = diagm([1.0/0.3;1.0/0.3;1.0/0.3])
obsiH2 = obsiH1
obsiH3 = obsiH1
obsc1 = [11.3;3.8;4.8]
obsc3 = [11.2;1.84;5.0]
obsc2 = [8.5;-0.15;5.0]
obsiH = cat(obsiH1,obsiH2,obsiH3;dims=3) # obsiH2
obsC  = cat(obsc1,obsc2,obsc3;dims=2) # obsc1
kozN  = 24
koz   = ISS_koz(kozN)
# kizN  = 6
# kiz   = ISS_kiz(kizN)

pars = FreeFlyerParameters(id_r,id_v,id_q,id_w,id_F,id_M,F_nrm_max,M_nrm_max,r_nrm_max,v_nrm_max,w_nrm_max,mass,inertia,radius,obsN,obsiH,obsC,kozN,koz)

# initial conditions
t0_min = 0.0
t0_max = 0.0
x0_min = [ 7.2; -0.4; 5.0; 0.035; 0.035; 0.0; 0.0; sqrt(3)/3; sqrt(3)/3; sqrt(3)/3; 0.0; 0.0; 0.0 ]
x0_max = [ 7.2; -0.4; 5.0; 0.035; 0.035; 0.0; 0.0; sqrt(3)/3; sqrt(3)/3; sqrt(3)/3; 0.0; 0.0; 0.0 ]
u0_min = [ NaN; NaN; NaN; NaN; NaN; NaN ]
u0_max = [ NaN; NaN; NaN; NaN; NaN; NaN ]

# final conditions
tf_min = 100
tf_max = 100
xf_min = [ 11.3; 6.0; 4.5; 0.0; 0.0; 0.0; 0.5; -0.5; 0.5; -0.5; 0.0; 0.0; 0.0 ]
xf_max = [ 11.3; 6.0; 4.5; 0.0; 0.0; 0.0; 0.5; -0.5; 0.5; -0.5; 0.0; 0.0; 0.0 ]
uf_min = [ NaN; NaN; NaN; NaN; NaN; NaN ]
uf_max = [ NaN; NaN; NaN; NaN; NaN; NaN ]

# linear path constraints
x_min = [ -12.0; -12.0; -6.0; -0.5; -0.5; -0.5; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0 ]
x_max = [  12.0;  12.0;  6.0;  0.5;  0.5;  0.5;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0 ]
u_min = [ -8e-2; -8e-2; -8e-2; -1e-2; -1e-2; -1e-2 ]
u_max = [  8e-2;  8e-2;  8e-2;  1e-2;  1e-2;  1e-2 ]

# sizes and tolerances
N = 50
Nsub = 15
nx = 13
nu = 6
iter_max = 30
wvc = 1e3
ρ_0 = 0.0
ρ_1 = 0.1
ρ_2 = 0.9
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
csm_plots_freeflyer(prob)
