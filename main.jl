include("scvx.jl")
include("model.jl")
include("csm_plots.jl")

# using PyCall, LaTeXStrings
# import PyPlot
# const plt = PyPlot

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
iter_max = 10
wvc = 1e2
ρ_0 = 0.0
ρ_1 = 0.1
ρ_2 = 0.9
α   = 2.0
β   = 2.0
tr  = 0.5
tr_lb = 0.001
tr_ub = 10.
cvrg_tol = 1e-2
feas_tol = 1e-2

init = ScvxInitBnds(t0_min,t0_max,x0_min,x0_max,u0_min,u0_max)
trgt = ScvxTrgtBnds(tf_min,tf_max,xf_min,xf_max,uf_min,uf_max)
path = ScvxPathBnds(x_min,x_max,u_min,u_max)
bnds = ScvxBnds(init,trgt,path)
ctrl = ScvxParameters(N,Nsub,nx,nu,iter_max,wvc,
                        ρ_0,ρ_1,ρ_2,α,β,tr_lb,tr_ub,
                        cvrg_tol,feas_tol,pars)

# initialize an Scvx problem
prob = scvx_initialize(bnds,ctrl,tr)

# solve the Scvx problem
flag = scvx_solve!(prob)

# plot the results
scvx_plot(prob)

# fig = plt.figure(figsize=(8,6))
# ax  = plt.gca()
#
# # Plot SCP solutions
# plt.plot(1,1,label="Initializer", linewidth=2)
# for iter = 2:10
#     plt.plot(iter,iter,label="Iterate $(iter - 1)", linewidth=2)
# end

# plt.grid(alpha=0.3)
# plt.draw()
