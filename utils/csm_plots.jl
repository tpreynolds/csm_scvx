using PyCall, LaTeXStrings, Colors
import PyPlot
const plt = PyPlot

struct CSMPlotCol
    red
    blue
    darkblue
    green
    gray
    purple
    gold
    darkgold
    magenta
    cyan
end
function CSMPlotCol()
    red = [234;61;37]/255
    blue = [0;32;91]/255
    darkblue = [4;28;44]/255
    green = [10;134;61]/255
    gray = [153;153;154]/255
    purple = [51;0;111]/255
    gold = [232;211;162]/255
    darkgold = [145;123;76]/255
    magenta = [1;0;1]
    cyan = [0;1;1]
    return CSMPlotCol(red,blue,darkblue,green,gray,purple,gold,darkgold,magenta,cyan)
end

struct CSMPlotFmt
    col::CSMPlotCol
    circle::Array{Float64,2}
    markersize::Integer
    gridalpha::Float64
    figsize::Tuple{Int64,Int64}
    dblwide::Tuple{Int64,Int64}
    lw::Integer
    labelsize::Integer
    fontsize::Integer
    titlesize::Integer
end
function CSMPlotFmt()
    col     = CSMPlotCol()
    angles  = LinRange(0,2*pi,100)
    circle  = zeros(2,100)
    for k = 1:100
        circle[:,k] = [cos(angles[k]);sin(angles[k])]
    end
    markersize  = 5
    gridalpha   = 0.3
    figsize     = (8,6)
    dblwide     = (12,6)
    linewidth   = 2
    labelsize   = 14
    fontsize    = 16
    titlesize   = 18
    return CSMPlotFmt(col,circle,markersize,gridalpha,
                        figsize,dblwide,linewidth,labelsize,fontsize,titlesize)
end

function csm_plots_quad(prob::ScvxProblem)
    # integrate nonlinear dynamics
    t_grid = LinRange(0,prob.new_sol.tf,prob.pars.N)
    T_grid = LinRange(0,prob.new_sol.tf,250)
    x0     = prob.new_sol.state[:,1]
    u      = prob.new_sol.control
    f(t,y) = dynamics(t,y,u,t_grid,prob.pars.mdl_pars)
    X      = rk4(f,T_grid,x0)

    fmt = CSMPlotFmt()
    csm_plot_quad_trj(prob,fmt,X)
    csm_plot_quad_tilt(prob,fmt,X)
    csm_plot_quad_accel(prob,fmt)
    csm_plot_quad_alltrjs(prob,fmt)
    return nothing
end

function csm_plots_freeflyer(prob::ScvxProblem)
    # integrate nonlinear dynamics
    t_grid = LinRange(0,prob.new_sol.tf,prob.pars.N)
    T_grid = LinRange(0,prob.new_sol.tf,250)
    x0     = prob.new_sol.state[:,1]
    u      = prob.new_sol.control
    f(t,y) = dynamics(t,y,u,t_grid,prob.pars.mdl_pars)
    X      = rk4(f,T_grid,x0)


    fmt = CSMPlotFmt()
    csm_plot_freeflyer_trj(prob,fmt,X)
    csm_plot_freeflyer_ctrl(prob,fmt)
    csm_plot_freeflyer_attitude(prob,fmt,X)
    # csm_plot_freeflyer_alltrjs(prob,fmt)
    return nothing
end

function csm_plot_quad_trj(prob::ScvxProblem,fmt::CSMPlotFmt,X)
    x   = prob.new_sol.state
    u   = prob.new_sol.control
    tf  = prob.new_sol.tf
    circle = fmt.circle
    scl = 0.3

    # create figure and sets of axes
    fig = plt.figure(figsize=fmt.figsize)
    ax  = plt.gca()
    axins1 = ax.inset_axes([0.725, 0.05, 0.25, 0.25])
    x1, x2, y1, y2 = 0.3, 0.8, 1.75, 2.5
    axins1.set_xlim(x1, x2)
    axins1.set_ylim(y1, y2)
    axins2 = ax.inset_axes([0.725, 0.375, 0.25, 0.25])
    x1, x2, y1, y2 = 2.2, 2.6, 4.4, 5.2
    axins2.set_xlim(x1, x2)
    axins2.set_ylim(y1, y2)

    # plot obstacles
    for i = 1:pars.obsN
        H = I(2)/pars.obsiH[1:2,1:2,i]
        c = pars.obsC[1:2,i]

        obs = H * circle .+ c
        ax.plot(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.8,
                linewidth=1,linestyle="-")
        axins1.plot(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.8,
                linewidth=1,linestyle="-")
        axins2.plot(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.8,
                linewidth=1,linestyle="-")
        if i==1
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-",label=L"Obstacle")
            else
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-")
        end
        axins1.fill_between(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.1,
                linewidth=1,linestyle="-")
        axins2.fill_between(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.1,
                linewidth=1,linestyle="-")
    end

    # plot discrete solution
    ax.plot(X[1,:],X[2,:],
            label="Integrated Trajectory",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(x[1,:],x[2,:],
            label="SCvx solution",
            marker="o",color=fmt.col.blue,linestyle="",
            markersize=fmt.markersize)
    axins1.plot(X[1,:],X[2,:],
            label="Integrated Trajectory",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    axins1.plot(x[1,:],x[2,:],
            label="SCvx solution",
            marker="o",color=fmt.col.blue,linestyle="",
            markersize=fmt.markersize)
    axins2.plot(X[1,:],X[2,:],
            label="Integrated Trajectory",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    axins2.plot(x[1,:],x[2,:],
            label="SCvx solution",
            marker="o",color=fmt.col.blue,linestyle="",
            markersize=fmt.markersize)
    # add thrust vectors
    udir = u[1:2,1]/norm(u[1:2,1])
    xs = [ x[1,1], x[1,1]+scl*udir[1] ]
    ys = [ x[2,1], x[2,1]+scl*udir[2] ]
    lines = Any[collect(zip(xs,ys))]
    for k = 2:prob.pars.N
        udir = u[1:2,k]/norm(u[1:2,k])
        xs = [ x[1,k], x[1,k]+scl*udir[1] ]
        ys = [ x[2,k], x[2,k]+scl*udir[2] ]
        push!(lines,collect(zip(xs,ys)))
    end
    horz_thrust_vecs = plt.matplotlib.collections.LineCollection(lines,
                                                    color=fmt.col.green,
                                                    label="Thrust Direction",
                                                    linewidth=fmt.lw)
    ax.add_collection(horz_thrust_vecs)
    # axins1.add_collection(horz_thrust_vecs)
    # axins2.add_collection(horz_thrust_vecs)

    # including these will add some gray lines that point to where the insets
    # come from
    # ax.indicate_inset_zoom(axins1)
    # ax.indicate_inset_zoom(axins2)

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    plt.xlim(-0.5,3)
    plt.ylim(-0.5,6.5)
    plt.grid(alpha=fmt.gridalpha)
    plt.title("Final Quadcopter Trajectory",fontsize=fmt.titlesize)
    plt.xlabel("E [m]",fontsize=fmt.fontsize)
    plt.ylabel("N [m]",fontsize=fmt.fontsize)
    plt.tight_layout()
    ax.legend(fontsize=fmt.fontsize)
    plt.show()

    return nothing
end

function csm_plot_quad_tilt(prob::ScvxProblem,fmt::CSMPlotFmt,X)
    u  = prob.new_sol.control
    tf = prob.new_sol.tf
    N  = prob.pars.N
    Ni = size(X,2)
    T  = LinRange(0,tf,N)
    Ti = LinRange(0,tf,Ni)
    tilt_max = rad2deg(prob.pars.mdl_pars.tilt_max)

    fig = plt.figure(figsize=fmt.figsize)
    ax  = plt.gca()

    # plot constraint
    ax.plot([0;tf],[tilt_max;tilt_max],linestyle="--",
                    linewidth=fmt.lw,
                    color=fmt.col.red)
    ax.fill_between([0;tf],[tilt_max;tilt_max],[tilt_max+10;tilt_max+10],
                    facecolor=fmt.col.red,alpha=0.1)
    # add discrete/continuous tilt angles
    tilt_d = zeros(N)
    tilt_c = zeros(Ni)
    for k = 1:N
        umag      = norm(u[1:3,k])
        tilt_d[k] = acosd(u[3,k]/umag)
    end
    for k = 1:Ni
        uk        = interp_vec(Ti[k],u,T)
        umag      = norm(uk[1:3])
        tilt_c[k] = acosd(uk[3]/umag)
    end
    ax.plot(Ti,tilt_c,color=fmt.col.blue,
                linestyle="-",linewidth=fmt.lw)
    ax.plot(T,tilt_d,color=fmt.col.blue,marker="o",
                markersize=fmt.markersize,
                linestyle="none")

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    plt.xlim(0,tf)
    plt.ylim(0,70)
    plt.grid(alpha=fmt.gridalpha)
    plt.title("Tilt Angle",fontsize=fmt.titlesize)
    plt.xlabel("Time [s]",fontsize=fmt.fontsize)
    plt.ylabel(L"$\theta$ [deg]",fontsize=fmt.fontsize)
    plt.tight_layout()
    plt.show()
    return nothing
end

function csm_plot_quad_accel(prob::ScvxProblem,fmt::CSMPlotFmt)
    u  = prob.new_sol.control
    tf = prob.new_sol.tf
    N  = prob.pars.N
    T  = LinRange(0,tf,N)
    u_max = prob.pars.mdl_pars.u_nrm_max
    u_min = prob.pars.mdl_pars.u_nrm_min

    umag = zeros(N)
    for k = 1:N
        umag[k] = norm(u[1:3,k])
    end

    fig = plt.figure(figsize=fmt.figsize)
    ax  = plt.gca()

    # plot acceleration limits
    ax.plot([0;tf],[u_max;u_max],linestyle="--",linewidth=fmt.lw,
                    color=fmt.col.red)
    ax.plot([0;tf],[u_min;u_min],linestyle="--",linewidth=fmt.lw,
                    color=fmt.col.red)
    ax.fill_between([0;tf],[u_max;u_max],[u_max+10;u_max+10],
                    facecolor=fmt.col.red,alpha=0.1)
    ax.fill_between([0;tf],[u_min;u_min],[0;0],
                    facecolor=fmt.col.red,alpha=0.1)

    # plot acceleration magnitude
    ax.plot(T,umag,linestyle="-",linewidth=fmt.lw,
                    color=fmt.col.blue,
                    marker="o",markersize=fmt.markersize)

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    plt.xlim(0,tf)
    plt.ylim(0,25)
    plt.grid(alpha=fmt.gridalpha)
    plt.title("Cmd. Acceleration",fontsize=fmt.titlesize)
    plt.xlabel("Time [s]",fontsize=fmt.fontsize)
    plt.ylabel(L"$\|\|u\|\|_2~\mathrm{[m/s^2]}$",fontsize=fmt.fontsize)
    plt.tight_layout()
    plt.show()
    return nothing
end

function csm_plot_quad_alltrjs(prob::ScvxProblem,fmt::CSMPlotFmt)
    x_all = prob.all_trj.state
    circle = fmt.circle
    col1 = fmt.col.cyan
    col2 = fmt.col.magenta
    cols = zeros(3,prob.solved)
    for i = 1:3
        cols[i,:] = LinRange(col1[i],col2[i],prob.solved)
    end

    fig = plt.figure(figsize=fmt.figsize)
    ax  = plt.gca()

    # plot obstacles
    for i = 1:pars.obsN
        H = I(2)/pars.obsiH[1:2,1:2,i]
        c = pars.obsC[1:2,i]

        obs = H * circle .+ c
        ax.plot(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.8,
                linewidth=1,linestyle="-")
        ax.fill_between(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.1,
                linewidth=1,linestyle="-")
    end

    # plot discrete solutions
    for iter = 1:prob.solved
        ax.plot(x_all[1,:,iter],x_all[2,:,iter],
                label="Iteration $(iter)",
                marker="o",color=cols[:,iter],linestyle="-",
                markersize=fmt.markersize)
    end

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    plt.xlim(-0.5,3)
    plt.ylim(-0.5,6.5)
    plt.grid(alpha=fmt.gridalpha)
    plt.title("Final Quadcopter Trajectory",fontsize=fmt.titlesize)
    plt.xlabel("E [m]",fontsize=fmt.fontsize)
    plt.ylabel("N [m]",fontsize=fmt.fontsize)
    plt.tight_layout()
    ax.legend(fontsize=fmt.fontsize)
    plt.show()
end

function csm_plot_freeflyer_trj(prob::ScvxProblem,fmt::CSMPlotFmt,X)
    x   = prob.new_sol.state
    u   = prob.new_sol.control
    tf  = prob.new_sol.tf
    circle = fmt.circle
    scl = 0.3

    fig = plt.figure(figsize=fmt.dblwide)
    ax = plt.subplot(121)
    ## Plot X-Y plane
    # plot obstacles
    for i = 1:pars.obsN
        H = I(2)/pars.obsiH[1:2,1:2,i]
        c = pars.obsC[1:2,i]

        obs = H * circle .+ c
        ax.plot(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.8,
                linewidth=1,linestyle="-")
        if i==1
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-",label="Obstacle")
            else
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-")
        end
    end

    # plot discrete solution
    ax.plot(X[1,:],X[2,:],
            label=L"Integrated\ Trajectory",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(x[1,:],x[2,:],
            label=L"SCvx\ solution",
            marker="o",color=fmt.col.blue,linestyle="",
            markersize=fmt.markersize)

    # add thrust vectors
    udir = u[:,1]/norm(u[:,1])
    xs = [ x[1,1], x[1,1]+scl*udir[1] ]
    ys = [ x[2,1], x[2,1]+scl*udir[2] ]
    lines = Any[collect(zip(xs,ys))]
    for k = 2:prob.pars.N
        udir = u[:,k]/norm(u[:,k])
        xs = [ x[1,k], x[1,k]+scl*udir[1] ]
        ys = [ x[2,k], x[2,k]+scl*udir[2] ]
        push!(lines,collect(zip(xs,ys)))
    end
    horz_thrust_vecs = plt.matplotlib.collections.LineCollection(lines,
                                                    color=fmt.col.green,
                                                    label=L"Thrust\ Direction",
                                                    linewidth=fmt.lw)
    ax.add_collection(horz_thrust_vecs)
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(6,12)
    ax.set_ylim(-2,8)
    ax.set_xlabel(L"X\ [m]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"Y\ [m]",fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=fmt.fontsize)

    ## Plot Y-Z plane
    ax = plt.subplot(122)
    # plot obstacles
    for i = 1:pars.obsN
        H = I(2)/pars.obsiH[2:3,2:3,i]
        c = pars.obsC[2:3,i]

        obs = H * circle .+ c
        ax.plot(obs[1,:],obs[2,:],
                color=fmt.col.red,alpha=0.8,
                linewidth=1,linestyle="-")
        if i==1
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-",label=L"Obstacle")
            else
            ax.fill_between(obs[1,:],obs[2,:],
                    color=fmt.col.red,alpha=0.1,
                    linewidth=1,linestyle="-")
        end
    end

    # plot discrete solution
    ax.plot(X[2,:],X[3,:],
            label=L"Integrated\ Trajectory",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(x[2,:],x[3,:],
            label=L"SCvx\ solution",
            marker="o",color=fmt.col.blue,linestyle="",
            markersize=fmt.markersize)

    # add thrust vectors
    udir = u[:,1]/norm(u[:,1])
    xs = [ x[2,1], x[2,1]+scl*udir[2] ]
    ys = [ x[3,1], x[3,1]+scl*udir[3] ]
    lines = Any[collect(zip(xs,ys))]
    for k = 2:prob.pars.N
        udir = u[:,k]/norm(u[:,k])
        xs = [ x[2,k], x[2,k]+scl*udir[2] ]
        ys = [ x[3,k], x[3,k]+scl*udir[3] ]
        push!(lines,collect(zip(xs,ys)))
    end
    horz_thrust_vecs = plt.matplotlib.collections.LineCollection(lines,
                                                    color=fmt.col.green,
                                                    label=L"Thrust\ Direction",
                                                    linewidth=fmt.lw)
    ax.add_collection(horz_thrust_vecs)
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(-2,8)
    ax.set_ylim(4,6)
    ax.set_xlabel(L"Y\ [m]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"Z\ [m]",fontsize=fmt.fontsize)
    # ax.legend(fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)

    fig.suptitle(L"Final\ FreeFlyer\ Trajectory",fontsize=fmt.titlesize)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    return nothing
end

function csm_plot_freeflyer_ctrl(prob::ScvxProblem,fmt::CSMPlotFmt)
    id_F = prob.pars.mdl_pars.id_F
    id_M = prob.pars.mdl_pars.id_M
    F = prob.new_sol.control[id_F,:]
    M = prob.new_sol.control[id_M,:]
    tf  = prob.new_sol.tf
    T  = LinRange(0,tf,prob.pars.N)

    F_nrm_max = prob.pars.mdl_pars.F_nrm_max
    M_nrm_max = prob.pars.mdl_pars.M_nrm_max

    # compute control input norms for plotting
    F_nrm = zeros(prob.pars.N)
    M_nrm = zeros(prob.pars.N)
    for k = 1:prob.pars.N
        F_nrm[k] = norm(F[:,k])
        M_nrm[k] = norm(M[:,k])
    end

    fig = plt.figure(figsize=fmt.figsize)
    ax = plt.subplot(121)
    ax.plot(T,F[1,:],label=L"x",
            color=fmt.col.red,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,F[2,:],label=L"y",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,F[3,:],label=L"z",
            color=fmt.col.green,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,F_nrm,label=L"||F||",
            color=[0;0;0],
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot([0;tf],[F_nrm_max;F_nrm_max],label=L"||F||_{2,max}",
            color=[1;0;0],
            linestyle="--",
            linewidth=fmt.lw)

    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,100)
    # ax.set_ylim(4,6)
    ax.set_xlabel(L"Time\ [s]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"Thrust\ [N]",fontsize=fmt.fontsize)
    # ax.legend(fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=fmt.fontsize)

    ## Plot moment/torque trajectory
    ax = plt.subplot(122)
    ax.plot(T,M[1,:],label=L"x",
            color=fmt.col.red,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,M[2,:],label=L"y",
            color=fmt.col.blue,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,M[3,:],label=L"z",
            color=fmt.col.green,linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,M_nrm,label=L"||M||",
            color=[0;0;0],
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot([0;tf],[M_nrm_max;M_nrm_max],label=L"||M||_{2,max}",
            color=[1;0;0],
            linestyle="--",
            linewidth=fmt.lw)
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,100)
    # ax.set_ylim(4,6)
    ax.set_xlabel(L"Time\ [s]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"Torque\ [N]",fontsize=fmt.fontsize)
    # ax.legend(fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=fmt.fontsize)

    fig.suptitle(L"Final\ FreeFlyer\ Controls",fontsize=fmt.titlesize)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    return nothing
end

function csm_plot_freeflyer_attitude(prob::ScvxProblem,fmt::CSMPlotFmt,X)
    id_q = prob.pars.mdl_pars.id_q
    id_w = prob.pars.mdl_pars.id_w
    quat_d = prob.new_sol.state[id_q,:]
    quat_c = X[id_q,:]
    wB_d   = prob.new_sol.state[id_w,:]
    wB_c   = X[id_w,:]
    rpy_d  = quat_2_rpy(quat_d)
    rpy_c  = quat_2_rpy(quat_c)

    # convert to degrees
    rad2deg_arr!(rpy_d)
    rad2deg_arr!(rpy_c)
    rad2deg_arr!(wB_d)
    rad2deg_arr!(wB_c)

    tf  = prob.new_sol.tf
    T   = LinRange(0,tf,prob.pars.N)
    Ti  = LinRange(0,tf,size(X,2))

    w_nrm_max = prob.pars.mdl_pars.w_nrm_max

    fig = plt.figure(figsize=fmt.figsize)
    ax = plt.subplot(121)
    ax.plot(T,rpy_d[1,:],
            color=fmt.col.red,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize)
    ax.plot(Ti,rpy_c[1,:],
            label=L"roll",
            color=fmt.col.red,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,rpy_d[2,:],
            color=fmt.col.green,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize)
    ax.plot(Ti,rpy_c[2,:],
            label=L"pitch",
            color=fmt.col.green,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,rpy_d[3,:],
            color=fmt.col.blue,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize)
    ax.plot(Ti,rpy_c[3,:],
            label=L"yaw",
            color=fmt.col.blue,
            linestyle="-",
            linewidth=fmt.lw)
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,100)
    ax.set_ylim(-180,180)
    ax.set_xlabel(L"Time\ [s]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"Attitude\ [deg]",fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=fmt.fontsize)

    ax = plt.subplot(122)
    ax.plot(T,wB_d[1,:],
            color=fmt.col.red,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize,
            linewidth=fmt.lw)
    ax.plot(Ti,wB_c[1,:],
            label=L"ω_x",
            color=fmt.col.red,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,wB_d[2,:],
            color=fmt.col.green,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize,
            linewidth=fmt.lw)
    ax.plot(Ti,wB_c[2,:],
            label=L"ω_y",
            color=fmt.col.green,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot(T,wB_d[3,:],
            color=fmt.col.blue,
            linestyle="none",
            marker="o",
            markersize=fmt.markersize,
            linewidth=fmt.lw)
    ax.plot(Ti,wB_c[3,:],
            label=L"ω_z",
            color=fmt.col.blue,
            linestyle="-",
            linewidth=fmt.lw)
    ax.plot([0;tf],[rad2deg(w_nrm_max);rad2deg(w_nrm_max)],label=L"||ω||_{∞,max}",
            color=[1;0;0],linestyle="--",
            linewidth=fmt.lw)
    ax.tick_params(axis="both", which="major", labelsize=fmt.labelsize)
    ax.set_xlim(0,100)
    # ax.set_ylim(-1,180)
    ax.set_xlabel(L"Time\ [s]",fontsize=fmt.fontsize)
    ax.set_ylabel(L"Angular Rate\ [deg/s]",fontsize=fmt.fontsize)
    ax.grid(alpha=fmt.gridalpha)
    ax.legend(fontsize=fmt.fontsize)

    fig.suptitle(L"Final\ FreeFlyer\ Attitude",fontsize=fmt.titlesize)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    return nothing
end
