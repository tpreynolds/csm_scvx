using PyCall, LaTeXStrings
import PyPlot
const plt = PyPlot

function scvx_plot(prob::ScvxProblem)
    fig = plt.figure(figsize=(8,6))
    ax  = plt.gca()

    # Plot SCP solutions
    plt.plot(1,1,label="Initializer", linewidth=2)
    for iter = 2:10
        ax.plot(iter,iter,label="Iterate $(iter - 1)", linewidth=2)
    end

    plt.grid(alpha=0.3)
    plt.draw()

    return fig
end
