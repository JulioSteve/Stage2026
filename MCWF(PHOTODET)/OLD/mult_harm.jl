using QuantumToolbox, Plots, LaTeXStrings, Statistics
theme(:dark)
palette = theme_palette(:dark)

χlist = [0.25:0.25:2; -1]

for (idx,χ) in enumerate(χlist)

    function N_BE(x)
        if x == -1
            return 0
        else
            return 1/(exp(x)-1)
        end
    end

    println("---------------------")

        # Simulation Parameters
    # ===========================================
    Ntraj = 10_000 # number of trajectories
    Ncut = 20 # Fock space truncation for the cavity mode
    n_input = 2.0 # Starting level for the initial coherent state

    titlechi(x) = x==-1 ? L"\frac{\hbar \omega_0}{k_B T}\rightarrow+\infty" : L"\frac{\hbar \omega_0}{k_B T}="*"$x"
    titlechiplot(x) = x==-1 ? "χ_zeroTemp" : "χ_$x"

    Nth = N_BE(χ) # mean thermal photon number in the cavity 
    ω_γ = 3.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_γ>>1)
    γ0 = 1.0 # unitless coupling system-environment
    τ = 2/(1+2*Nth) # dephasing time
    dt = τ/100
    Tmax = 10.0/(1+2*Nth) # 10 times the relaxation time of steady states
    println("dt = $(round(dt, digits=5)) Tmax = $(round(Tmax, digits=2))")

    tlist = 0:dt:Tmax # time list

    a = destroy(Ncut) # cavity mode annihilation operator
    H = ω_γ*a'*a # system Hamiltonian
    n_op = a'*a

    c_ops = [√(γ0*(1+Nth))*a, √(γ0*Nth)*a'] # operators of the master equation
    e_ops = [n_op] # measured operators

    ψ0 = coherent(Ncut, √n_input)

    sol_me = mesolve(H, ψ0, tlist, c_ops, e_ops=e_ops)
    me_mean = real.(sol_me.expect) # unconditional expectation value

    sol_mc = mcsolve(H, ψ0, tlist, c_ops, e_ops=e_ops, ntraj=Ntraj, keep_runs_results=Val(false))
    mc_trajs = real.(sol_mc.expect) # conditional expectation value

    PATH = "~/Bureau/Github/Stage2026/MCWF/HARMONIC/"
    condmean = mean(mc_trajs, dims=1)[1,:]
    P1 = plot()
    hline!(P1, [Nth], label="thermal limit", ls=:dot, lw=2, title=titlechi(χ))
    plot!(P1, tlist, me_mean[1,:], label="me_meanitional", legend=:outerright, lw=3)
    plot!(P1, tlist, condmean, label="conditional", lw=2, ls=:dash)

    savefig(P1, PATH*"$(idx)_$(titlechiplot(χ))_plot.pdf")
    closeall()
end