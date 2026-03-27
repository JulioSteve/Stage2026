using QuantumToolbox, Plots, LaTeXStrings, Statistics
theme(:dark)
palette = theme_palette(:dark)

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

χ = 1

    # Physical Parameters
# ===========================================
Nth = N_BE(χ) # mean thermal photon number in the cavity 
println("occupation thermique: $Nth")
ω_γ = 3.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_γ>>1)
γ0 = 1.0 # unitless coupling system-environment
τ = 2/(1+2*Nth) # dephasing time
dt = τ/1000
Tmax = 10.0/(1+2*Nth) # 10 times the relaxation time of steady states

tlist = 0:dt:Tmax # time list
# tlist_meas = tlist[1:end-1]

println("dt=$(round(dt, digits=5)), Tmax=$(round(Tmax, digits=2)), time steps:$(length(tlist))")

# operators
# -------------------------------------------
a = destroy(Ncut) # cavity mode annihilation operator
H = ω_γ*a'*a # system Hamiltonian
n_op = a'*a

c_ops = [√(γ0*(1+Nth))*a, √(γ0*Nth)*a'] # operators of the master equation
e_ops = [n_op] # measured operators

# initial state
# -------------------------------------------
ψ0 = coherent(Ncut, √n_input)

    # Simulation
# ===========================================
sol_me = mesolve(H, ψ0, tlist, c_ops, e_ops=e_ops)
me_mean = real.(sol_me.expect) # unconditional expectation value

sol_mc = mcsolve(H, ψ0, tlist, c_ops, e_ops=e_ops, ntraj=Ntraj, keep_runs_results=Val(true))
mc_trajs = real.(sol_mc.expect) # conditional expectation value

println("unconditional: ", shape(me_mean))
println("conditional: ", shape(mc_trajs))

    # Plotting
# ===========================================
function plotting(flag, uncond, condtrajs)
    PATH = "~/Bureau/Github/Stage2026/MCWF/HARMONIC/"
    condmean = mean(condtrajs, dims=2)[1,1,:]

    if flag
        P1 = plot()
        hline!(P1, [Nth], label="thermal limit", ls=:dot, lw=2)
        for i in Int.([1, Ntraj÷2, Ntraj])
            plot!(P1, tlist, condtrajs[1,i,:], label="$(i)"*L"^{th} trajectory", lw=1, alpha=0.4, ls=:dashdot)
        end
        plot!(P1, tlist, uncond[1,:], label="unconditional", legend=:outerright, lw=3)
        plot!(P1, tlist, condmean, label="conditional", lw=2, ls=:dash)

        savefig(P1, PATH*"plot.pdf")
        closeall()
    end
end

plotting(true, me_mean, mc_means)