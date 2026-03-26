using QuantumToolbox, Plots, LaTeXStrings
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
Ntraj = 1000 # number of trajectories
Ncut = 12 # Fock space truncation for the cavity mode

χ = 1 # ħω0/kT -> χ<<1 High temperature limit | χ>>1 Low temperature limit | χ=-1 for zero temperature case

    # Physical Parameters
# ===========================================
Nth = N_BE(χ) # mean thermal photon number in the cavity 
ω_γ = 3.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_γ>>1)
γ0 = 1.0 # unitless coupling system-environment
τ = 2/(1+2*Nth) # dephasing time
dt = τ/100
Tmax = 10.0/(1+2*Nth) # 10 times the relaxation time of steady states

tlist = 0:dt:Tmax # time list
tlist_meas = tlist[1:end-1]

# operators
# -------------------------------------------
a = destroy(Ncut) # cavity mode annihilation operator
H = ω_γ*a'*a # system Hamiltonian
n_op = a'*a
X(θ) = √γ0*(exp(im*θ)*a' + exp(-im*θ)*a) # homodyne measurement operator

sc_ops = [√(1+Nth)*a] # stochastic collapse operators
e_ops = [X(0), X(π/2), n_op] # measured operators
c_ops = [√(Nth)*a'] # absorption part

# initial state
# -------------------------------------------
ρ0 = coherent_dm(Ncut, √2)

    # Simulation
# ===========================================
sol_me = mesolve(H, ρ0, tlist, [sc_ops; c_ops], e_ops=e_ops)
uncond = sol_me.expect # unconditional expectation value

sol_sme = smesolve(H, ρ0, tlist, c_ops, sc_ops[1], e_ops=e_ops, ntraj=Ntraj, store_measurement=Val(true))
cond_meas = sol_sme.measurement # conditional measurement record
cond_mean = sol_sme.expect # conditional expectation value

# println("unconditional: ", shape(uncond))
# println("measure: ", shape(cond_meas))
# println("moy: ", shape(cond_mean))

    # Plotting
# ===========================================
function plotting(flag)
    PATH = "~/Bureau/Github/Stage2026/Homodyne/HARMONIC/"
    ε(i) = real.(uncond[i, :])-real.(sol_sme.expect[i, :])
    titlechi = χ == -1 ? L"\rightarrow +\infty" : "=$(round(χ, digits=2))"
    ek=Dict(:subplot => Dict(:legend_hfactor => 1.5))

    randidx = rand(1:Ntraj, 3) # random indices to chose 3 trajectories to plot
    traj1 = real.(cond_meas[1, randidx[1], :])
    traj2 = real.(cond_meas[1, randidx[2], :])
    traj3 = real.(cond_meas[1, randidx[3], :])
    moyenne3 = (traj1 + traj2 + traj3)/3

    moyenne100 = zeros(length(tlist_meas))
    for i in rand(1:Ntraj, 100)
        moyenne100 += real.(cond_meas[1, i, :])
    end
    moyenne100 /= 100

    if flag
        P_pop = plot(tlist, real.(uncond[3, :]), xlabel=L"\gamma_0 t", ylabel="Mean number of photon", label=L"\langle \hat{n} \rangle", lw=2, legend=:outerright, top_margin=10Plots.mm, extra_kwargs=Dict(:subplot => Dict(:legend_hfactor => 1.5)))
        plot!(P_pop, tlist, real.(sol_sme.expect[3, :]), label=L"\mathbb{E}\left[\langle \hat{n}^{(c)} \rangle\right]", lw=2, ls=:dot, title="Thermal state "*L"\frac{\hbar\omega_0}{k_B T}"*"$titlechi\n$Ntraj trajectories under homodyne detection\n")
        hline!(P_pop, [Nth], label="Thermal limit", ls=:dash, z_order=1, color=palette[4])
        plot!( # error plot
            P_pop, tlist, ε(3), xticks=nothing, tickfontsize=5, xaxis=false,
            inset=bbox(0.35, 0.3, 0.35, 0.35), subplot=2,lw=1, label=nothing, color=palette[3]
        )
        plot!(P_pop, [], [], color=palette[3], label="Error")

        l = @layout [a b c{0.2w} ; d{0.3h} e{0.3h} _]
        Pquad = plot(layout=l) 
        
        plot!(Pquad[1,1], tlist, real.(uncond[1, :]), xlabel=L"\gamma_0 t", ylabel="Quadrature "*L"\langle \hat{x}_\theta \rangle", label=nothing, lw=2, title=L"\hat{x}_{0} = \sqrt{k} ( a + a^\dagger )")
        plot!(Pquad[1,1], tlist, real.(sol_sme.expect[1, :]), label=nothing, lw=2, ls=:dot, top_margin=5Plots.mm)
        plot!(Pquad[2,1], tlist, ε(1), xticks=nothing, tickfontsize=5, xaxis=false, lw=1, label=nothing, color=palette[3], ylabel="Error")

        plot!(Pquad[1,3], [0], [0], label=L"\langle \hat{x}_\theta \rangle", lw=2, legend=:left, framestyle=:none, axis=false, ticks=nothing, background_color_inside=:transparent, legend_font_pointsize=10, legend_font_halign=:center)
        plot!(Pquad[1,3], [0], [0], label=L"\mathbb{E} \left[\langle \hat{x}_\theta^{(c)} \rangle\right]", lw=2, ls=:dot)

        plot!(Pquad[1,2], tlist, real.(uncond[2, :]), label=nothing, xlabel=L"\gamma_0 t", ylabel=nothing, lw=2, title=L"\hat{x}_{\pi/2} = i \sqrt{k} ( a^\dagger - a )")
        plot!(Pquad[1,2], tlist, real.(sol_sme.expect[2, :]), label=nothing, lw=2, ls=:dot, top_margin=5Plots.mm)
        plot!(Pquad[2,2], tlist, ε(2), xticks=nothing, tickfontsize=5, xaxis=false, lw=1, label=nothing, color=palette[3])

        P_sometraj = plot(tlist_meas, traj1, label=L"I^{(1)}", lw=0.5, alpha=0.7, marker=:circle, markersize=1, markerstrokewidth=0, extra_kwargs=Dict(:subplot => Dict(:legend_hfactor => 1.5)))
        # plot!(P_sometraj, tlist_meas, traj2, label=L"I^{(2)}", lw=0.5, alpha=0.7, marker=:circle, markersize=1, markerstrokewidth=0)
        # plot!(P_sometraj, tlist_meas, traj3, label=L"I^{(3)}", lw=0.5, alpha=0.7, marker=:circle, markersize=1, markerstrokewidth=0)
        plot!(P_sometraj, tlist_meas, moyenne3, label=L"\mathbb{E}_{(1)-(3)} [I]", lw=2)
        plot!(P_sometraj, tlist_meas, moyenne100, label=L"\mathbb{E}_{(1)-(100)} [I]", lw=2)
        plot!(P_sometraj, tlist, real.(uncond[1, :]), xlabel=L"\gamma_0 t", ylabel="Photocurrent I(t) & 0-Quadrature", label=L"\langle \hat{x}_0 \rangle", lw=2, title="Example of single trajectories", legend=:outerright)

        titlechi2 = χ == -1 ? "χZeroTemp" : "χ$(round(χ, digits=2))"
        savefig(P_pop, PATH*"harmonic_"*titlechi2*"_homodyne_pop.pdf")
        savefig(Pquad, PATH*"harmonic_"*titlechi2*"_homodyne_quad.pdf")
        savefig(P_sometraj, PATH*"harmonic_"*titlechi2*"_homodyne_traj.pdf")
        closeall()
    end
end

plotting(true)