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
χ = 2.0 # ħω0/kT -> χ<<1 High temperature limit | χ>>1 Low temperature limit | χ=-1 for zero temperature case

    # Physical Parameters
# ===========================================
const Nth = N_BE(χ) # mean thermal photon number in the cavity 
println("Nth = $(round(Nth, digits=2)), limite stationnaire = $(round(Nth/(1+2Nth), digits=2)), temps relaxation = $(round(1/(1+2Nth), digits=2))")
const ω_γ = 3.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_γ>>1)
const γ0 = 1.0 # unitless coupling system-environment
τ = 2/(1+2*Nth) # dephasing time
const dt = τ/100
println("dt = ", round(dt, digits=4)) 
const Tmax = 10.0/(1+2*Nth) # 10 times the relaxation time of steady states
println("Tmax = ", round(Tmax, digits=4))

tlist = 0:dt:Tmax # time list
tlist_meas = tlist[1:end-1]

# operators
# -------------------------------------------
σx = sigmax()
σy = sigmay()
σz = sigmaz()
σm = sigmam()
Pe = (σz+1)/2 # projector on the excited state -> gives excited population probability in average
H = (ω_γ/2)*σz # Hamiltonian
X(θ) = √γ0*(exp(im*θ)*σm' + exp(-im*θ)*σm) # homodyne measurement operator

sc_ops = [√(1+Nth)*σm] # stochastic collapse operators
e_ops = [X(0), X(π/2), Pe] # measured operators
c_ops = [√(Nth)*σm'] # absorption part

# initial state
# -------------------------------------------
ρ0 = (basis(2, 0)+basis(2,1))/√2

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
    PATH = "~/Bureau/Github/Stage2026/Homodyne/PLOTS/"
    ε(i) = real.(uncond[i, :])-real.(sol_sme.expect[i, :])
    titlechi = χ == -1 ? L"\rightarrow +\infty" : "=$(round(χ, digits=2))"

    if flag
        P_pop = plot(tlist, real.(uncond[3, :]), xlabel=L"\gamma_0 t", ylabel="Excited population", label=L"\rho_{ee}", lw=2, legend=:outerright, top_margin=10Plots.mm)
        plot!(P_pop, tlist, real.(sol_sme.expect[3, :]), label=L"\mathbb{E}\left[\rho^{(c)}_{ee}\right]", lw=2, ls=:dot, title="Thermal state "*L"\frac{\hbar\omega_0}{k_B T}"*"$titlechi\n$Ntraj trajectories under homodyne detection\n")
        annotate!(36, 0.4, text(L"|\psi_0\;\rangle = \frac{|e\rangle+|g\rangle}{\sqrt{2}}", 10, :white))
        plot!( # error plot
            P_pop, tlist, ε(3), xticks=nothing, tickfontsize=5, xaxis=false,
            inset=bbox(0.4, 0.3, 0.35, 0.35), subplots=2,lw=1, label="Error", color=palette[3]
        )

        l = @layout [a b c{0.1w} ; d{0.3h} e{0.3h} _]
        Pquad = plot(layout=l) 
        
        plot!(Pquad[1,1], tlist, real.(uncond[1, :]), xlabel=L"\gamma_0 t", ylabel="Quadrature "*L"\langle \hat{x}_\theta \rangle", label=nothing, lw=2, title=L"\hat{x}_{0} = \sqrt{k} \left( \sigma_- + \sigma_+ \right)")
        plot!(Pquad[1,1], tlist, real.(sol_sme.expect[1, :]), label=nothing, lw=2, ls=:dot, top_margin=5Plots.mm)
        plot!(Pquad[2,1], tlist, ε(1), xticks=nothing, tickfontsize=5, xaxis=false, lw=1, label=nothing, color=palette[3], ylabel="Error")

        plot!(Pquad[1,3], [0], [0], label=L"\langle \hat{x}_\theta \rangle", lw=2, legend=:inside, framestyle=:none, axis=false, ticks=nothing, background_color_inside=:transparent)
        plot!(Pquad[1,3], [0], [0], label=L"\mathbb{E} \left[\langle \hat{x}_\theta^{(c)} \rangle\right]", lw=2, ls=:dot)

        plot!(Pquad[1,2], tlist, real.(uncond[2, :]), label=nothing, xlabel=L"\gamma_0 t", ylabel=nothing, lw=2, title=L"\hat{x}_{\pi/2} = i \sqrt{k} \left( \sigma_+ - \sigma_- \right)")
        plot!(Pquad[1,2], tlist, real.(sol_sme.expect[2, :]), label=nothing, lw=2, ls=:dot, top_margin=5Plots.mm)
        plot!(Pquad[2,2], tlist, ε(2), xticks=nothing, tickfontsize=5, xaxis=false, lw=1, label=nothing, color=palette[3])
        
        titlechi2 = χ == -1 ? "0temp" : "χ$(round(χ, digits=2))"
        savefig(P_pop, PATH*"qubit_"*titlechi2*"_homodyne_pop.pdf")
        savefig(Pquad, PATH*"qubit_"*titlechi2*"_homodyne_quad.pdf")
        closeall()
    end
end

plotting(true)