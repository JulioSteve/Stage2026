using QuantumToolbox, Plots, LaTeXStrings
theme(:dark)
palette = theme_palette(:dark)

println("---------------------")
    # Simulation Parameters
# ===========================================
Ntraj = 1000 # number of trajectories
Ncut = 12 # Fock space truncation for the cavity
n_in = 5.0 # mean photon number in the input coherent state

    # Physical Parameters
# ===========================================
const ω0 = π # cavity frequency
const γ = ω0/3 # drive amplitude
const k = ω0/10 # coupling system-environment
const dt = 1/(10*ω0) #; println("dt = ", dt) # small dt compared to the smallest timescale (here 1/ω0)

tlist = 0:dt:30 # time list
tlist_meas = tlist[1:end-1]

# operators
# -------------------------------------------
a = destroy(Ncut) # cavity annihilation operator
H = ω0*a'*a + γ*(a+a') # Hamiltonian
X(θ) = √k*(exp(im*θ)*a + exp(-im*θ)*a') # homodyne measurement operator

c_ops = [] # collapse operators - needed by QuantumToolbox but not used in this simulation
sc_ops = [√k*a] # stochastic collapse operators
e_ops = [X(0), X(π/2), a'*a] # measured operators

# initial state
# -------------------------------------------
ρ0 = coherent(Ncut, √n_in) # coherent state

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


#     Plotting
#  ===========================================
function plotting(flag)
    PATH = "~/Bureau/Github/Stage2026/Homodyne/PLOTS/"
    if flag
        # P = plot(tlist_meas, real.(cond_meas[1, rand(1:Ntraj), :]), label="measurement record", lw=1)

        P_pop = plot(tlist, real.(uncond[3, :]), xlabel="t", ylabel=L"\langle \mathbb{P}_{ee} \rangle", label="mesolve", lw=2)
        plot!(P_pop, tlist, real.(sol_sme.expect[3, :]), label="smesolve", lw=2, ls=:dash, legend=:outerright)
        yticks!(0:6)
        plot!( # error plot
            P_pop, tlist, real.(uncond[3, :])-real.(sol_sme.expect[3, :]), 
            inset=bbox(0.4, 0.1, 0.35, 0.35), subplot=2, lw=1, label="Error", color=palette[3]
        )

        P_quadX = plot(tlist, real.(uncond[1, :]), xlabel="t", ylabel=L"\langle \hat{x}_{0} \;\rangle", label="mesolve", lw=2, legend=:outerright, title=L"\hat{x}_{0} = \sqrt{k} \left( \hat{a} + \hat{a}^\dagger \right)")
        plot!(P_quadX, tlist, real.(sol_sme.expect[1, :]), label="smesolve", lw=2, ls=:dash)
        # ylims!(-0.2,0.2)
        plot!( # error plot
            P_quadX, tlist, real.(uncond[1, :])-real.(sol_sme.expect[1, :]), 
            inset=bbox(0.5, 0.05, 0.25, 0.25), subplot=2, lw=1, label="Error", color=palette[3]
        )

        P_quadY = plot(tlist, real.(uncond[2, :]), xlabel="t", ylabel=L"\langle \hat{x}_{\pi/2} \;\rangle", label="mesolve", lw=2, legend=:outerright, title=L"\hat{x}_{\pi/2} = i \sqrt{k} \left( \hat{a}^dagger - \hat{a} \right)")
        plot!(P_quadY, tlist, real.(sol_sme.expect[2, :]), label="smesolve", lw=2, ls=:dash)
        # ylims!(-0.2,0.2)
        plot!( # error plot
            P_quadY, tlist, real.(uncond[1, :])-real.(sol_sme.expect[1, :]), 
            inset=bbox(0.5, 0.05, 0.25, 0.25), subplot=2, lw=1, label="Error", color=palette[3]
        )

        savefig(P_pop, PATH*"harmonic_homodyne_pop.pdf")
        savefig(P_quadX, PATH*"harmonic_homodyne_quadX.pdf")
        savefig(P_quadY, PATH*"harmonic_homodyne_quadY.pdf")
        closeall()
    end
end

plotting(true)