using QuantumToolbox, Plots, LaTeXStrings
theme(:dark)

println("---------------------")
    # Simulation Parameters
# ===========================================
Ntraj = 100 # number of trajectories
abs = false # whether to include absorption operators

    # Physical Parameters
# ===========================================
const ω0 = 2π
const k = 0.1 # coupling system-environment

tlist = 0:0.005:30 # time list
tlist_meas = tlist[1:end-1]

# operators
# -------------------------------------------
σx = sigmax()
σy = sigmay()
σz = sigmaz()
σm = sigmam()
Pe = (σz+1)/2 # projector on the excited state -> gives excited population probability in average
H = ω0*σz/2 # Hamiltonian
X(θ) = √k*(exp(im*θ)*σm' + exp(-im*θ)*σm) # homodyne measurement operator

sc_ops = [√k*σm] # stochastic collapse operators
e_ops = [X(0), X(π/2), Pe] # measured operators

if abs == true
    c_ops = [√k*σm'] # absorption operators
else
    c_ops = []
end

# initial state
# -------------------------------------------
ρ0 = basis(2, 0) # excited state

    # Simulation
# ===========================================
sol_me = mesolve(H, ρ0, tlist, [sc_ops; c_ops], e_ops=e_ops)
uncond = sol_me.expect # unconditional expectation value

sol_sme = smesolve(H, ρ0, tlist, c_ops, sc_ops[1], e_ops=e_ops, ntraj=Ntraj, store_measurement=Val(true))
cond_meas = sol_sme.measurement # conditional measurement record
cond_mean = sol_sme.expect # conditional expectation value

println("unconditional: ", shape(uncond))
println("measure: ", shape(cond_meas))
println("moy: ", shape(cond_mean))

    # Plotting
# ===========================================
function plotting(flag)
    PATH = "~/Bureau/Github/Stage2026/Homodyne/PLOTS/"
    if flag
        # P = plot(tlist_meas, real.(cond_meas[1, rand(1:Ntraj), :]), label="measurement record", lw=1)

        P_pop = plot(tlist, real.(uncond[3, :]), xlabel="t", ylabel=L"\langle \mathbb{P}_{ee} \rangle", label="mesolve", lw=2)
        plot!(P_pop, tlist, real.(sol_sme.expect[3, :]), label="smesolve", lw=2)
        ylims!(0,1)

        P_quadX = plot(tlist, real.(uncond[1, :]), xlabel="t", ylabel=L"\langle X \rangle", label="mesolve", lw=2)
        plot!(P_quadX, tlist, real.(sol_sme.expect[1, :]), label="smesolve", lw=2)
        ylims!(-0.2,0.2)

        P_quadY = plot(tlist, real.(uncond[2, :]), xlabel="t", ylabel=L"\langle Y \rangle", label="mesolve", lw=2)
        plot!(P_quadY, tlist, real.(sol_sme.expect[2, :]), label="smesolve", lw=2)
        ylims!(-0.2,0.2)

        savefig(P_pop, PATH*"main_homodyne_pop.pdf")
        savefig(P_quadX, PATH*"main_homodyne_quadX.pdf")
        savefig(P_quadY, PATH*"main_homodyne_quadY.pdf")
        closeall()
    end
end

plotting(true)