using QuantumToolbox, Plots, LaTeXStrings

println("---------------------")
    # Simulation Parameters
# ===========================================
Ntraj = 100 # number of trajectories
abs = true # whether to include absorption operators

    # Physical Parameters
# ===========================================
const ω0 = 2π
const k = 0.1 # coupling system-environment

tlist = 0:0.001:30 # time list
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
e_ops = [X(0), Pe] # measured operators

if abs == true
    c_ops = [√k*σm'] # absorption operators
else
    c_ops = []
end

# initial state
# -------------------------------------------
ψ0 = basis(2, 0) # excited state

    # Simulation
# ===========================================
sol_me = mesolve(H, ψ0, tlist, [sc_ops; c_ops], e_ops=e_ops)
uncond = sol_me.expect # unconditional expectation value

sol_sme = smesolve(H, ψ0, tlist, c_ops, sc_ops[1], e_ops=e_ops, ntraj=Ntraj, store_measurement=Val(true))
cond_meas = sol_sme.measurement # conditional measurement record
cond_mean = sol_sme.expect # conditional expectation value

println("unconditional: ", shape(uncond))
println("measure: ", shape(cond_meas))
println("moy: ", shape(cond_mean))

    # Plotting
# ===========================================
function plotting(flag)
    PATH = "~/Bureau/Github/Stage2026/Homodyne/"
    if flag
        # P = plot(tlist_meas, real.(cond_meas[1, rand(1:Ntraj), :]), label="measurement record", lw=1)

        P = plot(tlist, real.(uncond[2, :]), xlabel="t", ylabel=L"\langle \mathbb{P}_{ee} \rangle", label="mesolve", lw=2)
        plot!(P, tlist, real.(sol_sme.expect[2, :]), label="smesolve", lw=2)
        ylims!(-1,1)

        P2 = plot(tlist, real.(uncond[1, :]), xlabel="t", ylabel=L"\langle X_0 \rangle", label="mesolve", lw=2)
        plot!(P2,tlist, real.(sol_sme.expect[1, :]), label="smesolve", lw=2)
        ylims!(-1,1)

        savefig(P, PATH*"homodyne_pop.pdf")
        savefig(P2, PATH*"homodyne_quad.pdf")
        closeall()
    end
end

plotting(true)