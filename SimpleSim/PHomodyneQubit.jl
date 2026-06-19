ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)


Ω = 10 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the coupling constant k. Weak coupling implies Ω>>1.
dτ = 2π/(3000*Ω) # By definition T0=2π/ω0 and to work unitless: kT0=2πk/ω0. We choose to take 3000 points per oscillation.
# We choose a smaller dτ to reduce the numerical errors. The stochastic solver is greatly impacted (precision and convergence) by the size of dτ.
τlist = 0:dτ:10

σm = sigmam()
σz = sigmaz()
c(θ) = σm*exp(-im*θ)
quad(θ) = (σm*exp(-im*θ) + σm'*exp(im*θ))/√2
ψ0 = (fock(2, 0)+fock(2, 1))./√2
ρ0 = ψ0*ψ0' # ρ0 = |ψ0><ψ0|, where |ψ0> = (|e>+|g>)/√2
H = Ω*σz./2 # Hamiltonian of the system (unitless)

# Simulations
Ntraj = 500

homP = smesolve(H, ρ0, τlist, [], [c(π/2)] ; e_ops=[quad(π/2)], ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true))
# lindP = mesolve(H, ρ0, τlist, [c(π/2)] ; e_ops=[quad(π/2)], progress_bar=Val(true))

# Simulation outputs

Pmean_list = real.(homP.expect[1, :, :])
Pmean = vec(sum(Pmean_list, dims=1)./Ntraj)
# linPmean = real.(lindP.expect[1, :])
quadth(τ) = -sin(Ω*τ)*exp(-τ/2)/√2
### Plots

path = "HQ/"
mkpath(path)

quadPplot = plot(τlist, Pmean, label=L"\mathbb{E}_\mathrm{trajs}[\langle\hat{x}_{\pi/2}(\tau)\rangle]", xlabel=L"$\tau = kt$ (unitless)", ylabel=L"Mean $\frac{\pi}{2}$-quadrature $\propto \mathrm{Im}\{\rho_{eg}\}$", legend=:topright)
plot!(quadPplot, τlist, quadth.(τlist), ls=:dash, lw=3, label=L"\langle\hat{x}_{\pi/2}(\tau)\rangle_\mathrm{th} = \sqrt{2}\mathrm{Im}\{\rho_{eg}\}")

idxs = round.(Int, range(1, Ntraj, 3))
for (j,i) in enumerate(idxs)
    plot!(quadPplot, τlist, Pmean_list[i, :], label="Trajectory n°$(i)", alpha=0.5, ls=:dot, zorders=j)
end

savefig(quadPplot, path*"Pquadrature.pdf")


ErrorQuadPplot = plot(τlist, abs.(quadth.(τlist)-Pmean_list[end, :]), label="One trajectory-Lindblad", xlabel=L"$\tau = kt$ (unitless)", ylabel="Absolute difference", legend=:topright, color=palette[3], alpha=0.8)
plot!(ErrorQuadPplot, τlist, abs.(quadth.(τlist)-Pmean), label="Averaged trajectories-Lindblad", color=palette[1])
savefig(ErrorQuadPplot, path*"PquadSingleError.pdf")