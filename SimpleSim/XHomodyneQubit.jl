ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)


Ω = 10 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the couplign constant k. Weak coupling implies Ω>>1.
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
homX = smesolve(H, ρ0, τlist, [], [c(0)] ; e_ops=[quad(0)], ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true))
lindX = mesolve(H, ρ0, τlist, [c(0)] ; e_ops=[quad(0)], progress_bar=Val(true))

# Simulation outputs
Xmean_list = real.(homX.expect[1, :, :])
Xmean = vec(sum(Xmean_list, dims=1)./Ntraj) # SME solution
# linXmean = real.(lindX.expect[1, :]) # Lindblad-ME solution

quadth(τ) = cos(Ω*τ)*exp(-τ/2)/√2

path = "HQ/"
mkpath(path)

quadXplot = plot(τlist, Xmean, label=L"\mathbb{E}_\mathrm{trajs}[\langle\hat{x}_0(\tau)\rangle]", xlabel=L"$\tau = kt$ (unitless)", ylabel=L"Mean $0$-quadrature $\propto \mathrm{Re}\{\rho_{eg}\}$", legend=:topright)

plot!(quadXplot, τlist, quadth.(τlist), ls=:dash, lw=3, label=L"\langle\hat{x}_0(\tau)\rangle_\mathrm{th} = \sqrt{2}\mathrm{Re}\{\rho_{eg}\}")
idxs = round.(Int, range(1, Ntraj, 3))
for (j,i) in enumerate(idxs)
    plot!(quadXplot, τlist, Xmean_list[i, :], label="Trajectory n°$(i)", alpha=0.5, ls=:dot, zorders=j)
end
savefig(quadXplot, path*"Xquadrature.pdf")

ErrorQuadXplot = plot(τlist, abs.(quadth.(τlist)-Xmean_list[end, :]), label="One trajectory-Lindblad", xlabel=L"$\tau = kt$ (unitless)", ylabel="Absolute difference", legend=:topright, color=palette[3], alpha=0.8)
plot!(ErrorQuadXplot, τlist, abs.(quadth.(τlist)-Xmean), label="Averaged trajectories-Lindblad", color=palette[1])
savefig(ErrorQuadXplot, path*"XquadSingleError.pdf")