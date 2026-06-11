ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)
PATH = "MCQ/"
mkpath(PATH)

Ω = 10 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the couplign constant k. Weak coupling implies Ω>>1.
dτ = 2π/(50*Ω) # By definition T0=2π/ω0 and to work unitless: kT0=2πk/ω0. We choose to take 50 points per oscillation.
τlist = 0:dτ:10

ψ0 = (fock(2,0)+fock(2,1))./√2
# ψ0 = fock(2,0)
σm = sigmam()
σz = sigmaz()
Pe = fock(2,0)*fock(2,0)' # projector |e><e|. Note that <Pe>=ρee is the excited population.
H = Ω*σz/2 # unitless Hamiltonian

c_ops = [σm]

Ntraj = 5000
sim = mcsolve(H, ψ0, τlist, c_ops, e_ops=[Pe], ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true))
trajs = real.(sim.expect[1,:,:])
meantraj = sum(trajs, dims=1)./Ntraj

Pe0 = real(ψ0'*Pe*ψ0) # initial mean excited population
Pe_th(τ) = Pe0*exp(-τ) # evolution of the mean excited population

ekwg = Dict(:subplot=>Dict(:legend_hfactor=>1.25))
P = plot(legend=:topright, xlabel=L"$\tau = kt$ (unitless)", ylabel=L"Excited population $\rho_{ee}$", extra_kwargs=ekwg)
idxs = round.(Int, range(1, Ntraj, 3))
for i in idxs
    plot!(P, τlist, trajs[i, :], label="Trajectory n°$(i)", alpha=0.7, ls=:dot)
end
plot!(P, τlist, meantraj[1,:], label=L"\mathbb{E}[\rho_{ee}(\tau)]", lw=2)
plot!(P, τlist, Pe_th.(τlist), label=L"$\rho_{ee}(\tau)$ (analytical)", lw=3, ls=:dash)
plot!(P, τlist[1:(length(τlist)÷2)], abs.(Pe_th.(τlist)-meantraj[1,:])[1:(length(τlist)÷2)], title=L"|\mathbb{E}[\rho_{ee}(\tau)]-\rho_{ee}(\tau)|", color=palette[6], subplot=2, inset=(1, bbox(0.7, 0.55, 0.31, 0.31)), label=false)
savefig(P, PATH*"Pop_$(Ntraj).svg")

all_jumps = vcat(sim.col_times...) # Every jump times are aligned in a 1D vertical vector
Histo = histogram(all_jumps, bins=100, normalize=:pdf, label="Jump times", 
          xlabel=L"$\tau = kt$ (unitless)", ylabel="PDF", legend=:topright)
plot!(Histo, τlist, exp.(-τlist), label="Analytical???", lw=3, ls=:dash)
savefig(Histo, PATH*"histo.svg")