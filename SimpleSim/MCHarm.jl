ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)
PATH = "MCH/"
cd("/home/julio/Desktop/GitHub/Stage2026/SimpleSim")
mkpath(PATH)

Ω = 10 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the couplign constant k. Weak coupling implies Ω>>1.
dτ = 2π/(50*Ω) # By definition T0=2π/ω0 and to work unitless: kT0=2πk/ω0. We choose to take 50 points per oscillation.
τlist = 0:dτ:10

α0 = 2.0 # complex number to define the initial coherent state: a|α>=α|α>
Ncut = abs(α0)^2+4*abs(α0) # Upper bound μ+4σ of initial coherent state distribution.
prct = 50 # POURCENTAGE DE TRONCATURE EN PLUS
factor = 1.0+prct/100
Ncut = Int(ceil(Ncut*factor)) # We take +prct% of the maximal population and round it to the upper value

a = destroy(Ncut) # annhilation operator, the matrix is truncated at Ncut order
ψ0 = coherent(Ncut, α0)
H = Ω*a'*a # Hamiltonian of the system (unitless)

c_ops=[a]

Ntraj = 10000
sim = mcsolve(H, ψ0, τlist, c_ops, e_ops=[a'*a], ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true))
trajs = real.(sim.expect[1,:,:])
meantraj = sum(trajs, dims=1)./Ntraj

n_th(τ) = (abs(α0)^2)*exp(-τ) # evolution of the mean population (zero-temperature)

ekwg = Dict(:subplot=>Dict(:legend_hfactor=>1.25))
P = plot(legend=:topright, xlabel=L"$\tau = kt$ (unitless)", ylabel="Mean population", extra_kwargs=ekwg)
idxs = round.(Int, range(1, Ntraj, 3))
for i in idxs
    plot!(P, τlist, trajs[i, :], label="Trajectory n°$(i)", alpha=0.7, ls=:dot, lw=5)
end
plot!(P, τlist, meantraj[1,:], label=L"\mathbb{E}[\langle \hat{n}\,(\tau)\rangle]", lw=2)
plot!(P, τlist, n_th.(τlist), label=L"$\langle \hat{n}\,(\tau)\rangle$ (analytical)", lw=3, ls=:dash)
plot!(P, τlist, abs.(n_th.(τlist)-meantraj[1,:]), title=L"|\mathbb{E}[\langle \hat{n}\,(\tau)\rangle]-\langle \hat{n}\,(\tau)\rangle|", color=palette[6], subplot=2, inset=(1, bbox(0.69, 0.45, 0.31, 0.31)), label=false)
savefig(P, PATH*"Pop_$(Ntraj).pdf")

all_jumps = vcat(sim.col_times...) # Every jump times are aligned in a 1D vertical vector
Histo = histogram(all_jumps, bins=100, normalize=:pdf, label="Average jump times (normalized)", 
          xlabel=L"$\tau = kt$ (unitless)", ylabel="PDF", legend=:topright)

τmax = τlist[end]
# PDF analytique normalisée à 1 sur [0, τmax]
Pth(τ) = exp(-τ) / (1 - exp(-τmax))

plot!(Histo, τlist, Pth.(τlist), label="Analytical PDF", lw=3, ls=:dash)
savefig(Histo, PATH*"histo.pdf")