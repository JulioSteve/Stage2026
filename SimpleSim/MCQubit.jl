ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)

Ω = 10 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the couplign constant k. Weak coupling implies Ω>>1.
dτ = 2π/(50*Ω) # By definition T0=2π/ω0 and to work unitless: kT0=2πk/ω0. We choose to take 50 points per oscillation.
τlist = 0:dτ:10

ψ0 = (fock(2,0)+fock(2,1))./√2
σm = sigmam()
σz = sigmaz()
Pe = fock(2,0)*fock(2,0)' # projector |e><e|. <Pe>=ρee the excited population.
H = Ω*σz/2 # unitless Hamiltonian

c_ops = [σm]

Ntraj = 10
sim = mcsolve(H, ψ0, τlist, c_ops, e_ops=[Pe], ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true))


