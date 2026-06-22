ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)
cd("/home/julio/Desktop/GitHub/Stage2026/SimpleSim")
PATH = "MCQ/"
mkpath(PATH)

ψsup = (fock(2,0)+fock(2,1))./√2
ψe = fock(2,0)

function simulation(ψ0)
    Ω = 10 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the couplign constant k. Weak coupling implies Ω>>1.
    dτ = 2π/(50*Ω) # By definition T0=2π/ω0 and to work unitless: kT0=2πk/ω0. We choose to take 50 points per oscillation.
    τlist = 0:dτ:10

    σm = sigmam()
    σz = sigmaz()
    Pe = fock(2,0)*fock(2,0)' # projector |e><e|. Note that <Pe>=ρee is the excited population.
    H = Ω*σz/2 # unitless Hamiltonian

    c_ops = [σm]

    Ntraj = 10_000
    sim = mcsolve(H, ψ0, τlist, c_ops, ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true), reltol=1e-15, abstol=1e-15, normalize_states=Val(false))

    # Graphe évolution norme de l'état
    ridx = findfirst(t -> !isempty(t) && t[1] >= 3.0 && t[1] <= 3.5, sim.col_times)
    ψr = sim.states[ridx, :]

    Nall = norm.(sim.states)
    MeanN = sum(Nall, dims=1)/Ntraj

    α = fock(2,0)'*ψ0
    β = fock(2,1)'*ψ0
    N0_th(t) = √(exp(-t)*abs(α)^2+abs(β)^2) 

    ekwg = Dict(:subplot=>Dict(:legend_hfactor=>1.25))
    P = plot(ylabel="Norme of states", xlabel=L"$\tau = kt$ (unitless)", legend=:right, extra_kwargs=ekwg)
    plot!(P, τlist, norm.(ψr), label=L"$\sqrt{\langle \psi | \psi\rangle}$ (single traj.)", lw=3)
    plot!(P, τlist,  MeanN[1,:], label=L"\mathbb{E}\left[ \sqrt{\langle \psi | \psi\rangle} \right]", lw=3, zorders=1)
    plot!(P, τlist, N0_th.(τlist), label=L"$\sqrt{\langle \psi_0\, | \psi_0\rangle}$ (analytical)", lw=5, ls=:dot)

    if ψ0 == ψsup
        savefig(P, PATH*"NormeSuperposition.pdf")
    elseif ψ0 == ψe
        savefig(P, PATH*"NormeE.pdf")
    end
end

simulation(ψsup)
simulation(ψe)