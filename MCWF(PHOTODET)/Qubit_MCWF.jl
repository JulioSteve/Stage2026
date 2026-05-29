using Plots, LinearAlgebra, Statistics, LaTeXStrings, Format, PlotThemes
theme(:dark)

println("-------------------------------")
# -------- Paramètres --------------------------------------------------------
const dt = 0.01
const Tmax = 10.0
const Ntraj = 10_000

# -------- Constantes --------------------------------------------------------
const ω0 = 2π
const Γ = 0.1ω0#couplage
#const Ω = 0.3ω0

# -------- Opérateurs --------------------------------------------------------
const σm = [0 0 ; 1 0]
const σz = [1 0 ; 0 -1]

# -------- Hamiltonien et Etat de départ -------------------------------------
const ψ0::Vector{Complex{Float64}} = [1 ; 0] #start in excited state
const H = ω0/2*σz #+ ħ*Ω*(σm'+σm)
const Heff = H - 1im*Γ/2*(σm'*σm)

# -------- Code --------------------------------------------------------------
function singletraj(ψ::Vector{Complex{Float64}}, dt::Float64, Tmax::Float64)
    """
    Simule une trajectoire quantique unique selon la méthode de Monte Carlo Wave Function (MCWF) pour un système à deux niveaux soumis à une dissipation.

        Arguments:
        - `ψ`: état quantique initial du système (vecteur de dimension 2)
        - `dt`: pas de temps pour l'évolution continue
        - `Tmax`: temps total de simulation

        Return:
        - `JumpTimes`: vecteur contenant les instants de saut
        - `βlist`: vecteur contenant les valeurs de <σ+σ-> pour chaque instant
        - `Braketlist`: vecteur contenant les valeurs de <ψ|ψ> pour chaque instant

    """
    R::Float64 = rand()
    Tlist::Vector{Float64} = collect(0:dt:Tmax)
    nsteps = length(Tlist)
    JumpTimes::Vector{Float64} = Float64[]
    βlist = zeros(Float64, nsteps)
    Braketlist = zeros(Float64, nsteps)

    for k in 1:nsteps
        braket = real(ψ'*ψ) #calcul de <ψ|ψ>(t) pour la trajectoire en cours à l'instant t
        ψn = ψ/sqrt(braket) #normalisation de ψ pour suivre évolution normalisée
        βlist[k] = real(ψn'*σm'*σm*ψn) #calcul de <σ+σ-> à l'instant t pour la trajectoire en cours
        Braketlist[k] = braket

        # println("norme ≈ ", round(real(ψ'*ψ), digits=3), " R ≈ ", round(R, digits=2)," ⚠️ ", round(k*100/nsteps, digits=2), "%", " at to ", Tlist[k],"/",Tmax, " time units")
        
        if real(ψ'*ψ) < R #condition de saut quantique (mesure effective)
            push!(JumpTimes, Tlist[k])
            ψ = σm*ψ #application de l'opérateur de saut
            ψ = ψ/norm(ψ) #renormalisation après saut
            R = rand() #nouvelle valeur aléatoire pour la prochaine condition de saut
        end

        # Evolution continue selon l'équation de Schrödinger non linéaire (RK4)
        k1 = -im*Heff*ψ
        k2 = -im*Heff*(ψ+dt*k1/2)
        k3 = -im*Heff*(ψ+dt*k2/2)
        k4 = -im*Heff*(ψ+dt*k3)
        ψ = ψ + dt/6*(k1+2k2+2k3+k4)
    end

    if isempty(JumpTimes) #si aucun saut

        # println("⚠️ No jump ⚠️")

        return nothing, βlist, Braketlist #on retourne βlist et Braketlist pour pouvoir les moyenner même sans saut
    else
        return JumpTimes, βlist, Braketlist
    end
end

function multitraj(N::Int, ψ::Vector{Complex{Float64}}, dt::Float64, Tmax::Float64)
    """
    Simule N trajectoires quantiques selon la méthode de Monte Carlo Wave Function (MCWF) pour un système à deux niveaux soumis à une dissipation, et calcule les observables moyennées.

        Arguments:
        - `N`: nombre de trajectoires à simuler
        - `ψ`: état quantique initial du système (vecteur de dimension 2)
        - `dt`: pas de temps pour l'évolution continue
        - `Tmax`: temps total de simulation

        Return:
        - `JumpList`: vecteur contenant les instants de saut de toutes les trajectoires
        - `nojump`: nombre de trajectoires n'ayant subi aucun saut
        - `βmoy`: vecteur contenant les valeurs moyennées de <σ+σ-> pour chaque instant
        - `braketmoy`: vecteur contenant les valeurs moyennées de <ψ|ψ> pour chaque instant
    """
    JumpList::Vector{Float64} = Float64[]
    nojump::Int64 = 0
    nsteps = length(0:dt:Tmax)
    βsum = zeros(Float64, nsteps)
    braketsum = zeros(Float64, nsteps)

    for n in 1:N
        print("\r Calculating... ", round(n*100/N, digits=2),"%")
        flush(stdout)
        jumps, β, braket = singletraj(ψ, dt, Tmax)
        if isnothing(jumps) 
            nojump += 1
        else
            append!(JumpList, jumps)
        end
        βsum .+= β
        braketsum .+= braket
    end

    println("\n\n", N-nojump, " trajectories jumped, out of ", N)
    println("Jump rate: ", round((N-nojump)*100/N, digits=2), "%")
    return JumpList, nojump, βsum/N, braketsum/N
end

function dN(dt::Float64, Tmax::Float64, jumptime::Float64)
    """
    Fonction de test pour vérifier que la distribution des sauts quantiques suit bien une loi exponentielle. Retourne un vecteur de 0 partout sauf à l'instant du saut où il vaut 1.

        Arguments:
        - `dt`: pas de temps pour l'évolution continue
        - `Tmax`: temps total de simulation
        - `jumptime`: instant du saut quantique

        Return:
        - `Tlist`: vecteur contenant les instants de temps
        - `dNlist`: vecteur contenant les valeurs de dN/dt (0 partout sauf à jumptime où il vaut 1)
    """
    Tlist::Vector{Float64} = collect(0:dt:Tmax)
    return [t == jumptime ? 1.0 : 0.0 for t in Tlist]
end

# -------- Plots --------------------------------------------------------------
Tlist::Vector{Float64} = collect(0:dt:Tmax)
PATH = homedir()*"/Bureau/Github/Stage2026/MCWF/PLOTS/"
ext = ".png"
dim = (1920, 1080)

#plot style:
const style = (
    size=dim,
    titlefontsize=20,
    guidefontsize=20,
    tickfontsize=15,
    legendfontsize=30,
    lw=6,
    margin=20Plots.mm,
    top_margin=15Plots.mm
)

function plot_single(doplot::Bool, save::Bool)
    if doplot
        JT = 0.0
        βsingle = nothing
        braketsingle = nothing
        while JT < Tmax/4.0
            JTs, βsingle,braketsingle = singletraj(ψ0, dt, Tmax)
            isnothing(JTs) ? continue : JT=JTs[1]
        end
        PLOT_single = plot(Tlist, βsingle, title="\nSingle trajectory", xlabel="Time", ylabel="Quantum observables", label=L"\langle \sigma_-^\dagger \sigma_-\rangle_\psi (t)", legend=:outerright; style..., lw=10)
        plot!(Tlist, braketsingle, ls=:dash, label=L"\langle \tilde\psi \;\, | \tilde\psi \rangle (t)"; style..., lw=8)
        plot!(Tlist, exp.(-Γ.*Tlist), ls=:dot, label=L"e^{-\Gamma t}"; style...)
        if save
            savefig(PLOT_single, PATH*"SingleTrajPlot_N$(Ntraj)_dt$(dt)_T$(Tmax)"*ext)
        end
    end
end

function plotting(doplot::Bool, save::Bool)
    if doplot
        RES, miss, βmoy, braketmoy = multitraj(Ntraj, ψ0, dt, Tmax) 
        # ^ jumptimes of the (N-notconverged) trajectories, number of missed convergence and mean of sigma+sigma- over time
        
        PLOT_mean = plot(Tlist, βmoy, label=L"\mathbb{E} \left[\langle \sigma_-^\dagger \sigma_- \rangle_\psi (t)\right]", title="\n$(format(Ntraj-miss, commas=true))/$(format(Ntraj, commas=true)) trajectories averaged - dt=$dt", xlabel="Time", ylabel="Averaged Quantum Observables", legend=:outerright; style..., lw=10)
        plot!(Tlist, braketmoy, ls=:dash, label=L"\mathbb{E} \left[\langle \tilde\psi \;\, | \tilde\psi \;\rangle (t)\right]"; style...)
        hline!([1.0], ls=:dashdot, label=L"\mathbb{E} \left[\langle \psi | \psi \rangle \right]"; style...) 
        plot!(Tlist, exp.(-Γ.*Tlist), ls=:dot, label=L"e^{-\Gamma t}"; style...)

        P2 = histogram(RES, bins=200, normalize=:pdf, label=L"\mathbb{E} \left[ \frac{dN}{dt} \right]", xlabel="Time", ylabel="Quantum jump distribution", title="\nNormalized histogram of dN\n$(format(Ntraj-miss, commas=true))/$(format(Ntraj, commas=true)) jumps - dt=$dt", legend=:outerright; style..., lw=1)

        plot!(Tlist, Γ*exp.(-Γ.*Tlist), ls=:dot, label=L"\Gamma e^{-\Gamma t}"; style...)
        if save
            savefig(PLOT_mean, PATH*"Plot_N$(Ntraj)_dt$(dt)_T$(Tmax)"*ext)
            savefig(P2, PATH*"Histogram_N$(Ntraj)_dt$(dt)_T$(Tmax)"*ext)
        end
    end
    return nothing
end

plot_single(true, true)
plotting(true, true)