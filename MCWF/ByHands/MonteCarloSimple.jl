using Plots, LinearAlgebra, Statistics, LaTeXStrings

println("-------------------------------")
# -------- Paramètres --------------------------------------------------------
const dt = 1e-4
const Tmax = 5
const Ntraj = 10_000

# -------- Constantes --------------------------------------------------------
const ω0 = 2π
const Γ = 0.1ω0#couplage
#const Ω = 0.3ω0

# -------- Opérateurs --------------------------------------------------------
const σm = [0 0 ; 1 0]
const σz = [1 0 ; 0 -1]

# -------- Hamiltonien et Etat de départ -------------------------------------
const ψ0 = [1 ; 0] #start in excited state
const H = ω0/2*σz #+ ħ*Ω*(σm'+σm)
const Heff = H - 1im*Γ/2*(σm'*σm)

# -------- Code --------------------------------------------------------------
function singletraj(ψ, dt, Tmax)
    R = rand()
    # println("Nombre aléatoire: ", R)
    Tlist = 0:dt:Tmax
    nsteps = length(Tlist)
    res = nothing

    for k in 1:nsteps
        # println("norme ≈ ", round(real(ψ'*ψ), digits=3), " R ≈ ", round(R, digits=2)," ⚠️ ", round(k*100/nsteps, digits=2), "%", " at to ", Tlist[k],"/",Tmax, " time units")
        if real(ψ'*ψ) < R
            res = k
            # println("Trouvée au pas ", k)
            # println("\nJump at ", Tlist[res],"/", Tmax," time units")
            break
        end

        k1 = -im*Heff*ψ
        k2 = -im*Heff*(ψ+dt*k1/2)
        k3 = -im*Heff*(ψ+dt*k2/2)
        k4 = -im*Heff*(ψ+dt*k3)
        ψ = ψ + dt/6*(k1+2k2+2k3+k4)
    end

    if isnothing(res)
        # println("⚠️ The algorithm didn't converge ⚠️")
        return nothing
    else
        return Tlist[res]
    end
end

function multitraj(N, ψ, dt, Tmax)
    JumpList = Float64[]
    notconv = 0

    for n in 1:N
        print("\r Calculating... ", round(n*100/N, digits=2),"%")
        flush(stdout)
        res = singletraj(ψ, dt, Tmax)
        if isnothing(res) 
            notconv +=1 #comptabilise les non-convergences
        else
            push!(JumpList, res)
        end
    end
    println("\n\n",N-notconv, " trajectories converged to a solution, out of ", N)
    println("Success rate: ", round((N-notconv)*100/N, digits=2),"%")
    return JumpList, notconv
end

function Representation(dt, Tmax, jumptime)
    Tlist = 0:dt:Tmax
    return [t < jumptime ? 1.0 : 0.0 for t in Tlist]
end

function dN(dt, Tmax, jumptime)
    Tlist = 0:dt:Tmax
    return [t == jumptime ? 1.0 : 0.0 for t in Tlist]
end

# -------- Plots --------------------------------------------------------------
RES, miss = multitraj(Ntraj, ψ0, dt, Tmax) # jumptimes of the (N-notconverged) trajectories
Tlist = 0:dt:Tmax

function trajplot(doplot, save)
    if doplot
        Y = [Representation(dt, Tmax, jump) for jump in RES] # converts jumptime of one trajectory to values of excited population (1 -> 0)
        PATH = "MCWF/"
        rand_idx = rand(1:length(RES))
        PLOT_single = plot(Tlist, Y[rand_idx], lw=2, title="Single trajectory - dt=$dt", label=false, xlabel="Time", ylabel="Excited population", tickfontsize=12)
        yticks!([0,1], [L"|g\rangle", L"|e\rangle"])

        PLOT_mean = plot(Tlist, mean(Y), lw=2, label=L"\left< \rho_{ee} \right>", tickfontsize=12, title="$(Ntraj-miss)/$Ntraj trajectories averaged - dt=$dt", xlabel="Time", ylabel="Excited population")
        plot!(Tlist, exp.(-Γ*collect(Tlist)), lw=2, ls=:dot, label=L"e^{-\Gamma t}")
        yticks!([0,1], [L"|g\rangle", L"|e\rangle"])

        if save
            savefig(PLOT_single, PATH*"SingleTrajPlot.pdf")
            savefig(PLOT_mean, PATH*"$Ntraj"*"TrajsPlot.pdf")
        end
    end
    return nothing
end
trajplot(false, false)

function dNplot(doplot, save)
    if doplot
        P2 = histogram(RES, bins=100, normalize=:pdf, label="Quantum jump distribution", xlabel="Time", ylabel="Normalized count", tickfontsize=12)
        plot!(Tlist, Γ*exp.(-Γ*collect(Tlist)), lw=2, ls=:dash, label=L"\Gamma e^{-\Gamma t}", title="Histogram of dN\n$(Ntraj-miss)/$Ntraj trajectories averaged - dt=$dt")
        if save
            savefig(P2, PATH*"HistodN.pdf")
        end
    end
end
