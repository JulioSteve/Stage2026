using Plots, LinearAlgebra, Statistics, LaTeXStrings

println("-------------------------------")
# -------- Paramètres --------------------------------------------------------
const dt = 0.001
const Tmax = 10
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
    βlist = zeros(nsteps)

    for k in 1:nsteps
        βlist[k] = real(ψ'*σm'*σm*ψ)
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
        ψ = ψ + dt/6*(k1+2k2+2k3+k4) #RK4 Algorithm
    end

    if isnothing(res)
        # println("⚠️ The algorithm didn't converge ⚠️")
        return nothing, βlist
    else
        return Tlist[res], βlist
    end
end

function multitraj(N, ψ, dt, Tmax)
    JumpList = Float64[]
    notconv = 0
    nsteps = length(0:dt:Tmax)
    βsum = zeros(nsteps)

    for n in 1:N
        print("\r Calculating... ", round(n*100/N, digits=2),"%")
        flush(stdout)
        res, β = singletraj(ψ, dt, Tmax)
        if isnothing(res) 
            notconv +=1 #comptabilise les non-convergences
        else
            push!(JumpList, res)
        end
        βsum .+= β
    end
    println("\n\n",N-notconv, " trajectories jumped, out of ", N)
    println("Jump rate: ", round((N-notconv)*100/N, digits=2),"%")
    return JumpList, notconv, βsum/N
end

function dN(dt, Tmax, jumptime)
    Tlist = 0:dt:Tmax
    return [t == jumptime ? 1.0 : 0.0 for t in Tlist]
end

# -------- Plots --------------------------------------------------------------
Tlist = 0:dt:Tmax
PATH = homedir()*"/Bureau/Github/Stage2026/MCWF/ByHands/PLOTS/"
ext = ".png"
dim = (1920, 1080)
#plot style
const style = (
    size=dim,
    titlefontsize=20,
    guidefontsize=20,
    tickfontsize=15,
    legendfontsize=16,
    lw=4,
    margin=10Plots.mm
)
function plot_single(doplot, save)
    if doplot
        _, βsingle = singletraj(ψ0, dt, Tmax)
        PLOT_single = plot(Tlist, βsingle, title="Single trajectory", label=false, xlabel="Time", ylabel=L"\sigma_-^\dagger \sigma_- (t)"; style...)
        if save
            savefig(PLOT_single, PATH*"SingleTrajPlot"*ext)
        end
    end
end

function plotting(doplot, save)
    if doplot
        RES, miss, βmoy = multitraj(Ntraj, ψ0, dt, Tmax) 
        # ^ jumptimes of the (N-notconverged) trajectories, number of missed convergence and mean of sigma+sigma- over time
        
        PLOT_mean = plot(Tlist, βmoy, label=L"\left< \beta (t) \right>", title="$(Ntraj-miss)/$Ntraj trajectories averaged - dt=$dt", xlabel="Time", ylabel=L"<\sigma_-^\dagger \sigma_- (t)>"; style...)
        
        P2 = histogram(RES, bins=200, normalize=:pdf, label="Quantum jump distribution", xlabel="Time", ylabel="Normalized count", title="Histogram of dN\n$(Ntraj-miss)/$Ntraj jumps - dt=$dt"; style..., lw=1)
        
        if save
            savefig(PLOT_mean, PATH*"$(Ntraj)TrajsPlot_dt$dt"*ext)
            savefig(P2, PATH*"$(Ntraj)HistodN_dt$dt"*ext)
        end
    end
    return nothing
end

plot_single(false, false)
plotting(true, true)