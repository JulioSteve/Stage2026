ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings, Printf, Random, Distributions, Base.Threads, CurveFit
theme(:dao)
palette = theme_palette(:dao)

println("----------START----------")

function N_BE(x)
    if x == -1
        return 0.0 # χ→+∞ ≡ T→0 ⇒ N_th→0
    else
        return 1/(exp(x)-1)
    end
end

    # Simulation Parameters
# ===========================================
α0 = 2.0 # Starting level for the initial coherent state
χ = -1
Nth = N_BE(χ) # mean thermal photon number in the cavity
ω_k = 100.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_k>>1)


Ncut = Nth + 4*√(Nth*(1+Nth)) + abs(α0)^2+4*abs(α0) # Dimension cutoff in the operators: μ+4σ of thermal state (steady state) + μ+4σ of initial coherent state

prct = 0.0 # POURCENTAGE DE TRONCATURE EN PLUS
factor = 1.0+prct/100.0
Ncut = Int(ceil(Ncut*factor)) # We take +prct% of the maximal population and round it to the upper value
# println("⚠ Nth≈$(Int(ceil(Nth))) => Ncut=$Ncut")

# Operators
a = destroy(Ncut)

c_decay = √(1.0+Nth)*a
c_excite = √(Nth)*a'
c_ops = [c_decay, c_excite]

ρ0 = coherent_dm(Ncut, α0)

n_op = a'*a
e_ops=[n_op]
H = ω_k*a'*a 

dt = 2π/(100.0*ω_k) #1/100 of the unitless "frequency"
tlist = 0.0:dt:10.0

function quad(θ::Float64)
    return (a*exp(-im*θ)+a'*exp(im*θ))/sqrt(2)
end

function J(A::QuantumObject)
    return spre(A)+spost(A')
end

function dw_L(variance::Float64, length::Int64)
    σ=√(variance)
    d = Normal(0.0, σ)
    return rand(d, length)
end

function single_sim(θ::Float64)
    DW=dw_L(dt, length(tlist))

    # ρ_L = Vector{QuantumObject}(undef, length(tlist))
    meanquad_L = zeros(Float64, length(tlist))

    ρ = ρ0 #initial state
    for i in eachindex(tlist)
        # ρ_L[i] = ρ #saving the state
        dρ = -im*commutator(H,ρ)*dt+(lindblad_dissipator(c_decay)+lindblad_dissipator(c_excite))*ρ*dt+J(c_decay*exp(-im*θ))*ρ*DW[i]
        ρ = ρ + dρ # Evolution

        ###
        ρ_n = ρ/tr(ρ)
        meanquad_L[i] = real(tr(ρ_n*quad(θ)))
    end
    return meanquad_L
end

function sim(Ntraj::Int64, θ::Float64)
    results = [zeros(length(tlist)) for _ in 1:Ntraj]
    Threads.@threads for k in 1:Ntraj
        results[k] = single_sim(θ)
    end
    return sum(results) ./ Ntraj
end

function plotting()
    thsim = mesolve(H, ρ0, tlist, c_ops, e_ops=[quad(0.0)], progress_bar=Val(false))

    thquadx = real.(thsim.expect[1, :])

    P1 = plot(tlist, thquadx, label=L"\left\langle \hat{x}_0(\tau)\right\rangle_\mathrm{th}", xlabel=L"$\tau=kt$ (unitless)", ylabel="0-Quadrature", ylims=(-3,3))

    Y = sim(100_000, 0.0)

    plot!(P1, tlist, Y, label=L"\mathbb{E}[\langle \hat{x}_0(\tau)\rangle]")

    savefig(P1, "pop.svg")
end

ran=1:2:100

Time_to_exec = zeros(length(ran))
sim(1, 0.0) #warmup
for (i,ntraj) in enumerate(ran)
    t_iter = @elapsed begin
        sim(ntraj, 0.0)
    end
    Time_to_exec[i] = t_iter
end

prob = CurveFitProblem(ran, Time_to_exec)
sol = solve(prob, LinearCurveFitAlgorithm())
slope = round(sol.u[1], sigdigits=4)
intercept = round(sol.u[2], sigdigits=4)

Ptest = plot(ran, Time_to_exec, title=L"%$slope x+%$intercept", xlabel=L"Number of trajectories ($x$)", ylabel="Time (s)", lw=3)
plot!(Ptest, ran, sol.(ran), lw=1)
savefig(Ptest, "Timetoexec.svg")


println("----------STOP---------")