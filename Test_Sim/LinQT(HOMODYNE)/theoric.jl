ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings, Printf, Random, Base.Threads
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
    return √(variance).*randn(length)
end

function single_sim(θ::Float64, q::QuantumObject)
    DW=dw_L(dt, length(tlist))

    # ρ_L = Vector{QuantumObject}(undef, length(tlist))
    meanquad_L = zeros(Float64, length(tlist))

    ρ = ρ0 #initial state
    cθ = c_decay*exp(-im*θ)
    for i in eachindex(tlist)
        ρ = ρ/real(tr(ρ)) # Renormalization (else trace explodes)
        meanquad_L[i] = real(tr(ρ*q))

        # ρ_L[i] = ρ #saving the state
        Milstein = (J(cθ)*ρ)*(J(cθ)*eye(Ncut))*(DW[i]^2-dt)/2
        dρ = -im*commutator(H,ρ)*dt+(lindblad_dissipator(c_decay)+lindblad_dissipator(c_excite))*ρ*dt+J(cθ)*ρ*DW[i] + Milstein

        ρ = ρ + dρ # Evolution
    end
    return meanquad_L
end

function sim(Ntraj::Int64, θ::Float64)
    q = quad(θ)
    results = [zeros(length(tlist)) for _ in 1:Ntraj]
    Threads.@threads for k in 1:Ntraj
        results[k] = single_sim(θ, q)
    end
    moy = sum(results) ./ Ntraj
    moy_quad = sum(r .^2 for r in results)./Ntraj
    standdev = .√(moy_quad .-(moy .^2))
    return moy, standdev
end

function plotting(xlims::Tuple{Float64, Float64}, ylims::Tuple{Float64, Float64})
    thsim = mesolve(H, ρ0, tlist, c_ops, e_ops=[quad(0.0)], progress_bar=Val(false))

    thquadx = real.(thsim.expect[1, :])

    P1 = plot(tlist, thquadx, label=L"\left\langle \hat{x}_0(\tau)\right\rangle_\mathrm{th}", xlabel=L"$\tau=kt$ (unitless)", ylabel="0-Quadrature", ylims=ylims, xlims=xlims)

    # plot!(P1, tlist, QuadX, label=L"\mathbb{E}[\langle \hat{x}_0(\tau)\rangle]", ribbon=StdX)
    savefig(P1, "pop_th.svg")
end

mesolve(H, ρ0, tlist[1:5], c_ops, e_ops=[quad(0.0)], progress_bar=Val(false)) #warmup
single_sim(0.0, quad(0.0)) #warmup

plotting((0.0,1.0), (-3.0,3.0))
println("----------STOP---------")