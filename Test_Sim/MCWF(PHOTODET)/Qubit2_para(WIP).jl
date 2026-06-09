ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings, Printf, Base.Threads
theme(:dao)
palette = theme_palette(:dao)

#######################

χ=Inf
Ω=100.0 # Ω=ω0/k where ω0 is the natural frequency of the system and k is the coupling constant: Ω>>1 <=> k<<ω0.
dτ=2π/(50.0*Ω) # T0 = 2π/ω0 unitless: kT0, we choose dτ=kT0/10.
τlist = 0.0:dτ:10.0

#######################
BE(X) = 1.0/(exp(X)-1.0)
""" BE:
Distribution de Bose-Einstein avec la notation χ=ħω/(k_B T).
Le cas T->0 est équivalent à χ=+∞ et donc il faut évaluer la fonction en écrivant "Inf" (géré par Julia proprement).
"""
Nth=BE(χ)
σz = sigmaz() ; σp = sigmap()
H = Ω*σz/2
Emi = (1+Nth)/2*σp*σp'
Abs = Nth/2*σp'*σp
Heff = H-im*(Emi+Abs)
ketE = fock(2,0) # |e>
ψ0 = ketE # Initial state

ρee0 = abs(ψ0'*ketE)^2 # <e|ρ(t=0)|e>
ρee∞ = Nth/(1+2Nth) # Thermalization 
ρee(τ) = (ρee0-ρee∞)*exp(-(1+2Nth)*τ) + ρee∞ # analytical expression of ρee(t)=<e|ρ(t)|e>

function RK4(ψ::QuantumObject, Heff::QuantumObject)
    f(ϕ) = -im*Heff*ϕ
    k1 = f(ψ)
    k2 = f(ψ+dτ/2.0*k1)
    k3 = f(ψ+dτ/2.0*k2)
    k4 = f(ψ+dτ*k3)
    return ψ + (k1+2.0*k2+2.0*k3+k4)*dτ/6.0
end

function single_traj(; forceR=nothing)
    ψ=ψ0 # Initialization 
    Jumplist = zeros(Int, length(τlist)) # List of integers where 1 stands for a jump at corresponding τ
    Tr = zeros(Real, length(τlist))
    ρee_sim = zeros(Real, length(τlist))
    R = isnothing(forceR) ? rand() : forceR
    for i in 1:length(τlist)
        ψn = ψ/√(real(ψ'*ψ))
        Tr[i] = real(ψ'*ψ)
        ρee_sim[i] = abs(ψn'*ketE)^2
        if real(ψ'*ψ)<R
            Jumplist[i] = 1
            ψ = σp'*ψ
            ψ = ψ/√(real(ψ'*ψ))
        end
        ψ = RK4(ψ, Heff)
    end
    return Tr, ρee_sim, Jumplist
end

function sim(Ntraj::Integer; forceR=nothing)
    n = length(τlist)

    nt = nthreads() # one list for each thread:
    Tr_parts = [zeros(Float64, n) for _ in 1:nt]
    ρee_parts = [zeros(Float64, n) for _ in 1:nt]
    Jump_parts = [zeros(Int, n) for _ in 1:nt]

    @threads :static for _ in 1:Ntraj
        tid = threadid()
        Tr_t, ρee_t, Jump_t = single_traj()
        Tr_parts[tid] .+= Tr_t
        ρee_parts[tid] .+= ρee_t
        Jump_parts[tid] .+= Jump_t
    end

    # Reducing to average accross all threads
    Traveraged   = sum(Tr_parts)   ./ Ntraj
    ρeeaveraged  = sum(ρee_parts)  ./ Ntraj
    Jumpalltraj  = sum(Jump_parts)

    return Traveraged, ρeeaveraged, Jumpalltraj
end

function plotting(Ntraj::Integer, simres)
    Traveraged, ρeeaveraged, Jumpalltraj = simres
    ekw = Dict(:subplot => Dict(:legend_hfactor => 1.5))

    if Ntraj==1
        Psingle = plot()
        
        plot!(Psingle, τlist, ρeeaveraged, label=L"\rho_{ee}(\tau)=\langle\sigma_+\sigma_-\rangle_{\psi}(\tau)", color=palette[1], lw=5)

        plot!(Psingle, τlist, Traveraged, label=L"\mathrm{Tr}\;[\breve{\rho}\,(\tau)]=\langle\breve{\psi}\;\, |\breve{\psi}\,\rangle(\tau)", ls=:dash, color=palette[2], lw=3)

        plot!(Psingle, τlist, ρee.(τlist), label=L"\rho_{ee}^\mathrm{th}(\tau)", ls=:dot, color=palette[3], lw=2)

        savefig(Psingle, "single.svg")
    else
        Pav = plot(legend=:outerright, xlabel=L"$\tau = kt$ (unitless)", ylabel="Expectation value of averaged observables", extra_kwargs=ekw)
        
        plot!(Pav, τlist, ρeeaveraged, label=L"\mathbb{E}\left[\rho_{ee}(\tau)\right]=\mathbb{E}[\langle\sigma_+\sigma_-\rangle_{\psi}(\tau)]", color=palette[1], lw=5)

        plot!(Pav, τlist, Traveraged, label=L"\mathbb{E}\left[\mathrm{Tr}\;[\breve{\rho}\,(\tau)]\right]=\mathbb{E}\left[\langle\breve{\psi}\;\, |\breve{\psi}\,\rangle(\tau)\right]", ls=:dash, color=palette[2], lw=3)

        plot!(Pav, τlist, ρee.(τlist), label=L"\rho_{ee}^\mathrm{th}(\tau)", ls=:dot, color=palette[3], lw=2)
        
        savefig(Pav, "$(Ntraj)trajs.svg")
    end

end

####################################

Ntraj = 10_000
# RES = sim(1, forceR=0.05) # good visual for 1 single traj
RES = sim(Ntraj, forceR=nothing)
plotting(Ntraj, RES)

####################################