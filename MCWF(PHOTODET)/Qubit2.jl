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

function single_traj()
    ψ=ψ0 # Initialization 
    Jumplist = zeros(Int, length(τlist)) # List of integers where 1 stands for a jump at corresponding τ
    Tr = zeros(Real, length(τlist))
    ρee_sim = zeros(Real, length(τlist))
    R = rand()
    # if Ntraj==1 #demonstrative value for single graph
    #     R=0.05
    # end
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

function sim(Ntraj::Integer)
    Traveraged = zeros(Real, length(τlist))
    Jumpalltraj = zeros(Int, length(τlist))
    ρeeaveraged = zeros(Real, length(τlist))
    for _ in 1:Ntraj
        Res_single = single_traj()
        Traveraged += Res_single[1]
        ρeeaveraged += Res_single[2]
        Jumpalltraj += Res_single[3]
    end
    Traveraged = Traveraged./Ntraj
    ρeeaveraged = ρeeaveraged./Ntraj

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
        Pav = plot()
        
        plot!(Pav, τlist, ρeeaveraged, label=L"\mathbb{E}\left[\rho_{ee}(\tau)\right]=\mathbb{E}[\langle\sigma_+\sigma_-\rangle_{\psi}(\tau)]", color=palette[1], lw=5)
        plot!(Pav, τlist, Traveraged, label=L"\mathbb{E}\left[\mathrm{Tr}\;[\breve{\rho}\,(\tau)]\right]=\mathbb{E}\left[\langle\breve{\psi}\;\, |\breve{\psi}\,\rangle(\tau)\right]", ls=:dash, color=palette[2], lw=3)
        plot!(Pav, τlist, ρee.(τlist), label=L"\rho_{ee}^\mathrm{th}(\tau)", ls=:dot, color=palette[3], lw=2)
        savefig(Pav, "$(Ntraj)trajs.svg")
        
        Histo = histogram(Jumpalltraj, bins=200, normalize=:pdf, legend=:outerright, xlabel=L"$\tau = kt$ (unitless)", ylabel="Quantum jump distribution", extra_kwargs=ekw, label="")
        savefig(Histo, "histo.svg")
    end

end

####################################

function timetest()
    Times = []
    for ntraj in 10:1:100
        titer = @elapsed begin
            sim(ntraj);
        end
        push!(Times, titer)
    end

    times = plot(10:100, Times)
    savefig(times, "Timetoexec.svg")
end

# timetest()

####################################

Ntraj = 100_000
RES = sim(Ntraj)
plotting(Ntraj, RES)