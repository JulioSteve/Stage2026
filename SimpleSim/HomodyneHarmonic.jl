ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)
cd("/home/julio/Desktop/GitHub/Stage2026/SimpleSim/")

Ω = 10 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the couplign constant k. Weak coupling implies Ω>>1.
dτ = 2π/(50*Ω)/10 # By definition T0=2π/ω0 and to work unitless: kT0=2πk/ω0. We choose to take 500 points per oscillation.
τlist = 0:dτ:10

α0 = 2.0 # complex number to define the initial coherent state: a|α>=α|α>
Ncut = abs(α0)^2+4*abs(α0) # Upper bound μ+4σ of initial coherent state distribution.
prct = 200 # POURCENTAGE DE TRONCATURE EN PLUS
factor = 1.0+prct/100
Ncut = Int(ceil(Ncut*factor)) # We take +prct% of the maximal population and round it to the upper value

a = destroy(Ncut) # annhilation operator, the matrix is truncated at Ncut order
c(θ) = a*exp(-im*θ)
quad(θ) = (a*exp(-im*θ) + a'*exp(im*θ))/√2
ρ0 = coherent_dm(Ncut, α0)
H = Ω*a'*a # Hamiltonian of the system (unitless)

meanpop(t) = (abs(α0)^2)*exp(-t)

# Simulations
Ntraj = 100
homX = smesolve(H, ρ0, τlist, [], [c(0)] ; e_ops=[quad(0)], ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true))
# lindX = mesolve(H, ρ0, τlist, [c(0)] ; e_ops=[quad(0)], progress_bar=Val(true))

homP = smesolve(H, ρ0, τlist, [], [c(π/2)] ; e_ops=[quad(π/2)], ntraj=Ntraj, progress_bar=Val(true), keep_runs_results=Val(true))
# lindP = mesolve(H, ρ0, τlist, [c(π/2)] ; e_ops=[quad(π/2)], progress_bar=Val(true))

# Simulation outputs
Xmean_list = real.(homX.expect[1, :, :])
Xmean = vec(sum(Xmean_list, dims=1)./Ntraj)
# linXmean = real.(lindX.expect[1, :])

Pmean_list = real.(homP.expect[1, :, :])
Pmean = vec(sum(Pmean_list, dims=1)./Ntraj)
# linPmean = real.(lindP.expect[1, :])

quadth(τ,θ) = √2*abs(α0)*exp(-τ/2)*cos(Ω*τ+θ)

function plotting()
    path = "HH/"
    mkpath(path)
    quadXplot = plot(τlist, Xmean, label=L"\mathbb{E}[\langle\hat{a}_0(τ)\rangle]", xlabel=L"$\tau = kt$ (unitless)", ylabel="Mean value of quantum operators", legend=:topright)
    plot!(quadXplot, τlist, quadth.(τlist,0), ls=:dash, lw=3, label=L"\langle\hat{a}_0(τ)\rangle_\mathrm{th}")
    savefig(quadXplot, path*"Xquadrature.pdf")

    ErrorQuadXplot = plot(τlist, abs.(quadth.(τlist,0)-Xmean_list[end, :]), label="One trajectory-Lindblad", xlabel=L"$\tau = kt$ (unitless)", ylabel="Absolute difference", legend=:topright, color=palette[3], alpha=0.8)
    plot!(ErrorQuadXplot, τlist, abs.(quadth.(τlist,0)-Xmean), label="Averaged trajectories-Lindblad", color=palette[1])
    savefig(ErrorQuadXplot, path*"XquadSingleError.pdf")

    quadPplot = plot(τlist, Pmean, label=L"\mathbb{E}[\langle\hat{a}_{\pi/2} (τ)\rangle]", xlabel=L"$\tau = kt$ (unitless)", ylabel="Mean value of quantum operators", legend=:topright)
    plot!(quadPplot, τlist, quadth.(τlist,π/2), ls=:dash, lw=3, label=L"\langle\hat{a}_{\pi/2} (τ)\rangle_\mathrm{th}")
    savefig(quadPplot, path*"Pquadrature.pdf")

    ErrorQuadPplot = plot(τlist, abs.(quadth.(τlist,π/2)-Pmean_list[end, :]), label="One trajectory-Lindblad", xlabel=L"$\tau = kt$ (unitless)", ylabel="Absolute difference", legend=:topright, color=palette[3], alpha=0.8)
    plot!(ErrorQuadPplot, τlist, abs.(quadth.(τlist,π/2)-Pmean), label="Averaged trajectories-Lindblad", color=palette[1])
    savefig(ErrorQuadPplot, path*"PquadSingleError.pdf")

    open(path*"settings.txt", "w") do file
        println(file, "Ntraj = $(Ntraj)")
        println(file, "ω = $Ω×k")
        println(file, "Number of points per oscillation: $( round(2π/(Ω*dτ), sigdigits=2) )")
    end
end

plotting()