ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings
theme(:dao)
palette = theme_palette(:dao)


Ω = 100 # Ω=ω0/k the ratio between natural frequency of the system ω0 and the couplign constant k. Weak coupling implies Ω>>1.
dτ = 2π/(10*Ω) # By definition T0=2π/ω0 and to work unitless: kT0=2πk/ω0. We take 1/10 of this oscillation period as "infinitesimal" element. 
τlist = 0:dτ:10

α0 = 2.0 # complex number to define the initial coherent state: a|α>=α|α>
Ncut = abs(α0)^2+4*abs(α0) # Upper bound μ+4σ of initial coherent state distribution.
prct = 0 # POURCENTAGE DE TRONCATURE EN PLUS
factor = 1.0+prct/100
Ncut = Int(ceil(Ncut*factor)) # We take +prct% of the maximal population and round it to the upper value

a = destroy(Ncut) # annhilation operator, the matrix is truncated at Ncut order
c(θ) = a*exp(-im*θ)
quad(θ) = (a*exp(-im*θ) + a'*exp(im*θ))/√2
ρ0 = coherent_dm(Ncut, α0)
H = Ω*a'*a # Hamiltonian of the system (unitless)

meanpop(t) = (abs(α0)^2)*exp(-t)

Ntraj = 100
homX = smesolve(H, ρ0, τlist, [], [c(0)] ; e_ops=[quad(0)], ntraj=Ntraj, progress_bar=Val(true))
lindX = mesolve(H, ρ0, τlist, [c(0)] ; e_ops=[quad(0)], progress_bar=Val(true))

# homP = smesolve(H, ρ0, τlist, [], [c(π/2)] ; e_ops=[quad(π/2)], ntraj=Ntraj, progress_bar=Val(true))
# lindP = mesolve(H, ρ0, τlist, [c(π/2)] ; e_ops=[quad(π/2)], progress_bar=Val(true))

# Ysim = real.(sol_lindblad.expect[1,:])
# Yth = meanpop.(τlist)
# pop = plot(τlist, Ysim, lw=2, label="sim")
# plot!(pop, τlist, Yth, ls=:dash, lw=3, label="th")
# plot!(pop, τlist, abs.(Yth.-Ysim), color=palette[3], label="", lw=2;
#     inset=bbox(0.45, 0.2, 0.4, 0.4), subplot=2)
# savefig(pop, "pop.svg")

Xmean = real.(homX.expect[1, :])
linXmean = real.(lindX.expect[1, :])
# Pmean = homP.expect[1, :]

quadXplot = plot(τlist, Xmean)
plot!(quadXplot, τlist, linXmean, ls=:dash)
savefig(quadXplot, "test.svg")