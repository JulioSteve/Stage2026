ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings, Printf, Base.Threads
theme(:dao)
palette = theme_palette(:dao)

N_BE(x) = 1.0/(exp(x)-1.0)

println("---------------------")

    # Simulation Parameters
# ===========================================
α0 = 2.0 # Starting level for the initial coherent state
χ = Inf
Nth = N_BE(χ) # mean thermal photon number in the cavity 
tlist = range(0,10,1000)
ω_k = 100 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_k>>1)


Ncut = Nth + 4*√(Nth*(1+Nth)) + abs(α0)^2+4*abs(α0) # Dimension cutoff in the operators: μ+4σ of thermal state (steady state) + μ+4σ of initial coherent state

prct = 50 # POURCENTAGE DE TRONCATURE EN PLUS
factor = 1.0+prct/100
Ncut = Int(ceil(Ncut*factor)) # We take +prct% of the maximal population and round it to the upper value
# println("⚠ Nth≈$(Int(ceil(Nth))) => Ncut=$Ncut")

# Operators
a = destroy(Ncut)

c_decay = √(1+Nth)*a
c_excite = √(Nth)*a'
c_ops = [c_decay, c_excite]

ψ0 = coherent(Ncut, α0)

n_op = a'*a
e_ops=[n_op]
H = ω_k*a'*a 

sol = mesolve(H, ψ0, tlist, c_ops; e_ops=e_ops, progress_bar=Val(true)) # solver
ρf = sol.states[end]
pop_sim = real.(diag(ρf))

function pop_th(level,χ)
    if χ==-1 && level==0
        return 1.0
    elseif χ==-1 && level!=0
        return 0.0
    else
        return (1-exp(-χ))*exp(-level*χ)
    end
end

Levels = 0:Ncut-1
pop_thv = pop_th.(Levels,χ)
yticksval = range(minimum(pop_sim),maximum(pop_sim), 7)
ytickslab = ["$(round(v*100, digits=2))%" for v in yticksval]
ytickslab = ["$(round(v, sigdigits=1))" for v in yticksval]

pop = plot(Levels, pop_sim, label="pop sim", lw=5, title="Population of each level after the evolution\n"*L"\sum pop_{th}\approx %$(round(sum(pop_thv), digits=5))\qquad \sum pop_{sim}\approx %$(round(sum(pop_sim), digits=10))", xlabel="Levels", ylabel="Population", yticks=(yticksval, ytickslab), legend=:bottomleft)
plot!(pop, Levels, pop_thv, label="pop th", ls=:dash)
plot!(pop, Levels, abs.(pop_thv.-pop_sim), color=palette[3], label="", title="Absolute error", lw=2 ; inset=bbox(0.5, 0.2, 0.4, 0.4), subplot=2)
savefig(pop, "populations$(ω_k).svg")

n_anal(t) = α0^2*exp(-t)+Nth*(1-exp(-t)) # analytical solution
N = n_anal.(tlist)
Y = real.(sol.expect[1,:])
Y2 = abs.(N.-Y)

p = plot(title=L"ħω=%$(χ)\times k_{B}T \qquad N_{cut}=%$(Ncut)"*"(+$(prct)%)", legend=:outerright)
# Graph of the population
plot!(p, tlist, Y,label="Simulation",lw=3, color=palette[1]);
plot!(p, tlist, N, lw=2, label="Analytical model", color=palette[2], ls=:dash);
xlabel!(p, L"$\tau=k t$ (unitless)")
ylabel!(p, L"\langle \hat{n}(\tau)\rangle")

# Error graph (inset inside the previous graph)
yticks_val = range(minimum(Y2),maximum(Y2),5)
yticks_label = [@sprintf("%.1e", v) for v in yticks_val]
yticks = (yticks_val, yticks_label)

plot!(p, tlist, Y2, color=palette[3], label="", yticks=yticks, title="Absolute error", lw=2 ; inset=bbox(0.3, 0.4, 0.4, 0.4), subplot=2)
savefig(p, "$(prct)%$(ω_k).svg")

P_coh(n,a) = exp(-abs(a)^2)*abs(a)^(2n)/factorial(big(n))
yticksv = 10.0 .^(.-(1:6))

distrib = plot(Levels, pop_thv, title="Population distribution (log scale on the y-scale)", xlabel=L"n", ylabel=L"\mathbb{P}(n)", label="Thermal", legend=:topright, z_order=:front, lw=3, yscale=:log10, yticks=:native, ylims=(1e-7,0.1), ls=:dash)
plot!(distrib, Levels, pop_sim, label="Simulation population", z_order=:back, lw=3)
savefig(distrib, "distrib$(ω_k).svg")

println("---------------------")