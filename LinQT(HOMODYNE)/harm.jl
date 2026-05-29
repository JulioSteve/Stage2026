using Plots, LaTeXStrings, QuantumToolbox
theme(:dao)
palette = theme_palette(:dao)

function N_BE(x)
    if x == 0
        return 0
    elseif x>0
        return 1/(exp(x)-1)
    else
        throw(ArgumentError("x must be positive, you have written x=$x.")) 
    end
end

println("---------------------")

    # Simulation Parameters
# ===========================================
# Ntraj = 3 # number of trajectories
α0 = 2.0 # Starting level for the initial coherent state
χ = 0 # ħω/k_BT - if it's 0 it inputs the 0-temperature case, if positive it uses the Bose-Einstein distribution
    # Physical Parameters
# ===========================================
Nth = N_BE(χ) # initial mean thermal photon number in the cavity 
tlist = range(0,10,1000)
Ncut = Nth + 4*√(Nth*(1+Nth)) + abs(α0)^2+4*abs(α0) 
# |-> Dimension cutoff in the operators: μ+4σ of thermal state (steady state) + μ+4σ of initial coherent state
Ncut = Int(ceil(Ncut*1.5)) # We take +50% of the maximal population and round it to the upper value to be safe
println("⚠ Nth≈$(Int(ceil(Nth))) => Ncut=$Ncut")
ω_k = 1.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_k>>1)

# Operators
a = destroy(Ncut)

c_decay = √(1+Nth)*a
c_excite = √(Nth)*a'
c_ops = [c_decay, c_excite]

ψ0 = coherent(Ncut, α0)

n_op = a'*a
e_ops=[n_op]
p = plot(title=L"ħω=%$(χ)k_{B}T \qquad N_{cut}=%$(Ncut)", size=(1920,1080), legend=:topright, titlefont=20, xtickfont=11, ytickfont=11, legendfontsize=20);

H = ω_k*a'*a
sol = mesolve(H, ψ0, tlist, c_ops; e_ops=e_ops, progress_bar=Val(true))
Y = real.(sol.expect[1,:])
plot!(p, tlist, Y,label=L"ω/k=%$(ω_k)",lw=10, color=palette[1]);
n(t) = α0^2*exp(-t)+Nth*(1-exp(-t))
N = n.(tlist)
plot!(p, tlist, N, lw=3, label="model", color=palette[2]);

Y2 = abs.(N.-Y)
ytickslog = 10 .^round.(range(log10(minimum(Y2)),log10(maximum(Y2)),5), digits=2)
yticks = round.(range(log10(minimum(Y2)),log10(maximum(Y2)),5), digits=2)

# plot!(p, tlin, Ylin, color=palette[3], label="", yticks=ceil(minimum(Ylin)) ; inset=bbox(0.6, 0.4, 0.4, 0.4), subplot=2)
plot!(p, tlist, Y2, color=palette[3], label="", yticks=yticks ; inset=bbox(0.6, 0.4, 0.4, 0.4), subplot=2)
savefig(p, "test.svg");
