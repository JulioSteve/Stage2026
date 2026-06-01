ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings, Printf
theme(:dao)
palette = theme_palette(:dao)

function N_BE(x)
    if x == -1
        return 0 # χ→+∞ ≡ T→0 ⇒ N_th→0
    else
        return 1/(exp(x)-1)
    end
end

println("---------------------")

    # Simulation Parameters
# ===========================================
α0 = 2.0 # Starting level for the initial coherent state
χ = 5
Nth = N_BE(χ) # mean thermal photon number in the cavity 
tlist = range(0,10,1000)
ω_k = 1.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_k>>1)


Ncut = Nth + 4*√(Nth*(1+Nth)) + abs(α0)^2+4*abs(α0) # Dimension cutoff in the operators: μ+4σ of thermal state (steady state) + μ+4σ of initial coherent state

prct = 100 # POURCENTAGE DE TRONCATURE EN PLUS
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



println("---------------------")